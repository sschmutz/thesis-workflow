# This function should infer the class labels of the contigs created through
# metagenome assembly
# Its output is the basis to create an overview on how the undetermined reads
# are distributed

# inputs are:
# threshold (in percent) for class label inference (at least x% of the reads which map to
# a contig have to be labeled by VirMet) -> calculate fractions, remove unclassified
# and pick the label with the highest fraction

library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(here)


read_virmet_classes <- function(sample_name_short, virmet_class){
  class_label_path <- here("classification", paste0(sample_name_short, "_", virmet_class, ".lst.gz"))
  
  # only returns first part of class label ("bacterial_1" gives "bacterial")
  virmet_class_short <- str_split(virmet_class, pattern = "_")[[1]][1]
  
  virmet_class_labels <-
    read_csv(class_label_path, col_names = c("read_name")) %>%
    mutate(virmet_class = virmet_class_short)
  
  return(virmet_class_labels)
}

label_contigs <- function(sample_name_short, threshold){
  
  alignment_path <- here("metagenome_assembly_read_mapping", paste0(sample_name_short, "_aln.tsv.gz"))
  
  # caution, some reads might have mapped to more than one contig and are listed multiple times
  alignment <-
    read_delim(alignment_path, delim = "\t", col_names = c("read_name", "contig_name"), na = "*") %>%
    # remove multiple rows when a single reads was mapped multiple times to the same contig
    unique() %>%
    separate(contig_name, into = c(NA, "contig_number"), convert = TRUE)
  

  virmet_classes <- c("human",
                      "bacterial_1", "bacterial_2", "bacterial_3", "bacterial_4", "bacterial_5",
                      "fungal",
                      "viral")
  
  virmet_class_labels <-
    map_dfr(virmet_classes, ~read_virmet_classes(sample_name_short, .x)) %>%
    # viral read names contain an additional part, remove that
    separate(read_name, into = "read_name", sep = " ", extra = "drop")
  
  alignment_undetermined_reads <-
    anti_join(alignment, virmet_class_labels, by = "read_name")
  
  combined <-
    full_join(alignment, virmet_class_labels, by = "read_name") %>%
    mutate(virmet_class = replace_na(virmet_class, "undetermined"))
  
  contig_labels <-
    combined %>%
    filter(!is.na(contig_number)) %>%
    group_by(contig_number, virmet_class) %>%
    count() %>%
    group_by(contig_number) %>%
    mutate(n_total = sum(n)) %>%
    ungroup() %>%
    mutate(n_fraction = n/n_total) %>%
    # select top label (undetermined excluded) which is above threshold
    # to label the contigs
    filter(virmet_class != "undetermined",
           n_fraction > (threshold/100)) %>%
    group_by(contig_number) %>%
    slice_max(order_by = n_fraction, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  unclassified_in_contig <-
    alignment_undetermined_reads %>%
    filter(!is.na(contig_number)) %>%
    left_join(contig_labels, by = "contig_number")
  
  n_unclassified_not_in_contig <-
    alignment_undetermined_reads %>%
    filter(is.na(contig_number)) %>%
    nrow()
  
  classes <- tibble(class = c("human", "bacterial", "fungal", "viral", "unclassified_in_contig", "unclassified_not_in_contig"))
    
  class_order <- classes$class
  
  n_summary <-
    unclassified_in_contig %>%
    count(virmet_class) %>%
    rename(class = virmet_class) %>%
    mutate(class = replace_na(class, "unclassified_in_contig")) %>%
    add_row(class = "unclassified_not_in_contig", n = n_unclassified_not_in_contig) %>%
    full_join(classes, by = "class") %>%
    mutate(n = replace_na(n, 0),
           class = factor(class, levels = class_order),
           sample = sample_name_short,
           treshold_percent = threshold) %>%
    arrange(class)

  
  return(n_summary)
}

write_undetermined_class_label <- function(sample_name_short, threshold){
  
  output_path <- here("undetermined_class_label", paste0(sample_name_short, "_", threshold, ".csv"))
  
  label_contigs(sample_name_short, threshold) %>%
    write_csv(output_path)
  
}

list_human <- snakemake@input[["list_human"]]

# extract the sample name by looking at everything which is inbetween "/" and "_"
sample_name_short <- str_extract(list_human, pattern = "(?<=/)(.*)(?=_)")

threshold <- as.numeric(snakemake@params[["threshold"]])

write_undetermined_class_label(sample_name_short, threshold)





