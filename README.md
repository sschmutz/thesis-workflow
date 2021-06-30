# Characterizing undetermined sequences from metagenomic sequencing
> Snakemake Workflow of Master's Thesis

## Setup
Get the contents of this repository by cloning it from GitHub by running  
`git clone https://github.com/sschmutz/thesis-workflow.git`.  

Edit the `Snakefile` by adding all names of the samples which should be processed and determine a threshold (in percent) which is applied to label the contigs in rule "infer_class_labels".

### Data Required
*Note: The data is not part of this repository as it contains sensitive information.*  

Prepare a folder within the cloned git repository named `data/sequencing_files/` which should containing the compressed raw- and unclassified sequencing files (`{sample}.fastq.gz` and `{sample}_unclassified-reads.fastq.gz`) for each sample which will be analysed. Creating a symbolic link instead of copying the files is also possible. Remaining the same folder and filename structure suggested here allows that the `Snakefile` doesn't need to be adapted beyond specifying the sample names and threshold value.  
In another folder named `data/virmet_dbs` the fasta files of all database sequences which were used for VirMet need to be present.  
The virmet output itself has also to be prepared in a separate folder `data/virmet_output` where the results of each sample needs to be present in a subfolder named after the sample. Required are all cram files created during the decontamination steps of VirMet and `viral_reads.fastq.gz` which contains all sequencing reads which were assigned to a virus.

### Software required
Listed below are all primary software tools (without dependencies) and their version which was used for the thesis. All could be installed (including required dependencies) via the package manager Conda.  
A complete list can be found in the [`environment.yml`](environment.yml) document.

Python v3.9.5, van Rossum (1995)  
Snakemake v6.4.1, Mölder et al. (2021)  
seqtk v1.3, Li (2021)  
prinseq v0.20.4, Schmieder and Edwards (2011)  
seqkit v0.16.1, Wei et al. (2016)  
megahit v1.2.9, Li et al. (2015)  
minimap2 v2.17, Li et al. (2018)  
samtools v1.12, Danecek et al. (2021)  
fastp v0.20.1, Chen et al. (2018)  
R v4.0.3, R Core Team (2020)  
R - tidyverse v1.3.1, Wickham et al. (2019)  
R - here v1.0.1, Müller (2020)  


## Rules
The workflow is divided in six different steps, also called rules in Snakemake. And a pseudo-rule (all) which combines all rules.  
Following figure shows a directed acyclic graph (DAG) which helps visualising the rules and how they are connected.  


![DAG of all rules](dag.svg)

### All
This pseudo-rule lists all files that should be created and combines all rules explained below.  
Running this rule executes all processes at once for the samples and threshold value specified at the top of the `Snakefile`. It is the simplest way to execute the whole workflow and can be run using the following command:  
```
snakemake --cores 4
```

The following output files are required to replicate the findings presented in the Results section of the thesis and are thus listed in this final rule:  
- **Quality profiling** for all samples and types (*classified-reads*, *unclassified-reads*, *unclassified-in-contig-reads* and *unclassified-not-in-contig-reads*) written to `data/quality_metrics/{sample}_{type}.json`.  

- **Assembly statistics** (contig name, flag, multi and length) for each sample is written to `data/metagenome_assembly/{sample}_final.contigs.lst`.  

- The **outcome of the Semi-supervised approach** to classify previously unclassified reads is written for each sample and chosen threshold to `data/undetermined_class_label/{sample}_{threshold}.csv`.

### Get quality filtered sequences
Since there is no one file containing all classified sequences among the VirMet output, it needs to be created first.  
There are however files containing all raw sequences and all unclassified sequences respectively (prepared in section [Data Required](#data-required)). It is therefore possible to get the classified sequences as follows:  

raw sequences - sequences removed by qc - unclassified sequences = classified sequences  

And since there's also no information about the reads which failed VirMets QC, the quality filtering steps need to be repeated again.  
Next to some intermediate temporary files (which are automatically deleted after finishing this rule) the only output is a compressed fastq file named `data/sequencing_files/{sample}_classified-reads.fastq.gz` which contains all classified sequencing reads for each sample.  

To therefore run this rule for a sample, the following command can be used (replace `{sample}` with an actual sample name):  
```
snakemake --cores 4 data/sequencing_files/{sample}_classified-reads.fastq.gz
```

### Metagenome assembly
All reads per sample which passed the VirMet quality filtering steps (`data/sequencing_files/{sample}_classified-reads.fastq.gz` and `data/sequencing_files/{sample}_unclassified-reads.fastq.gz`) are used for a metagenome assembly using megahit with the default parameter. Only the memory (`-m`) and thread (`-t`) usage are specified.  
From the metagenome assembly, only the minimum which is needed to proceed is kept and compressed to optimize storage space required. These are the final.contigs (`data/metagenome_assembly/{sample}_final.contigs.fasta.gz`) and the assembly stats (`data/metagenome_assembly/{sample}_final.contigs.lst`) which are the sequencing headers.  
If the intermediate contigs should be kept, remove the `temp()` tag of the metagenome assembly output.  

To run this rule for a sample, the following command can be used (replace {sample} with an actual sample name):  
```
snakemake --cores 4 data/metagenome_assembly/{sample}
```

### Read mapping
Since during the metagenome assembly the information of which read is part of which contig is lost, all sequencing reads used for the assembly have to be mapped back to those contigs.  
The read mapping is done using minimap2 with the settings for short accurate genomic reads (`-ax sr`). Again, to save storage space, only the minimum of information required (query sequence name, target sequence name) is stored in a table (`data/metagenome_assembly_read_mapping/{sample}_aln.tsv.gz`). For the mapped reads, the target sequence name is the contig name it mapped to, while for the unmapped reads it is an asterisk.  

To run this rule for a sample, the following command can be used (replace {sample} with an actual sample name):  
```
snakemake --cores 4 data/metagenome_assembly_read_mapping/{sample}_aln.tsv.gz
```

### Get VirMet labels
Using VirMet, the sequencing reads were classified into different classes. The main ones are: human, bacterial, fungal and viral.  
For the Semi-supervised approach, the labels of all classified reads are required. With the exception of the viral results, VirMet stores the alignment information in CRAM (compressed SAM) files.  
To access the information of which read was mapped against a read of which reference database, one needs to uncompress the CRAM files. This can be done using samtools. It requires the CRAM files itself and the reference database (fasta file) which was used for the alignment by VirMet. All mapped reads can then be written to a fastq file and by extracting just the sequence ids, one gets a list of all reads which were aligned to a sequence of each of the classes.  
It is a bit more straight forward for the viral reads, because those are stored in `viral_reads.fastq.gz` for each sample. Only the second part of what is described in the previous section, extracting the sequence ids, needs to be done to recover the viral ids.  

The output are compressed lists containing all sequence ids assigned by VirMet to the different classes.  

To run this rule for a sample, the following command can be used (replace {sample} with an actual sample name):  
```
snakemake --cores 4 data/classification/{sample}_viral.lst.gz
```

### Infer class labels
Given a threshold value in percent as described in section [Setup](#setup), the previously unclassified reads are put into classes.  
This happens using the R script `infer_class_labels.R` which uses the VirMet class labels and the contig information of each sequencing read (previous two rules) to infer class labels of previously unclassified reads.  
The summary statistics are written to the `undetermined_class_label` folder and a list of sequencing reads for the two classes "unclassified_in_contig" and "unclassified_not_in_contig" are written to the `classification` folder.  

To run this rule for a sample, the following command can be used (replace {sample} with an actual sample name and define a {threshold}):  
```
snakemake --cores 4 data/undetermined_class_label/{sample}_{threshold}.csv
```

### Get quality metrics
The bioinformatic tool fastp is used to get multiple quality metrics of the "classified" and "unclassified" reads separate for each sample. The "unclassified" fraction is additionally further divided into "unclassified-in-contig" an "unclassified-not-in-contig" depending if they could be mapped to a longer contig or not. The output is written to json files in the `quality_measures/` folder.  

To run this rule for a sample, the following command can be used (replace {sample} with an actual sample name):  
```
snakemake --cores 4 data/quality_metrics/{sample}_classified-reads.json`
```

## Snakemake usage and documentation
There are different ways to execute parts or the whole snakemake workflow (see [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/)).  
