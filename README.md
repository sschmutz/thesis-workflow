# Characterizing undetermined sequences from metagenomic sequencing
> Snakemake Workflow of Master's Thesis

## Setup

### Sequencing Data
*Note: The sequencing data is not in this repository as it contains sensitive data.*  

Prepare a folder containing the compressed raw- and unclassified sequencing files (`*.fastq.gz`) for each sample which should be analysed. Creating a symbolic link instead of copying the files is also possible.  
In another folder named `virmet_dbs` the fasta files of all database sequences which were used for VirMet need to be present.  
The virmet output itself has also to be prepared in a separate folder `virmet_output` where the results of each sample is copied in a subfolder named after the sample.

## Rules
The workflow is divided in different steps, also called rules in Snakemake.

To visualise the rules, a directed acyclic graph (DAG) can be written with the following command:

`snakemake --rulegraph metagenome_assembly/1000580287-AR-RNA quality_measures/1000580287-AR-RNA_classified-reads.json classification/1000580287-AR-RNA_viral.lst | dot -Tsvg > dag.svg`

![DAG of all rules](dag.svg)

### Get classified sequences
Since there is no one file containing all classified sequences of VirMet, it needs to be created first.  
There are however files containing all raw sequences and all unclassified sequences (prepared in section [Sequencing Data](sequencing-data)). It is therefore possible to get the classified sequences as follows:  

raw sequences - sequences removed by qc - unclassified sequences = classified sequences  

Next to some intermediate temporary files (which are deleted after finishing this rule) the only output is a compressed fastq file with the ending `*_classified-reads.fastq.gz` which contains all classified sequencing reads of each sample.


`snakemake --cores 8 data/sequencing_files/1000580287-AR-RNA_classified-reads.fastq.gz`

### Get quality measures
The tool fastp is used to get multiple quality measures of the unclassified and classified reads separate for each sample. The output is written to a json and html file in the `quality_measures/` folder.

`snakemake --cores 8 data/quality_measures/1000580287-AR-RNA_classified-reads.json`

### Get sequence labels
Using VirMet, the sequencing reads were classified into different classes. The main ones are: human, bacterial, fungal and viral.  
For the semisupervised approach, the labels of all reads are required. With the exception of the viral results, VirMet stores the alignment information in CRAM (compressed SAM) files.  
To access the information of which read was mapped against a read of which reference database, one needs to uncompress the CRAM files. This can be done using cramtools. It requires the CRAM files itself and the reference database (fasta file) which was used for the alignment. All mapped reads can then be written to a fastq file and by extracting just the sequence ids, one gets a list of all reads which were aligned to a sequence of each of the classes.  
It is a bit more straight forward for the viral reads, because those are saved in a `viral_reads.fastq.gz`. Only the second part of what is described in the previous section, extracting the sequence ids, needs to be done for those.

`snakemake --cores 8 data/classification/1000580287-AR-RNA_viral.lst.gz`

### Metagenome assembly
All reads per sample which passed the VirMet quality filtering steps are used for a metagenome assembly using megahit.

`snakemake --cores 8 data/metagenome_assembly/1000580287-AR-RNA`

### Remove intermediate contigs
From the metagenome assembly, only the minimum which is needed to proceed (final.contigs) are kept and compressed to optimize space required. If the intermediate contigs should be kept, remove the `temp()` tag of the metagenome assembly output.

`snakemake --cores 8 data/metagenome_assembly/1000580287-AR-RNA_final.contigs.fasta.gz`

### Read mapping
To get the information of which read is part of which contig, the sequencing reads have to be mapped back to those contigs.

`snakemake --cores 8 data/metagenome_assembly_read_mapping/1000580287-AR-RNA_aln.tsv.gz`

### Label contigs
Given a threshold value in percent, the previously unclassified reads are put into classes. The summary statistics are written to the `undetermined_class_label` folder and a list of sequencing reads for the two classes "unclassified_in_contig" and "unclassified_not_in_contig" are written to the `classification` folder.

`snakemake --cores 8 data/undetermined_class_label/1000580287-AR-RNA_5.csv`

### Get quality measures of subcategories
The previously used category "unclassified" can now be further divided into "unclassified-in-contig" and "unclassified-not-in-contig".
Those sequencing read files are created with this rule and a json quality report is created using fastp as previously done.

`snakemake --cores 8 data/quality_measures/1000580287-AR-RNA_unclassified-not-in-contig-reads_5.json`

## Usage example
There are different ways to execute parts or the whole snakemake workflow (see [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/)).  
Described here is an example to run the full workflow of one sample.

`snakemake --cores 8 metagenome_assembly_read_mapping/{1000576042-AR-DNA,1000576042-AR-RNA,1000580287-AR-DNA,1000580287-AR-RNA}_aln.tsv.gz`
