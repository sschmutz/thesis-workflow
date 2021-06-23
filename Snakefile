rule get_classified_sequences:
    input:
        raw_reads="data/sequencing_files/{sample}.fastq.gz",
        unclassified_reads="data/sequencing_files/{sample}_unclassified-reads.fastq.gz"
    output:
        intermediate_reads=temp("data/sequencing_files/{sample}_intermediate.fastq"),
        good_reads=temp("data/sequencing_files/{sample}_good.fastq"),
        bad_reads=temp("data/sequencing_files/{sample}_bad.fastq"),
        good_reads_list=temp("data/sequencing_files/{sample}_good.lst"),
        unclassified_reads_list=temp("data/sequencing_files/{sample}_unclassified.lst"),
        good_classified_reads_list=temp("data/sequencing_files/{sample}_good-classified.lst"),
        classified_reads="data/sequencing_files/{sample}_classified-reads.fastq.gz"
    run:
        # trim and filter for sequences > 75nt
        shell("seqtk trimfq {input.raw_reads} | seqtk seq -L 75 - > {output.intermediate_reads}")

        # filter low quality sequences and sequences with low entropy
        shell("prinseq-lite.pl -fastq data/sequencing_files/{wildcards.sample}_intermediate.fastq -lc_method entropy -lc_threshold 70 -min_qual_mean 20 -out_good data/sequencing_files/{wildcards.sample}_good -out_bad data/sequencing_files/{wildcards.sample}_bad")

        # get sequence ids from all good and unclssified reads
        shell('grep "@M0" data/sequencing_files/{wildcards.sample}_good.fastq | sort | sed "s/@//" > {output.good_reads_list}')
        shell('zgrep "@M0" {input.unclassified_reads} | sort | sed "s/@//" > {output.unclassified_reads_list}')

        # get only classified sequences (difference between filtered and undetermined)
        shell("comm -23 data/sequencing_files/{wildcards.sample}_good.lst data/sequencing_files/{wildcards.sample}_unclassified.lst > {output.good_classified_reads_list}")
        shell("seqkit grep -n -f data/sequencing_files/{wildcards.sample}_good-classified.lst data/sequencing_files/{wildcards.sample}_good.fastq | gzip > {output.classified_reads}")


rule get_quality_measures:
    input:
        unclassified_reads="data/sequencing_files/{sample}_unclassified-reads.fastq.gz",
        classified_reads="data/sequencing_files/{sample}_classified-reads.fastq.gz"
    output:
        unclassified_reads_quality="data/quality_measures/{sample}_unclassified-reads.json",
        classified_reads_quality="data/quality_measures/{sample}_classified-reads.json"
    run:
        shell("fastp -i {input.unclassified_reads} --overrepresentation_analysis --low_complexity_filter -j {output.unclassified_reads_quality} -h /dev/null")
        shell("fastp -i {input.classified_reads} --overrepresentation_analysis --low_complexity_filter -j {output.classified_reads_quality} -h /dev/null")


rule get_sequence_labels:
    input:
        reference_human="data/virmet_dbs/GRCh38.fasta",
        cram_human="data/virmet_output/{sample}/good_humanGRCh38.cram",
        reference_bacterial_1="data/virmet_dbs/bact1.fasta",
        cram_bacterial_1="data/virmet_output/{sample}/good_humanGRCh38_bact1.cram",
        reference_bacterial_2="data/virmet_dbs/bact2.fasta",
        cram_bacterial_2="data/virmet_output/{sample}/good_humanGRCh38_bact1_bact2.cram",
        reference_bacterial_3="data/virmet_dbs/bact3.fasta",
        cram_bacterial_3="data/virmet_output/{sample}/good_humanGRCh38_bact1_bact2_bact3.cram",
        reference_bacterial_4="data/virmet_dbs/bact4.fasta",
        cram_bacterial_4="data/virmet_output/{sample}/good_humanGRCh38_bact1_bact2_bact3_bact4.cram",
        reference_bacterial_5="data/virmet_dbs/bact5.fasta",
        cram_bacterial_5="data/virmet_output/{sample}/good_humanGRCh38_bact1_bact2_bact3_bact4_bact5.cram",
        reference_fungal="data/virmet_dbs/fungi1.fasta",
        cram_fungal="data/virmet_output/{sample}/good_humanGRCh38_bact1_bact2_bact3_bact4_bact5_fungi1.cram",
        fastq_viral="data/virmet_output/{sample}/viral_reads.fastq.gz"
    output:
        list_human="data/classification/{sample}_human.lst.gz",
        list_bacterial_1="data/classification/{sample}_bacterial_1.lst.gz",
        list_bacterial_2="data/classification/{sample}_bacterial_2.lst.gz",
        list_bacterial_3="data/classification/{sample}_bacterial_3.lst.gz",
        list_bacterial_4="data/classification/{sample}_bacterial_4.lst.gz",
        list_bacterial_5="data/classification/{sample}_bacterial_5.lst.gz",
        list_fungal="data/classification/{sample}_fungal.lst.gz",
        list_viral="data/classification/{sample}_viral.lst.gz"
    run:
        # get human sequence ids
        shell('samtools fastq --reference {input.reference_human} {input.cram_human} | grep "@M0" | sed "s/@//" | gzip > {output.list_human}')

        # get bacterial sequence ids
        shell('samtools fastq --reference {input.reference_bacterial_1} {input.cram_bacterial_1} | zgrep "@M0" | sed "s/@//" | gzip > {output.list_bacterial_1}')
        shell('samtools fastq --reference {input.reference_bacterial_2} {input.cram_bacterial_2} | zgrep "@M0" | sed "s/@//" | gzip > {output.list_bacterial_2}')
        shell('samtools fastq --reference {input.reference_bacterial_3} {input.cram_bacterial_3} | zgrep "@M0" | sed "s/@//" | gzip > {output.list_bacterial_3}')
        shell('samtools fastq --reference {input.reference_bacterial_4} {input.cram_bacterial_4} | zgrep "@M0" | sed "s/@//" | gzip > {output.list_bacterial_4}')
        shell('samtools fastq --reference {input.reference_bacterial_5} {input.cram_bacterial_5} | zgrep "@M0" | sed "s/@//" | gzip > {output.list_bacterial_5}')

        # get fungal sequence ids
        shell('samtools fastq --reference {input.reference_fungal} {input.cram_fungal} | zgrep "@M0" | sed "s/@//" | gzip > {output.list_fungal}')

        # get viral sequence ids
        shell('zgrep "@M0" {input.fastq_viral} | sed "s/@//" | gzip > {output.list_viral}')


rule metagenome_assembly:
    input:
        classified_reads="data/sequencing_files/{sample}_classified-reads.fastq.gz",
        unclassified_reads="data/sequencing_files/{sample}_unclassified-reads.fastq.gz"
    output:
        assembly_folder=temp(directory("data/metagenome_assembly/{sample}"))
    shell:
        "megahit -r {input.classified_reads},{input.unclassified_reads} -m 0.5 -t 4 -o {output.assembly_folder}"


rule remove_intermediate_contigs:
    input:
        assembly_folder="data/metagenome_assembly/{sample}"
    output:
        assembly_stats=temp("data/metagenome_assembly/{sample}_final.contigs.lst"),
        assembly="data/metagenome_assembly/{sample}_final.contigs.fasta.gz"
    run:
        shell('grep ">" data/metagenome_assembly/{wildcards.sample}/final.contigs.fa | sed "s/>//" > data/metagenome_assembly/{wildcards.sample}_final.contigs.lst')
        shell("cp data/metagenome_assembly/{wildcards.sample}/final.contigs.fa data/metagenome_assembly/{wildcards.sample}_final.contigs.fasta")
        shell("gzip data/metagenome_assembly/{wildcards.sample}_final.contigs.fasta")


rule read_mapping:
    input:
        unclassified_reads="data/sequencing_files/{sample}_unclassified-reads.fastq.gz",
        classified_reads="data/sequencing_files/{sample}_classified-reads.fastq.gz",
        assembly_fasta="data/metagenome_assembly/{sample}_final.contigs.fasta.gz"
    output:
        "data/metagenome_assembly_read_mapping/{sample}_aln.tsv.gz"
    shell:
        "minimap2 -ax sr {input.assembly_fasta} <(cat {input.unclassified_reads} {input.classified_reads}) | samtools view | cut -f 1,3 | gzip > {output}"

rule label_contigs:
    input:
        list_human="data/classification/{sample}_human.lst.gz",
        list_bacterial_1="data/classification/{sample}_bacterial_1.lst.gz",
        list_bacterial_2="data/classification/{sample}_bacterial_2.lst.gz",
        list_bacterial_3="data/classification/{sample}_bacterial_3.lst.gz",
        list_bacterial_4="data/classification/{sample}_bacterial_4.lst.gz",
        list_bacterial_5="data/classification/{sample}_bacterial_5.lst.gz",
        list_fungal="data/classification/{sample}_fungal.lst.gz",
        list_viral="data/classification/{sample}_viral.lst.gz"
    params:
        # the threshold value is in percent
        threshold="{threshold}"
    output:
        "data/undetermined_class_label/{sample}_{threshold}.csv",
        temp("data/classification/{sample}_unclassified-in-contig_{threshold}.lst"),
        temp("data/classification/{sample}_unclassified-not-in-contig_{threshold}.lst")
    script:
        "scripts/label_contigs.R"


rule get_quality_measures_subcategory:
    input:
        unclassified_reads="data/sequencing_files/{sample}_unclassified-reads.fastq.gz",
        unclassified_in_contig_reads_list="data/classification/{sample}_unclassified-in-contig_{threshold}.lst",
        unclassified_not_in_contig_reads_list="data/classification/{sample}_unclassified-not-in-contig_{threshold}.lst"
    output:
        unclassified_in_contig_reads="data/sequencing_files/{sample}_unclassified-in-contig_{threshold}.fastq.gz",
        unclassified_in_contig_reads_quality="data/quality_measures/{sample}_unclassified-in-contig-reads_{threshold}.json",
        unclassified_not_in_contig_reads="data/sequencing_files/{sample}_unclassified-not-in-contig_{threshold}.fastq.gz",
        unclassified_not_in_contig_reads_quality="data/quality_measures/{sample}_unclassified-not-in-contig-reads_{threshold}.json"
    run:
        # extract sequences listed in unclassified_in_contig and unclassified_not_in_contig lists
        shell("seqkit grep -f {input.unclassified_in_contig_reads_list} {input.unclassified_reads} | gzip > {output.unclassified_in_contig_reads}")
        shell("seqkit grep -f {input.unclassified_not_in_contig_reads_list} {input.unclassified_reads} | gzip > {output.unclassified_not_in_contig_reads}")

        # run fastp
        shell("fastp -i {output.unclassified_in_contig_reads} --overrepresentation_analysis --low_complexity_filter -j {output.unclassified_in_contig_reads_quality} -h /dev/null")
        shell("fastp -i {output.unclassified_not_in_contig_reads} --overrepresentation_analysis --low_complexity_filter -j {output.unclassified_not_in_contig_reads_quality} -h /dev/null")
