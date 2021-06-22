rule get_classified_sequences:
    input:
        raw_reads="sequencing_files/{sample}.fastq.gz",
        unclassified_reads="sequencing_files/{sample}_unclassified-reads.fastq.gz"
    output:
        intermediate_reads=temp("sequencing_files/{sample}_intermediate.fastq"),
        good_reads=temp("sequencing_files/{sample}_good.fastq"),
        bad_reads=temp("sequencing_files/{sample}_bad.fastq"),
        good_reads_list=temp("sequencing_files/{sample}_good.lst"),
        unclassified_reads_list=temp("sequencing_files/{sample}_unclassified.lst"),
        good_classified_reads_list=temp("sequencing_files/{sample}_good-classified.lst"),
        classified_reads="sequencing_files/{sample}_classified-reads.fastq.gz"
    run:
        # trim and filter for sequences > 75nt
        shell("seqtk trimfq {input.raw_reads} | seqtk seq -L 75 - > {output.intermediate_reads}")

        # filter low quality sequences and sequences with low entropy
        shell("prinseq-lite.pl -fastq sequencing_files/{wildcards.sample}_intermediate.fastq -lc_method entropy -lc_threshold 70 -min_qual_mean 20 -out_good sequencing_files/{wildcards.sample}_good -out_bad sequencing_files/{wildcards.sample}_bad")

        # get sequence ids from all good and unclssified reads
        shell('grep "@M0" sequencing_files/{wildcards.sample}_good.fastq | sort | sed "s/@//" > {output.good_reads_list}')
        shell('zgrep "@M0" {input.unclassified_reads} | sort | sed "s/@//" > {output.unclassified_reads_list}')

        # get only classified sequences (difference between filtered and undetermined)
        shell("comm -23 sequencing_files/{wildcards.sample}_good.lst sequencing_files/{wildcards.sample}_unclassified.lst > {output.good_classified_reads_list}")
        shell("seqkit grep -n -f sequencing_files/{wildcards.sample}_good-classified.lst sequencing_files/{wildcards.sample}_good.fastq | gzip > {output.classified_reads}")


rule get_quality_measures:
    input:
        unclassified_reads="sequencing_files/{sample}_unclassified-reads.fastq.gz",
        classified_reads="sequencing_files/{sample}_classified-reads.fastq.gz"
    output:
        unclassified_reads_quality="quality_measures/{sample}_unclassified-reads.json",
        unclassified_reads_quality_report="quality_measures/{sample}_unclassified-reads.html",
        classified_reads_quality="quality_measures/{sample}_classified-reads.json",
        classified_reads_quality_report="quality_measures/{sample}_classified-reads.html",
    run:
        shell("fastp -i {input.unclassified_reads} --overrepresentation_analysis --low_complexity_filter -j {output.unclassified_reads_quality} -h {output.unclassified_reads_quality_report}")
        shell("fastp -i {input.classified_reads} --overrepresentation_analysis --low_complexity_filter -j {output.classified_reads_quality} -h {output.classified_reads_quality_report}")


rule get_sequence_labels:
    input:
        reference_human="virmet_dbs/GRCh38.fasta",
        cram_human="virmet_output/{sample}/good_humanGRCh38.cram",
        reference_bacterial_1="virmet_dbs/bact1.fasta",
        cram_bacterial_1="virmet_output/{sample}/good_humanGRCh38_bact1.cram",
        reference_bacterial_2="virmet_dbs/bact2.fasta",
        cram_bacterial_2="virmet_output/{sample}/good_humanGRCh38_bact1_bact2.cram",
        reference_bacterial_3="virmet_dbs/bact3.fasta",
        cram_bacterial_3="virmet_output/{sample}/good_humanGRCh38_bact1_bact2_bact3.cram",
        reference_bacterial_4="virmet_dbs/bact4.fasta",
        cram_bacterial_4="virmet_output/{sample}/good_humanGRCh38_bact1_bact2_bact3_bact4.cram",
        reference_bacterial_5="virmet_dbs/bact5.fasta",
        cram_bacterial_5="virmet_output/{sample}/good_humanGRCh38_bact1_bact2_bact3_bact4_bact5.cram",
        reference_fungal="virmet_dbs/fungi1.fasta",
        cram_fungal="virmet_output/{sample}/good_humanGRCh38_bact1_bact2_bact3_bact4_bact5_fungi1.cram",
        fastq_viral="virmet_output/{sample}/viral_reads.fastq.gz"
    output:
        fastq_human=temp("classification/{sample}_human.fastq"),
        list_human="classification/{sample}_human.lst",
        fastq_bacterial_1=temp("classification/{sample}_bacterial_1.fastq"),
        list_bacterial_1="classification/{sample}_bacterial_1.lst",
        fastq_bacterial_2=temp("classification/{sample}_bacterial_2.fastq"),
        list_bacterial_2="classification/{sample}_bacterial_2.lst",
        fastq_bacterial_3=temp("classification/{sample}_bacterial_3.fastq"),
        list_bacterial_3="classification/{sample}_bacterial_3.lst",
        fastq_bacterial_4=temp("classification/{sample}_bacterial_4.fastq"),
        list_bacterial_4="classification/{sample}_bacterial_4.lst",
        fastq_bacterial_5=temp("classification/{sample}_bacterial_5.fastq"),
        list_bacterial_5="classification/{sample}_bacterial_5.lst",
        fastq_fungal=temp("classification/{sample}_fungal.fastq"),
        list_fungal="classification/{sample}_fungal.lst",
        list_viral="classification/{sample}_viral.lst"
    run:
        # get human sequence ids
        shell("samtools fastq --reference {input.reference_human} {input.cram_human} | gzip > {output.fastq_human}")
        shell('zgrep "@M0" {output.fastq_human} | sed "s/@//" > {output.list_human}')

        # get bacterial sequence ids
        shell("samtools fastq --reference {input.reference_bacterial_1} {input.cram_bacterial_1} | gzip > {output.fastq_bacterial_1}")
        shell('zgrep "@M0" {output.fastq_bacterial_1} | sed "s/@//" > {output.list_bacterial_1}')
        shell("samtools fastq --reference {input.reference_bacterial_2} {input.cram_bacterial_2} | gzip > {output.fastq_bacterial_2}")
        shell('zgrep "@M0" {output.fastq_bacterial_2} | sed "s/@//" > {output.list_bacterial_2}')
        shell("samtools fastq --reference {input.reference_bacterial_3} {input.cram_bacterial_3} | gzip > {output.fastq_bacterial_3}")
        shell('zgrep "@M0" {output.fastq_bacterial_3} | sed "s/@//" > {output.list_bacterial_3}')
        shell("samtools fastq --reference {input.reference_bacterial_4} {input.cram_bacterial_4} | gzip > {output.fastq_bacterial_4}")
        shell('zgrep "@M0" {output.fastq_bacterial_4} | sed "s/@//" > {output.list_bacterial_4}')
        shell("samtools fastq --reference {input.reference_bacterial_5} {input.cram_bacterial_5} | gzip > {output.fastq_bacterial_5}")
        shell('zgrep "@M0" {output.fastq_bacterial_5} | sed "s/@//" > {output.list_bacterial_5}')

        # get fungal sequence ids
        shell("samtools fastq --reference {input.reference_fungal} {input.cram_fungal} | gzip > {output.fastq_fungal}")
        shell('zgrep "@M0" {output.fastq_fungal} | sed "s/@//" > {output.list_fungal}')

        # get viral sequence ids
        shell('zgrep "@M0" {input.fastq_viral} | sed "s/@//" > {output.list_viral}')


rule metagenome_assembly:
    input:
        classified_reads="sequencing_files/{sample}_classified-reads.fastq.gz",
        unclassified_reads="sequencing_files/{sample}_unclassified-reads.fastq.gz"
    output:
        "metagenome_assembly/{sample}/",
        assembly_fasta="metagenome_assembly/{sample}/final.contigs.fa",
        assembly_fastg="metagenome_assembly/{sample}/final.contigs.fastg",
        assembly_list="metagenome_assembly/{sample}/final.contigs.lst"
    run:
        shell("megahit -r {input.classified_reads},{input.unclassified_reads} -m 0.5 -t 4 -o metagenome_assembly/{wildcards.sample}")
        shell("megahit_core contig2fastg 141 {output.assembly_fasta} > {output.assembly_fastg}")
        shell('grep ">" {output.assembly_fasta} | sed "s/>//" > {output.assembly_list}')
