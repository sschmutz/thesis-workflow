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
