"""
Alignment Rules: STAR genome index build + read alignment.
"""

rule star_genome_index:
    """
    Build STAR genome index from reference FASTA and GTF.
    Only needs to run once. Output directory is checked as sentinel.
    Runtime: ~45–60 min with 8 threads, 32 GB RAM.
    """
    input:
        fasta = config["genome_fasta"],
        gtf   = config["gtf_file"],
    output:
        directory(config["genome_dir"]),
    params:
        genome_dir   = config["genome_dir"],
        overhang     = 99,       # ReadLength - 1; adjust if reads differ from 100bp
    threads: 8
    resources:
        mem_mb = 40000,
    conda:
        "../envs/alignment.yaml"
    log:
        f"{RESULTS}/logs/star_index/star_genome_index.log"
    shell:
        """
        mkdir -p {params.genome_dir}
        STAR --runMode genomeGenerate \
             --runThreadN {threads} \
             --genomeDir {params.genome_dir} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.overhang} \
             2> {log}
        """


rule star_align:
    """
    Align trimmed paired-end reads to the reference genome with STAR.
    Outputs coordinate-sorted BAM files ready for downstream quantification.
    """
    input:
        r1         = "{results}/trimmed/{sample}_R1_trimmed.fastq.gz",
        r2         = "{results}/trimmed/{sample}_R2_trimmed.fastq.gz",
        genome_dir = config["genome_dir"],
    output:
        bam        = "{results}/alignments/{sample}.Aligned.sortedByCoord.out.bam",
        log_final  = "{results}/alignments/{sample}.Log.final.out",
        sj_tab     = "{results}/alignments/{sample}.SJ.out.tab",
    params:
        prefix         = "{results}/alignments/{sample}.",
        genome_dir     = config["genome_dir"],
        genome_load    = config["star"]["genome_load"],
        out_sam_type   = config["star"]["out_sam_type"],
        out_sam_attrs  = config["star"]["out_sam_attrs"],
    threads: 8
    resources:
        mem_mb = 36000,
    conda:
        "../envs/alignment.yaml"
    log:
        "{results}/logs/star_align/{sample}.log"
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {params.genome_dir} \
             --genomeLoad {params.genome_load} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outSAMtype {params.out_sam_type} \
             --outSAMattributes {params.out_sam_attrs} \
             --outFileNamePrefix {params.prefix} \
             --outSAMunmapped Within \
             --quantMode GeneCounts \
             2> {log}
        """


rule samtools_index:
    """Index BAM file for downstream tools."""
    input:
        bam = "{results}/alignments/{sample}.Aligned.sortedByCoord.out.bam",
    output:
        bai = "{results}/alignments/{sample}.Aligned.sortedByCoord.out.bam.bai",
    conda:
        "../envs/alignment.yaml"
    log:
        "{results}/logs/samtools_index/{sample}.log"
    shell:
        "samtools index {input.bam} 2> {log}"


rule samtools_flagstat:
    """Generate alignment statistics per sample."""
    input:
        bam = "{results}/alignments/{sample}.Aligned.sortedByCoord.out.bam",
    output:
        stats = "{results}/qc/{sample}_flagstat.txt",
    conda:
        "../envs/alignment.yaml"
    log:
        "{results}/logs/flagstat/{sample}.log"
    shell:
        "samtools flagstat {input.bam} > {output.stats} 2> {log}"
