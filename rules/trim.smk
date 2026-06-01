"""
Trim Rule: Adapter trimming and quality filtering with Trimmomatic.
"""

rule trimmomatic:
    """
    Trim adapters and low-quality bases from paired-end reads.
    Parameters are set in config/config.yaml under trimmomatic:
    """
    input:
        r1 = get_fastq_r1,
        r2 = get_fastq_r2,
    output:
        r1_paired   = "{results}/trimmed/{sample}_R1_trimmed.fastq.gz",
        r2_paired   = "{results}/trimmed/{sample}_R2_trimmed.fastq.gz",
        r1_unpaired = "{results}/trimmed/{sample}_R1_unpaired.fastq.gz",
        r2_unpaired = "{results}/trimmed/{sample}_R2_unpaired.fastq.gz",
    params:
        adapters       = config["trimmomatic"]["adapters"],
        leading        = config["trimmomatic"]["leading"],
        trailing       = config["trimmomatic"]["trailing"],
        slidingwindow  = config["trimmomatic"]["slidingwindow"],
        minlen         = config["trimmomatic"]["minlen"],
    threads: 4
    resources:
        mem_mb = 8000,
    conda:
        "../envs/alignment.yaml"
    log:
        "{results}/logs/trimmomatic/{sample}.log"
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            {input.r1} {input.r2} \
            {output.r1_paired} {output.r1_unpaired} \
            {output.r2_paired} {output.r2_unpaired} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:{params.leading} \
            TRAILING:{params.trailing} \
            SLIDINGWINDOW:{params.slidingwindow} \
            MINLEN:{params.minlen} \
            2> {log}
        """
