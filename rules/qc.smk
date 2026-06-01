"""
QC Rules: FastQC (per-sample) + MultiQC (aggregate report)
"""

rule fastqc_raw:
    """Run FastQC on raw paired-end FASTQ files."""
    input:
        r1 = get_fastq_r1,
        r2 = get_fastq_r2,
    output:
        html_r1 = "{results}/qc/{sample}_R1_fastqc.html",
        html_r2 = "{results}/qc/{sample}_R2_fastqc.html",
        zip_r1  = "{results}/qc/{sample}_R1_fastqc.zip",
        zip_r2  = "{results}/qc/{sample}_R2_fastqc.zip",
    params:
        outdir = "{results}/qc",
    threads: 2
    resources:
        mem_mb = 4000,
    conda:
        "../envs/alignment.yaml"
    log:
        "{results}/logs/fastqc/{sample}.log"
    shell:
        """
        fastqc --threads {threads} \
               --outdir {params.outdir} \
               {input.r1} {input.r2} \
               2> {log}
        """


rule fastqc_trimmed:
    """Run FastQC on trimmed FASTQ files (post-Trimmomatic)."""
    input:
        r1 = "{results}/trimmed/{sample}_R1_trimmed.fastq.gz",
        r2 = "{results}/trimmed/{sample}_R2_trimmed.fastq.gz",
    output:
        html_r1 = "{results}/qc/{sample}_trimmed_R1_fastqc.html",
        html_r2 = "{results}/qc/{sample}_trimmed_R2_fastqc.html",
    params:
        outdir = "{results}/qc",
    threads: 2
    resources:
        mem_mb = 4000,
    conda:
        "../envs/alignment.yaml"
    log:
        "{results}/logs/fastqc_trimmed/{sample}.log"
    shell:
        """
        fastqc --threads {threads} \
               --outdir {params.outdir} \
               {input.r1} {input.r2} \
               2> {log}
        """


rule multiqc:
    """Aggregate all FastQC and alignment reports into a single MultiQC HTML."""
    input:
        # Collect all raw FastQC zip files
        expand("{results}/qc/{sample}_R1_fastqc.zip",
               results=RESULTS, sample=SAMPLE_IDS),
        expand("{results}/qc/{sample}_R2_fastqc.zip",
               results=RESULTS, sample=SAMPLE_IDS),
        # Include STAR alignment logs
        expand("{results}/alignments/{sample}.Log.final.out",
               results=RESULTS, sample=SAMPLE_IDS),
    output:
        html   = "{results}/qc/multiqc_report.html",
        data   = directory("{results}/qc/multiqc_data"),
    params:
        indir  = "{results}",
        outdir = "{results}/qc",
    conda:
        "../envs/alignment.yaml"
    log:
        "{results}/logs/multiqc.log"
    shell:
        """
        multiqc {params.indir} \
                --outdir {params.outdir} \
                --filename multiqc_report.html \
                --force \
                2> {log}
        """
