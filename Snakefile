"""
RNA-seq Differential Expression Pipeline
Master Snakefile

Description: End-to-end pipeline from raw FASTQ to DEGs and pathway enrichment.
             Manages QC, trimming, alignment, quantification, and statistical analysis.
"""

import pandas as pd
from pathlib import Path

localrules: all

wildcard_constraints:
    sample = r"[A-Za-z0-9_\-]+"

# ─── Configuration ─────────────────────────────────────────────────────────────

configfile: "config/config.yaml"

SAMPLES = pd.read_csv(config["samples"], sep="\t", index_col="sample_id")
SAMPLE_IDS = list(SAMPLES.index)
RESULTS = config["results_dir"]


# ─── Helper functions ──────────────────────────────────────────────────────────

def get_fastq_r1(wildcards):
    return SAMPLES.loc[wildcards.sample, "fastq_r1"]

def get_fastq_r2(wildcards):
    return SAMPLES.loc[wildcards.sample, "fastq_r2"]


# ─── Target rule ───────────────────────────────────────────────────────────────

rule all:
    input:
        # QC outputs
        expand("{results}/qc/{sample}_R1_fastqc.html", results=RESULTS, sample=SAMPLE_IDS),
        expand("{results}/qc/{sample}_R2_fastqc.html", results=RESULTS, sample=SAMPLE_IDS),
        f"{RESULTS}/qc/multiqc_report.html",
        # Alignment outputs
        expand("{results}/alignments/{sample}.Aligned.sortedByCoord.out.bam",
               results=RESULTS, sample=SAMPLE_IDS),
        # Count matrix
        f"{RESULTS}/tables/raw_counts.tsv",
        # DESeq2 outputs
        f"{RESULTS}/tables/deseq2_results.tsv",
        f"{RESULTS}/tables/normalized_counts.tsv",
        # Figures
        f"{RESULTS}/figures/volcano_plot.png",
        f"{RESULTS}/figures/pca_plot.png",
        f"{RESULTS}/figures/heatmap_top50.png",
        f"{RESULTS}/figures/ma_plot.png",
        f"{RESULTS}/figures/gsea_dotplot.png",


# ─── Include rule modules ──────────────────────────────────────────────────────

include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/quantify.smk"
include: "rules/deseq2.smk"
