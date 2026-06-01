# RNA-seq Differential Expression Pipeline

**End-to-end, reproducible workflow from raw sequencing reads to differential expression and pathway analysis.**

This repository is a **portfolio-ready Snakemake template** for RNA-seq differential expression. It wires together industry-standard tools (FastQC, Trimmomatic, STAR, featureCounts, DESeq2, clusterProfiler) in a modular, config-driven workflow you can run locally or scale on an HPC cluster.

The default manifest includes **four example samples** (2 tumor, 2 normal). Add rows to `config/samples.tsv` to scale to larger cohorts without changing the rule logic.

**Intended application:** TCGA-BRCA breast cancer vs. normal tissue ([GEO: GSE62944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944)) on GRCh38 / GENCODE v44. You supply FASTQ and reference files; the repository ships workflow code, not sequencing data.

[![Snakemake](https://img.shields.io/badge/Snakemake-≥7.0-brightgreen)](https://snakemake.readthedocs.io)
[![Python](https://img.shields.io/badge/Python-3.10-blue)](https://python.org)
[![R](https://img.shields.io/badge/R-≥4.2-blue)](https://r-project.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow)](LICENSE)

---

## Validation Status

| Check | Status |
|---|---|
| Python scripts (`scripts/*.py`) | **Validated** — syntax checked with `ast.parse` |
| R scripts (`scripts/*.R`) | **Not run in CI on author machine** — review before first R execution; run `Rscript -e "parse(file='scripts/deseq2_analysis.R')"` locally |
| Snakemake dry-run | **Requires conda on Linux, macOS, or WSL** — use `bash setup.sh`, then `snakemake --dry-run --cores 1 --use-conda`. Native Windows without WSL is not supported for bioconda tool chains. |
| Full pipeline execution | **Not included in this repo** — requires user-downloaded paired-end FASTQ files and a GRCh38 reference (FASTA + GTF + STAR index). Outputs are written to `results/` (git-ignored). |
| GitHub Actions | Optional workflow in `.github/workflows/snakemake_dry_run.yml` runs a dry-run on Ubuntu with dummy FASTQ paths when you push to GitHub. |

This project demonstrates **workflow engineering and reproducibility**, not pre-computed biological results. After you run the pipeline, summarize DEG counts with `notebooks/exploratory_analysis.ipynb` or `results/tables/deseq2_results.tsv`.

---

## Portfolio Relevance

| Role focus | What this project shows |
|---|---|
| **Bioinformatics analyst** | End-to-end NGS pipeline design: QC → trim → align → quantify → statistics |
| **Computational biologist** | Splice-aware RNA-seq alignment, gene-level counting, DESeq2 contrasts, GSEA interpretation |
| **Data science (life sciences)** | Reproducible workflows, conda environments, config-driven parameters, R/Python visualization |

**Skills demonstrated:** Snakemake · conda environment management · STAR · featureCounts · DESeq2 (apeglm LFC shrinkage) · clusterProfiler GSEA · ggplot2 / pheatmap · scientific documentation · git hygiene (`.gitignore` for large omics data)

---

## Overview

### What the pipeline does

1. **QC** — FastQC per sample; MultiQC aggregate report  
2. **Trim** — Trimmomatic adapter and quality trimming  
3. **Align** — STAR to hg38 (GENCODE annotation)  
4. **Quantify** — featureCounts gene-level count matrix  
5. **Differential expression** — DESeq2 (volcano, PCA, heatmap, MA plots)  
6. **Pathway analysis** — GSEA with MSigDB Hallmark sets (clusterProfiler)

### Scaling beyond four samples

- Add one row per sample to `config/samples.tsv` (same columns: `sample_id`, `condition`, `fastq_r1`, `fastq_r2`, optional `batch`).  
- Snakemake expands rules automatically from the manifest.  
- For large cohorts, use more cores or a cluster profile (SLURM example in [Quick Start](#quick-start)).

---

## Repository Structure

```
rna-seq-pipeline/
├── Snakefile                      # Master workflow
├── setup.sh                       # One-time Snakemake conda environment
├── config/
│   ├── config.yaml                # Paths, thresholds, tool parameters
│   └── samples.tsv                # Sample manifest (4 examples by default)
├── rules/
│   ├── qc.smk                     # FastQC + MultiQC
│   ├── trim.smk                   # Trimmomatic
│   ├── align.smk                  # STAR index + alignment + SAMtools
│   ├── quantify.smk               # featureCounts
│   └── deseq2.smk                 # DESeq2 + GSEA
├── scripts/
│   ├── deseq2_analysis.R
│   ├── pathway_enrichment.R
│   └── clean_count_matrix.py
├── envs/                          # Per-stage conda environments
├── notebooks/
│   └── exploratory_analysis.ipynb # Post-run EDA (after pipeline completes)
├── .gitignore
├── LICENSE
└── README.md
```

Generated at run time (not committed): `results/`, `data/`, `.snakemake/`.

---

## Quick Start

### 1. Clone and create the Snakemake environment

Use Linux, macOS, or WSL with [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Mamba installed.

```bash
cd rna-seq-pipeline   # or clone from your GitHub remote, then enter the repo
bash setup.sh
conda activate snakemake
```

### 2. Add data and references

See [Dataset Instructions](#dataset-instructions). You need:

- Paired-end FASTQ files under `data/raw/` matching paths in `config/samples.tsv`
- `data/genome/hg38.fa` and `data/genome/gencode.v44.annotation.gtf`
- STAR index (built by the pipeline’s `star_genome_index` rule on first run)

### 3. Validate the workflow (dry run)

```bash
snakemake --dry-run --cores 1 --use-conda
```

A successful dry run lists jobs without executing them. Fix missing paths or conda issues before a full run.

### 4. Run the pipeline

```bash
snakemake --cores 8 --use-conda
```

**HPC (SLURM) example:**

```bash
snakemake --cluster "sbatch --mem={resources.mem_mb}M --cpus-per-task={threads}" \
          --jobs 50 --use-conda
```

### 5. Explore results

```bash
jupyter notebook notebooks/exploratory_analysis.ipynb
```

---

## Dataset Instructions

Data are **not** bundled in this repository (see `.gitignore`).

### Recommended first run: four samples

The default `config/samples.tsv` defines two tumor and two normal samples. Download four paired FASTQ files from [GEO GSE62944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944) (or SRA-derived accessions), place them under `data/raw/` using the filenames in the manifest, or edit the manifest to match your file names.

**SRA Toolkit example** (replace accessions with your chosen samples):

```bash
conda install -c bioconda sra-tools
prefetch SRR_ACCESSION
fasterq-dump --split-files --gzip SRR_ACCESSION
# Move/rename outputs to match config/samples.tsv
```

### Scaling to a full cohort

Add rows to `config/samples.tsv` and ensure FASTQ paths in `data/raw/` are correct. Disk and RAM requirements grow with sample count and read depth.

### Reference genome

```bash
mkdir -p data/genome
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz \
     -O data/genome/hg38.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz \
     -O data/genome/gencode.v44.annotation.gtf.gz
gunzip data/genome/*.gz
```

STAR indexing runs automatically on the first pipeline execution (`rules/align.smk`).

---

## Pipeline Steps

```
FASTQ (raw reads)
    │
    ▼
[FastQC] ──────────────────► Per-sample QC (HTML)
    │
    ▼
[Trimmomatic] ─────────────► Trimmed FASTQ
    │
    ▼
[MultiQC] ─────────────────► Aggregate QC report
    │
    ▼
[STAR] ────────────────────► Coordinate-sorted BAM
    │
    ▼
[featureCounts] ───────────► Gene-level count matrix
    │
    ▼
[DESeq2] ──────────────────► DEG table + normalized counts + plots
    │
    ▼
[clusterProfiler GSEA] ────► Pathway enrichment table + dotplot
    │
    ▼
results/figures/  results/tables/  results/qc/
```

---

## Output Files (after a successful run)

| File | Description |
|---|---|
| `results/qc/multiqc_report.html` | Aggregate QC |
| `results/tables/raw_counts.tsv` | Gene-level counts |
| `results/tables/deseq2_results.tsv` | DEG statistics and direction labels |
| `results/tables/normalized_counts.tsv` | VST-normalized matrix |
| `results/tables/gsea_results.tsv` | GSEA pathway table |
| `results/figures/volcano_plot.png` | Volcano plot |
| `results/figures/pca_plot.png` | Sample PCA |
| `results/figures/heatmap_top50.png` | Top DEG heatmap (or placeholder if few DEGs) |
| `results/figures/ma_plot.png` | MA plot |
| `results/figures/gsea_dotplot.png` | GSEA dotplot |

---

## Methods (workflow specification)

When executed, the pipeline applies the following methods (versions pinned in `envs/*.yaml`):

- **QC:** FastQC; MultiQC aggregation  
- **Trimming:** Trimmomatic (ILLUMINACLIP, LEADING/TRAILING, SLIDINGWINDOW, MINLEN per `config.yaml`)  
- **Alignment:** STAR to GRCh38 with GENCODE v44 GTF  
- **Quantification:** featureCounts (paired-end; strand mode set in `config.yaml`)  
- **Differential expression:** DESeq2 Wald test; \|log2FC\| and FDR thresholds in `config.yaml`; apeglm LFC shrinkage when available  
- **Enrichment:** GSEA via clusterProfiler using MSigDB Hallmark gene sets (`msigdbr` downloads sets at runtime)

---

## Requirements

- Linux, macOS, or **WSL** (recommended on Windows)
- conda or mamba
- Snakemake ≥ 7.0 (installed via `setup.sh`)
- ~32 GB RAM recommended for STAR genome indexing
- Disk: reference genome (~3 GB) plus FASTQ storage (dataset-dependent)

---

## Author

**Mohamed Elsaid**

---

## License

MIT License — see [LICENSE](LICENSE) (Copyright © 2025 Mohamed Elsaid).
