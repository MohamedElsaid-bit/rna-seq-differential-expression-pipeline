#!/usr/bin/env bash
# setup.sh — One-time environment setup for rna-seq-pipeline
# Run this once before your first pipeline execution.
# Requires conda or mamba to be installed.

set -euo pipefail

CONDA_CMD="conda"
command -v mamba &>/dev/null && CONDA_CMD="mamba"
echo "Using: $CONDA_CMD"

echo "Creating Snakemake environment..."
$CONDA_CMD create -n snakemake -c conda-forge -c bioconda \
    snakemake=7.32 python=3.10 pandas=2.0 --yes

echo ""
echo "Setup complete."
echo "Activate with:  conda activate snakemake"
echo "Dry run with:   snakemake --dry-run --cores 1 --use-conda"
echo "Full run with:  snakemake --cores 8 --use-conda"
