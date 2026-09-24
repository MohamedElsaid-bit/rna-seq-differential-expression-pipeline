#!/usr/bin/env bash
# download_data.sh: fetch the scoped airway dataset (GEO GSE52778, ENA PRJNA229998).
#
# Reads: the first N_PAIRS read pairs of each run are streamed from ENA and the
#   download is stopped early, so the full FASTQ files are never stored. The head of
#   a FASTQ is not a random sample (limited flowcell tiles); this is documented in
#   the README.
# Reference: GRCh38 chromosomes and the GENCODE v44 annotation restricted to those
#   chromosomes, so STAR indexing fits in a few GB of RAM.
#
# Usage: bash scripts/download_data.sh [config/samples.tsv]

set -euo pipefail

SAMPLES="${1:-config/samples.tsv}"
N_PAIRS="${N_PAIRS:-5000000}"
CHROMS="${CHROMS:-chr5 chr6 chr17}"
RAW_DIR="data/raw"
GENOME_DIR="data/genome"
ENA="https://ftp.sra.ebi.ac.uk/vol1/fastq"
UCSC="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"

mkdir -p "$RAW_DIR" "$GENOME_DIR"

# ── Reference ────────────────────────────────────────────────────────────────
if [[ ! -s "$GENOME_DIR/hg38_subset.fa" ]]; then
  : > "$GENOME_DIR/hg38_subset.fa"
  for chrom in $CHROMS; do
    echo "Fetching $chrom"
    curl -fsSL "$UCSC/$chrom.fa.gz" | gzip -dc >> "$GENOME_DIR/hg38_subset.fa"
  done
fi

if [[ ! -s "$GENOME_DIR/gencode.v44.subset.gtf" ]]; then
  echo "Fetching GENCODE v44 annotation"
  pattern=$(echo "$CHROMS" | tr ' ' '|')
  curl -fsSL "$GTF_URL" | gzip -dc \
    | awk -F'\t' -v pat="^($pattern)$" '/^#/ || $1 ~ pat' \
    > "$GENOME_DIR/gencode.v44.subset.gtf"
fi

# ── Reads ────────────────────────────────────────────────────────────────────
fetch_subset() {
  # $1 = URL, $2 = output path. Keep the first N_PAIRS reads (4 lines each).
  # Retries because ENA connections occasionally time out mid-stream; a truncated
  # file is never kept (gzip integrity and exact line count are both checked).
  local url="$1" out="$2" tmp="$2.part" attempt lines
  for attempt in 1 2 3 4; do
    # head closes the pipe early, so curl and gzip exit with SIGPIPE by design.
    set +o pipefail
    curl -fsSL --connect-timeout 30 --speed-limit 10000 --speed-time 60 "$url" \
      | gzip -dc | head -n $((N_PAIRS * 4)) | gzip -c > "$tmp" 2>/dev/null
    set -o pipefail
    lines=$(gzip -dc "$tmp" 2>/dev/null | wc -l || true)
    if [[ "$lines" -eq $((N_PAIRS * 4)) ]] && gzip -t "$tmp" 2>/dev/null; then
      mv "$tmp" "$out"
      return 0
    fi
    echo "Attempt $attempt failed for $out ($lines of $((N_PAIRS * 4)) lines), retrying" >&2
    rm -f "$tmp"
  done
  echo "ERROR: could not fetch $out after 4 attempts" >&2
  exit 1
}

tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r sample_id condition donor srr r1 r2; do
  prefix="${srr:0:6}"
  last="${srr: -1}"
  for mate in 1 2; do
    out="$RAW_DIR/${sample_id}_R${mate}.fastq.gz"
    if [[ -s "$out" ]]; then echo "Skip $out"; continue; fi
    echo "Streaming $srr mate $mate -> $out"
    fetch_subset "$ENA/$prefix/00$last/$srr/${srr}_${mate}.fastq.gz" "$out"
  done
done

echo "Done."
