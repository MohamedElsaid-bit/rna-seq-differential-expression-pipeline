#!/usr/bin/env python3
"""
clean_count_matrix.py
--------------------
Cleans featureCounts output to strip full BAM paths from column headers,
leaving only clean sample IDs matching config/samples.tsv.

Usage:
    python scripts/clean_count_matrix.py results/tables/raw_counts.tsv
"""

import sys
import re
import pandas as pd
from pathlib import Path


def clean_column_name(col: str) -> str:
    """
    Extract sample ID from a full BAM path.
    e.g. 'results/alignments/TCGA_BR_001.Aligned.sortedByCoord.out.bam'
         -> 'TCGA_BR_001'
    """
    stem = Path(col).stem  # remove .bam
    # Remove STAR suffix
    stem = re.sub(r"\.Aligned\.sortedByCoord\.out$", "", stem)
    return stem


def main(counts_file: str) -> None:
    print(f"Cleaning count matrix: {counts_file}")

    df = pd.read_csv(counts_file, sep="\t", comment="#", index_col=0)

    # featureCounts first 5 cols after gene_id are metadata — keep as-is
    meta_cols = ["Chr", "Start", "End", "Strand", "Length"]
    meta_mask = df.columns.isin(meta_cols)

    # Rename BAM path columns
    new_columns = []
    for col in df.columns:
        if col in meta_cols:
            new_columns.append(col)
        else:
            new_columns.append(clean_column_name(col))

    df.columns = new_columns

    # Drop featureCounts metadata columns — keep only count columns
    count_cols = [c for c in df.columns if c not in meta_cols]
    df_counts = df[count_cols]

    df_counts.to_csv(counts_file, sep="\t")
    print(f"  Cleaned columns: {list(df_counts.columns)}")
    print(f"  Shape: {df_counts.shape[0]} genes x {df_counts.shape[1]} samples")
    print("Done.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python clean_count_matrix.py <counts.tsv>")
        sys.exit(1)
    main(sys.argv[1])
