#!/usr/bin/env bash
# check_strandedness.sh: infer library strandedness by counting assigned fragments
# under each featureCounts strand setting on two aligned BAMs.
#
# Unstranded libraries assign about equal numbers of fragments under -s 1 and -s 2,
# and more under -s 0 (which accepts either orientation). A stranded library assigns
# most fragments under exactly one of -s 1 or -s 2. Set featurecounts.strand in
# config/config.yaml from the result: a wrong value silently skews the counts.
#
# Usage: bash scripts/check_strandedness.sh BAM1 BAM2 GTF OUT.tsv
# Requires featureCounts (subread) on PATH.

set -euo pipefail

BAM1="$1"; BAM2="$2"; GTF="$3"; OUT="$4"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

printf 'strand_setting\tassigned_fragments\tunassigned_no_features\tunassigned_ambiguity\n' > "$OUT"
for s in 0 1 2; do
  featureCounts -T 4 -p --countReadPairs -s "$s" -Q 10 -t exon -g gene_id \
    -a "$GTF" -o "$TMP/s$s.txt" "$BAM1" "$BAM2" 2> /dev/null
  # Sum the two BAM columns of each summary row
  awk -v s="$s" -F'\t' '
    $1 == "Assigned"              { a = $2 + $3 }
    $1 == "Unassigned_NoFeatures" { n = $2 + $3 }
    $1 == "Unassigned_Ambiguity"  { m = $2 + $3 }
    END { printf "%s\t%d\t%d\t%d\n", s, a, n, m }' "$TMP/s$s.txt.summary" >> "$OUT"
done
cat "$OUT"
