"""
Quantification Rule: Gene-level read counting with featureCounts.
"""

rule featurecounts:
    """
    Count reads per gene across all samples simultaneously.
    Produces a single count matrix TSV used as DESeq2 input.
    """
    input:
        bams = expand(
            "{results}/alignments/{sample}.Aligned.sortedByCoord.out.bam",
            results=RESULTS,
            sample=SAMPLE_IDS,
        ),
        gtf  = config["gtf_file"],
    output:
        counts   = f"{RESULTS}/tables/raw_counts.tsv",
        summary  = f"{RESULTS}/tables/raw_counts.tsv.summary",
        annotation = f"{RESULTS}/tables/gene_annotation.tsv",
    params:
        strand    = config["featurecounts"]["strand"],
        min_mapq  = config["featurecounts"]["min_mapq"],
        # subread 2.0.2 and later needs --countReadPairs to count fragments, not mates
        pair_end  = "-p --countReadPairs" if config["featurecounts"]["pair_end"] else "",
    threads: 8
    resources:
        mem_mb = 4000,
    conda:
        "../environment.yml"
    log:
        f"{RESULTS}/logs/featurecounts/featurecounts.log"
    shell:
        """
        featureCounts \
            -T {threads} \
            {params.pair_end} \
            -s {params.strand} \
            -Q {params.min_mapq} \
            -a {input.gtf} \
            -t exon -g gene_id --extraAttributes gene_name \
            -o {output.counts} \
            {input.bams} \
            2> {log}

        # Rename BAM columns to clean sample IDs for downstream R scripts
        python scripts/clean_count_matrix.py {output.counts}
        """
