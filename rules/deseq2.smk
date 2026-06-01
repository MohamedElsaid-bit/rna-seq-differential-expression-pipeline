"""
DESeq2 Rules: Differential expression analysis + GSEA pathway enrichment.
"""

rule deseq2:
    """
    Run DESeq2 differential expression analysis.
    Produces DEG table, normalized counts, and core visualizations.
    """
    input:
        counts   = f"{RESULTS}/tables/raw_counts.tsv",
        metadata = config["samples"],
    output:
        deseq2_results    = f"{RESULTS}/tables/deseq2_results.tsv",
        normalized_counts = f"{RESULTS}/tables/normalized_counts.tsv",
        volcano_plot      = f"{RESULTS}/figures/volcano_plot.png",
        pca_plot          = f"{RESULTS}/figures/pca_plot.png",
        heatmap           = f"{RESULTS}/figures/heatmap_top50.png",
        ma_plot           = f"{RESULTS}/figures/ma_plot.png",
    params:
        ref_level     = config["deseq2"]["reference_level"],
        lfc_threshold = config["deseq2"]["lfc_threshold"],
        fdr_threshold = config["deseq2"]["fdr_threshold"],
        min_counts    = config["deseq2"]["min_counts"],
        results_dir   = f"{RESULTS}",
    conda:
        "../envs/r_analysis.yaml"
    log:
        f"{RESULTS}/logs/deseq2/deseq2.log"
    script:
        "../scripts/deseq2_analysis.R"


rule pathway_enrichment:
    """
    Run GSEA with clusterProfiler using MSigDB Hallmark gene sets.
    Input: ranked gene list from DESeq2 (by stat or log2FC * -log10p).
    """
    input:
        deseq2_results = f"{RESULTS}/tables/deseq2_results.tsv",
    output:
        gsea_dotplot   = f"{RESULTS}/figures/gsea_dotplot.png",
        gsea_results   = f"{RESULTS}/tables/gsea_results.tsv",
    params:
        gene_sets     = config["gsea"]["gene_sets"],
        organism      = config["gsea"]["organism"],
        pvalue_cutoff = config["gsea"]["pvalue_cutoff"],
        min_gs_size   = config["gsea"]["min_gs_size"],
        max_gs_size   = config["gsea"]["max_gs_size"],
        results_dir   = f"{RESULTS}",
    conda:
        "../envs/r_analysis.yaml"
    log:
        f"{RESULTS}/logs/pathway_enrichment/gsea.log"
    script:
        "../scripts/pathway_enrichment.R"
