#!/usr/bin/env Rscript
# =============================================================================
# Pathway Enrichment Analysis (GSEA)
# Input:  DESeq2 results table (ranked by stat)
# Output: GSEA dotplot + results table (MSigDB Hallmark gene sets)
# =============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(ggplot2)
  library(dplyr)
  library(enrichplot)
})

# ── Snakemake params ──────────────────────────────────────────────────────────
deseq2_file   <- snakemake@input[["deseq2_results"]]
pvalue_cutoff <- as.numeric(snakemake@params[["pvalue_cutoff"]])
gene_sets     <- snakemake@params[["gene_sets"]]          # "H" = Hallmark
organism      <- snakemake@params[["organism"]]            # "hsa"

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, append = TRUE, type = "output")
sink(log_con, append = TRUE, type = "message")


# ── 1. Build ranked gene list ─────────────────────────────────────────────────

message("Loading DESeq2 results...")
res_df <- read.delim(deseq2_file)

# Rank by DESeq2 Wald statistic (preferred over log2FC for GSEA)
ranked_genes <- res_df %>%
  filter(!is.na(stat)) %>%
  arrange(desc(stat)) %>%
  dplyr::select(gene_id, stat) %>%
  deframe()

message(sprintf("Ranked gene list: %d genes", length(ranked_genes)))


# ── 2. Load gene sets ─────────────────────────────────────────────────────────

message(sprintf("Loading MSigDB gene sets (category: %s, organism: %s)...",
                gene_sets, organism))

msig_df <- msigdbr(species = "Homo sapiens", category = gene_sets)
msig_t2g <- msig_df %>%
  dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame()

message(sprintf("Gene sets loaded: %d pathways", length(unique(msig_t2g$gs_name))))


# ── 3. Run GSEA ───────────────────────────────────────────────────────────────

if (length(ranked_genes) < 10) {
  stop("Ranked gene list has fewer than 10 genes. Check DESeq2 output and gene ID format.")
}

message("Running GSEA...")
set.seed(42)

gsea_result <- GSEA(
  geneList    = ranked_genes,
  TERM2GENE   = msig_t2g,
  pvalueCutoff = pvalue_cutoff,
  pAdjustMethod = "BH",
  minGSSize   = snakemake@params[["min_gs_size"]],
  maxGSSize   = snakemake@params[["max_gs_size"]],
  seed        = TRUE,
  verbose     = FALSE
)

n_sig <- nrow(as.data.frame(gsea_result))
message(sprintf("Significant pathways (FDR < %g): %d", pvalue_cutoff, n_sig))

if (nrow(as.data.frame(gsea_result)) == 0) {
  message("WARNING: No significant pathways found at current p-value cutoff.")
  p_empty <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "No significant pathways found.\nConsider relaxing pvalue_cutoff in config.yaml.",
             size = 4, color = "grey40", hjust = 0.5) +
    theme_void()
  ggsave(snakemake@output[["gsea_dotplot"]], plot = p_empty,
         width = 7, height = 4, dpi = 150)
  write.table(data.frame(), file = snakemake@output[["gsea_results"]],
              sep = "\t", row.names = FALSE)
  quit(save = "no", status = 0)
}


# ── 4. Save GSEA results table ────────────────────────────────────────────────

gsea_df <- as.data.frame(gsea_result) %>%
  arrange(p.adjust) %>%
  dplyr::select(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue)

write.table(gsea_df,
            file      = snakemake@output[["gsea_results"]],
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)


# ── 5. GSEA dotplot ───────────────────────────────────────────────────────────

message("Generating GSEA dotplot...")

# Clean pathway names for display
gsea_df_plot <- gsea_df %>%
  slice_min(p.adjust, n = 20) %>%
  mutate(
    Description = gsub("HALLMARK_", "", Description),
    Description = gsub("_", " ", Description),
    Direction   = ifelse(NES > 0, "Enriched in Tumor", "Enriched in Normal")
  )

p_gsea <- ggplot(gsea_df_plot,
                 aes(x    = NES,
                     y    = reorder(Description, NES),
                     size = setSize,
                     color = p.adjust)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_gradient(low = "#D85A30", high = "#B4B2A9",
                       name = "Adjusted\np-value") +
  scale_size_continuous(name = "Gene set\nsize", range = c(3, 9)) +
  labs(
    title    = "GSEA: Hallmark Pathway Enrichment",
    subtitle = "Tumor vs. Normal (TCGA-BRCA)",
    x        = "Normalized Enrichment Score (NES)",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "grey40"),
    panel.grid.minor = element_blank(),
    axis.text.y     = element_text(size = 9)
  )

ggsave(snakemake@output[["gsea_dotplot"]],
       plot = p_gsea, width = 9, height = 7, dpi = 300)


message("Pathway enrichment analysis complete.")
sink()
sink(type = "message")
