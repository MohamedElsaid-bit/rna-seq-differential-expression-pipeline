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
  library(tibble)      # deframe()
  library(enrichplot)
})

# ── Snakemake params ──────────────────────────────────────────────────────────
deseq2_file   <- snakemake@input[["deseq2_results"]]
pvalue_cutoff <- as.numeric(snakemake@params[["pvalue_cutoff"]])
gene_sets     <- snakemake@params[["gene_sets"]]          # "H" = Hallmark
organism      <- snakemake@params[["organism"]]            # "hsa"
ref_level     <- snakemake@params[["ref_level"]]
treat_level   <- snakemake@params[["treat_level"]]

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, append = TRUE, type = "output")
sink(log_con, append = TRUE, type = "message")


# ── 1. Build ranked gene list ─────────────────────────────────────────────────

message("Loading DESeq2 results...")
res_df <- read.delim(deseq2_file)

# Rank by DESeq2 Wald statistic (preferred over log2FC for GSEA)
# MSigDB gene sets are keyed by gene symbol, so rank by symbol. Where several Ensembl
# genes share a symbol, keep the one with the largest absolute statistic so that
# each symbol appears once in the ranked list (GSEA requires unique names).
ranked_genes <- res_df %>%
  filter(!is.na(stat), !is.na(gene_name), gene_name != "") %>%
  group_by(gene_name) %>%
  slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(stat)) %>%
  dplyr::select(gene_name, stat) %>%
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
  # Keep every tested gene set (cutoff 1); significance is flagged below using the
  # configured FDR cutoff, so the table stays informative when nothing passes.
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize   = snakemake@params[["min_gs_size"]],
  maxGSSize   = snakemake@params[["max_gs_size"]],
  seed        = TRUE,
  verbose     = FALSE
)

n_tested <- nrow(as.data.frame(gsea_result))
n_sig    <- sum(as.data.frame(gsea_result)$p.adjust < pvalue_cutoff)
message(sprintf("Gene sets tested: %d; significant (FDR < %g): %d",
                n_tested, pvalue_cutoff, n_sig))

if (n_tested == 0) {
  message("WARNING: No gene sets met the size limits, nothing was tested.")
  p_empty <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "No gene sets were testable.\nCheck min_gs_size and gene symbol overlap.",
             size = 4, color = "grey40", hjust = 0.5) +
    theme_void()
  ggsave(snakemake@output[["gsea_dotplot"]], plot = p_empty,
         width = 7, height = 4, dpi = 150, bg = "white")
  write.table(data.frame(), file = snakemake@output[["gsea_results"]],
              sep = "\t", row.names = FALSE)
  quit(save = "no", status = 0)
}


# ── 4. Save GSEA results table ────────────────────────────────────────────────

gsea_df <- as.data.frame(gsea_result) %>%
  arrange(p.adjust) %>%
  mutate(significant = p.adjust < pvalue_cutoff) %>%
  dplyr::select(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue,
                significant)

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
    Direction   = ifelse(NES > 0, paste("Enriched in", treat_level),
                         paste("Enriched in", ref_level))
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
  scale_y_discrete(expand = expansion(add = 1)) +   # room so the largest dots are not clipped
  labs(
    title    = "GSEA: Hallmark Pathway Enrichment",
    subtitle = sprintf("%s vs. %s  |  %d of %d gene sets significant (FDR < %g)",
                       treat_level, ref_level, n_sig, n_tested, pvalue_cutoff),
    x        = "Normalized Enrichment Score (NES)",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(color = "grey40"),
    panel.grid.minor = element_blank(),
    axis.text.y     = element_text(size = 9),
    plot.margin     = margin(10, 15, 15, 10)
  )

ggsave(snakemake@output[["gsea_dotplot"]],
       plot = p_gsea, width = 9, height = 7, dpi = 300, bg = "white")


message("Pathway enrichment analysis complete.")
sink()
sink(type = "message")
