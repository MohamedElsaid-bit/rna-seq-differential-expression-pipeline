#!/usr/bin/env Rscript
# =============================================================================
# DESeq2 Differential Expression Analysis
# Pipeline: RNA-seq Differential Expression Pipeline
# Input:  raw_counts.tsv + sample metadata
# Output: DEG table, normalized counts, volcano, PCA, heatmap, MA plot
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(dplyr)
  library(tibble)
})

# ── Snakemake passes params via the snakemake object ─────────────────────────
counts_file    <- snakemake@input[["counts"]]
metadata_file  <- snakemake@input[["metadata"]]
ref_level      <- snakemake@params[["ref_level"]]
lfc_threshold  <- as.numeric(snakemake@params[["lfc_threshold"]])
fdr_threshold  <- as.numeric(snakemake@params[["fdr_threshold"]])
min_counts     <- as.integer(snakemake@params[["min_counts"]])
results_dir    <- snakemake@params[["results_dir"]]

log_file       <- snakemake@log[[1]]
log_con        <- file(log_file, open = "wt")
sink(log_con, append = TRUE, type = "output")
sink(log_con, append = TRUE, type = "message")


# ── 1. Load data ─────────────────────────────────────────────────────────────

message("Loading count matrix...")
counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
# Drop featureCounts metadata columns (Chr, Start, End, Strand, Length)
counts <- counts[, !colnames(counts) %in% c("Chr", "Start", "End", "Strand", "Length")]
counts <- as.matrix(counts)

message("Loading sample metadata...")
metadata <- read.delim(metadata_file, row.names = 1)
metadata$condition <- factor(metadata$condition,
                             levels = c(ref_level,
                                        setdiff(unique(metadata$condition), ref_level)))

# Ensure column order matches metadata row order
counts <- counts[, rownames(metadata)]
message(sprintf("Loaded %d genes x %d samples", nrow(counts), ncol(counts)))


# ── 2. Create DESeqDataSet ────────────────────────────────────────────────────

message("Building DESeqDataSet...")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ condition
)

# Pre-filter: remove genes with very low total counts
keep <- rowSums(counts(dds)) >= min_counts
dds  <- dds[keep, ]
message(sprintf("After filtering: %d genes retained", nrow(dds)))


# ── 3. Run DESeq2 ─────────────────────────────────────────────────────────────

message("Running DESeq2 (Wald test)...")
dds <- DESeq(dds)

res <- results(dds,
               contrast       = c("condition", "tumor", ref_level),
               alpha          = fdr_threshold,
               lfcThreshold   = 0)

# Shrink LFC estimates with apeglm for accurate visualization
coef_name <- resultsNames(dds)[grepl("condition_", resultsNames(dds)) &
                               !grepl("Intercept", resultsNames(dds))][1]
message(sprintf("lfcShrink coef: %s", coef_name))

res_shrunk <- lfcShrink(dds,
                        coef = coef_name,
                        type = "apeglm",
                        res  = res)

summary(res_shrunk)

# Convert to data frame and add gene significance labels
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj) %>%
  mutate(
    significant = !is.na(padj) &
                  padj < fdr_threshold &
                  abs(log2FoldChange) > lfc_threshold,
    direction   = case_when(
      significant & log2FoldChange >  lfc_threshold ~ "Upregulated",
      significant & log2FoldChange < -lfc_threshold ~ "Downregulated",
      TRUE                                           ~ "Not significant"
    )
  )

n_up   <- sum(res_df$direction == "Upregulated",   na.rm = TRUE)
n_down <- sum(res_df$direction == "Downregulated", na.rm = TRUE)
message(sprintf("DEGs: %d upregulated, %d downregulated", n_up, n_down))


# ── 4. Save result tables ─────────────────────────────────────────────────────

message("Writing result tables...")
write.table(res_df,
            file      = snakemake@output[["deseq2_results"]],
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

# VST-normalized counts for visualization
vst_counts <- vst(dds, blind = FALSE)
norm_mat    <- assay(vst_counts)
write.table(as.data.frame(norm_mat) %>% rownames_to_column("gene_id"),
            file      = snakemake@output[["normalized_counts"]],
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)


# ── 5. Volcano plot ───────────────────────────────────────────────────────────

message("Generating volcano plot...")

# Label top 15 DEGs by significance
top_genes <- res_df %>%
  filter(significant) %>%
  slice_min(padj, n = 15) %>%
  pull(gene_id)

pal <- c("Upregulated" = "#D85A30", "Downregulated" = "#378ADD", "Not significant" = "#888780")

p_volcano <- ggplot(res_df %>% filter(!is.na(padj)),
                    aes(x = log2FoldChange,
                        y = -log10(padj),
                        color = direction)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_text_repel(
    data   = filter(res_df, gene_id %in% top_genes),
    aes(label = gene_id),
    size   = 2.8,
    color  = "black",
    box.padding = 0.4,
    max.overlaps = 20
  ) +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = -log10(fdr_threshold),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = pal) +
  labs(
    title    = "Differential Gene Expression: Tumor vs. Normal",
    subtitle = sprintf("%d upregulated  |  %d downregulated  (|LFC| > %g, FDR < %g)",
                       n_up, n_down, lfc_threshold, fdr_threshold),
    x        = expression(log[2]~"Fold Change"),
    y        = expression(-log[10]~"Adjusted p-value"),
    color    = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey40", size = 10),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(snakemake@output[["volcano_plot"]],
       plot = p_volcano, width = 8, height = 6, dpi = 300)


# ── 6. PCA plot ───────────────────────────────────────────────────────────────

message("Generating PCA plot...")

pca_data  <- plotPCA(vst_counts, intgroup = c("condition", "batch"),
                     returnData = TRUE)
pct_var   <- round(100 * attr(pca_data, "percentVar"), 1)

p_pca <- ggplot(pca_data,
                aes(x = PC1, y = PC2,
                    color = condition, shape = batch)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = name), size = 2.5, color = "grey30") +
  scale_color_manual(values = c("tumor" = "#D85A30", "normal" = "#378ADD")) +
  labs(
    title  = "PCA of VST-Normalized Expression",
    x      = sprintf("PC1: %g%% variance", pct_var[1]),
    y      = sprintf("PC2: %g%% variance", pct_var[2]),
    color  = "Condition",
    shape  = "Batch"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(snakemake@output[["pca_plot"]],
       plot = p_pca, width = 7, height = 6, dpi = 300)


# ── 7. Heatmap (top 50 DEGs) ─────────────────────────────────────────────────

message("Generating heatmap...")

n_top <- min(50, sum(res_df$significant, na.rm = TRUE))
top50_genes <- res_df %>%
  filter(significant) %>%
  slice_min(padj, n = n_top) %>%
  pull(gene_id)

if (length(top50_genes) < 3) {
  message("WARNING: Fewer than 3 significant DEGs — skipping heatmap.")
  png(snakemake@output[["heatmap"]], width = 600, height = 400, res = 100)
  plot.new()
  text(0.5, 0.5, "Insufficient DEGs for heatmap.\nTry relaxing FDR or LFC thresholds.",
       cex = 1.2, col = "grey40")
  dev.off()
} else {
  heatmap_mat <- norm_mat[top50_genes, , drop = FALSE]
  # Scale rows (z-score) for visualization
  heatmap_mat <- t(scale(t(heatmap_mat)))

  ann_col <- data.frame(
    Condition = metadata$condition,
    row.names = rownames(metadata)
  )
  ann_colors <- list(Condition = c("tumor" = "#D85A30", "normal" = "#378ADD"))

  png(snakemake@output[["heatmap"]], width = 900, height = 1100, res = 150)
  pheatmap(heatmap_mat,
           annotation_col  = ann_col,
           annotation_colors = ann_colors,
           color           = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           cluster_rows    = TRUE,
           cluster_cols    = TRUE,
           show_rownames   = TRUE,
           show_colnames   = TRUE,
           fontsize_row    = 7,
           fontsize_col    = 9,
           main            = "Top 50 Differentially Expressed Genes (z-score)")
  dev.off()
}


# ── 8. MA plot ────────────────────────────────────────────────────────────────

message("Generating MA plot...")

p_ma <- ggplot(res_df %>% filter(!is.na(padj)),
               aes(x = log10(baseMean + 1),
                   y = log2FoldChange,
                   color = significant)) +
  geom_point(size = 0.6, alpha = 0.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(-lfc_threshold, lfc_threshold),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = c("TRUE" = "#D85A30", "FALSE" = "#888780"),
                     labels = c("TRUE" = "Significant DEG", "FALSE" = "Not significant")) +
  labs(
    title  = "MA Plot: Mean Expression vs. Fold Change",
    x      = expression(log[10]~"(Mean Expression + 1)"),
    y      = expression(log[2]~"Fold Change"),
    color  = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "top")

ggsave(snakemake@output[["ma_plot"]],
       plot = p_ma, width = 7, height = 5, dpi = 300)


message("DESeq2 analysis complete.")
print(sessionInfo())
sink()
sink(type = "message")
