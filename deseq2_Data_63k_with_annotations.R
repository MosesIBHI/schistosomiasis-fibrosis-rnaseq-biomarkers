#!/usr/bin/env Rscript

# ============================================================
# DESeq2 differential expression analysis
# Contrast: Sm_Hf vs n.c
# Input adapted for: Data_63k.csv
# Expected structure:
#   Column 1 = Geneid
#   Columns 2:41 = raw count samples
#   Groups inferred from column names:
#     Sm_Hf, Sm_Hf.1 ... Sm_Hf.9
#     Sm,    Sm.1    ... Sm.9
#     Hf,    Hf.1    ... Hf.9
#     n.c,   n.c.1   ... n.c.9
#
# Additional annotation module:
#   Ensembl Gene ID ENSG... -> gene symbol / gene name / Entrez ID
#   Creates additional annotated tables and annotated plots
#   Does NOT overwrite the previous main result files.
# ============================================================

# -----------------------------
# 0. Required packages
# -----------------------------
required_pkgs <- c(
  "DESeq2",
  "ggplot2",
  "pheatmap",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "ggrepel"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install them before running this script. Example:\n",
    "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
    "BiocManager::install(c('DESeq2', 'AnnotationDbi', 'org.Hs.eg.db'))\n",
    "install.packages(c('ggplot2', 'pheatmap', 'ggrepel'))"
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggrepel)
})

# -----------------------------
# 1. User-defined input/output
# -----------------------------
# The uploaded file is named Data_63k.csv.
# The script also checks a few common variants to avoid path/name errors.
candidate_files <- c("Data_63k.csv", "data_63.csv")
input_counts <- candidate_files[file.exists(candidate_files)][1]

if (is.na(input_counts)) {
  stop(
    "Input count file not found.\n",
    "Expected one of: ", paste(candidate_files, collapse = ", "), "\n",
    "Place the CSV file in the working directory or update 'candidate_files'."
  )
}

output_dir <- "results_deseq2_Data_63k"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Read count matrix
# -----------------------------
counts_df <- read.csv(
  input_counts,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(counts_df) == 0 || ncol(counts_df) < 2) {
  stop("The input file is empty or does not contain count columns.")
}

# Detect gene identifier column
gene_col <- colnames(counts_df)[1]

if (!gene_col %in% c("Geneid", "gene", "Gene", "gene_id", "GeneID")) {
  message("Using first column as gene identifier: ", gene_col)
}

if (any(is.na(counts_df[[gene_col]]) | counts_df[[gene_col]] == "")) {
  stop("The gene identifier column contains empty or NA values.")
}

# If duplicated gene IDs exist, aggregate counts by gene ID.
# This prevents row.names duplication errors.
if (anyDuplicated(counts_df[[gene_col]]) > 0) {
  message("Duplicated gene IDs detected. Aggregating counts by gene ID using sum().")
  count_cols <- setdiff(colnames(counts_df), gene_col)
  counts_df[count_cols] <- lapply(counts_df[count_cols], function(x) as.numeric(as.character(x)))
  counts_df <- aggregate(
    counts_df[count_cols],
    by = list(Geneid = counts_df[[gene_col]]),
    FUN = sum,
    na.rm = TRUE
  )
  gene_col <- "Geneid"
}

rownames(counts_df) <- counts_df[[gene_col]]
counts_df[[gene_col]] <- NULL

# Convert all sample columns to numeric counts
counts_df[] <- lapply(counts_df, function(x) as.numeric(as.character(x)))

if (any(is.na(as.matrix(counts_df)))) {
  stop("Count matrix contains NA values after numeric conversion. Please check the input file.")
}

counts <- as.matrix(counts_df)

if (any(counts < 0)) {
  stop("Count matrix contains negative values, which are invalid for DESeq2.")
}

if (any(counts != round(counts))) {
  stop("Count matrix contains non-integer values. DESeq2 requires raw integer counts.")
}

storage.mode(counts) <- "integer"

if (ncol(counts) == 0 || nrow(counts) == 0) {
  stop("Count matrix is empty after reading the file.")
}

# ------------------------------------------
# 3. Derive sample groups from real column names
# ------------------------------------------
sample_names <- colnames(counts)

# Remove suffixes like .1, .2, ..., .9 to recover group labels
condition_raw <- sub("\\.[0-9]+$", "", sample_names)

valid_groups <- c("Sm_Hf", "Sm", "Hf", "n.c")

if (!all(condition_raw %in% valid_groups)) {
  stop(
    "Unexpected group names detected in column names.\n",
    "Detected groups: ", paste(sort(unique(condition_raw)), collapse = ", "), "\n",
    "Expected groups: ", paste(valid_groups, collapse = ", "), "\n",
    "Current sample names: ", paste(sample_names, collapse = ", ")
  )
}

condition <- factor(
  condition_raw,
  levels = c("n.c", "Hf", "Sm", "Sm_Hf")
)

coldata <- data.frame(
  row.names = sample_names,
  condition = condition
)

# Check sample/metadata alignment
if (!identical(colnames(counts), rownames(coldata))) {
  stop("Sample names in counts and coldata are not aligned.")
}

# Check replicate numbers
group_table <- table(condition)
message("Detected sample groups:")
message(paste(capture.output(print(group_table)), collapse = "\n"))

if (any(group_table < 2)) {
  stop("At least one group has fewer than 2 replicates. DESeq2 requires biological replication.")
}

# For Data_63k.csv, this should be 40 samples: 10 per group.
expected_group_sizes <- c("n.c" = 10, "Hf" = 10, "Sm" = 10, "Sm_Hf" = 10)

if (!all(names(expected_group_sizes) %in% names(group_table)) ||
    !all(group_table[names(expected_group_sizes)] == expected_group_sizes)) {
  warning(
    "The detected group sizes differ from the expected Data_63k.csv structure.\n",
    "Detected: ", paste(names(group_table), as.integer(group_table), sep = "=", collapse = ", "), "\n",
    "Expected: ", paste(names(expected_group_sizes), expected_group_sizes, sep = "=", collapse = ", ")
  )
}

# -------------------------------------
# 4. Filter low-count genes
# -------------------------------------
# Keep genes with total raw counts > 10 across all samples
counts_filtered <- counts[rowSums(counts) > 10, , drop = FALSE]

if (nrow(counts_filtered) == 0) {
  stop("No genes remaining after filtering with rowSums(counts) > 10.")
}

# -------------------------------------
# 5. Build DESeq2 dataset and run model
# -------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)

# Variance stabilizing transformation for exploratory plots
vsd <- vst(dds, blind = FALSE)

# -------------------------------------
# 6. Main contrast: Sm_Hf vs n.c
# -------------------------------------

res <- results(
  dds,
  contrast = c("condition", "Sm_Hf", "n.c"),
  independentFiltering = FALSE
)

res <- res[order(res$padj), ]

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]

sig_df <- subset(res_df, !is.na(padj) & padj < 0.2)
sig_lfc_df <- subset(sig_df, abs(log2FoldChange) >= 1)

# Add direction column
res_df$direction <- "Not_significant"
res_df$direction[!is.na(res_df$padj) & res_df$padj < 0.2 & res_df$log2FoldChange > 0] <- "Up_in_Sm_Hf"
res_df$direction[!is.na(res_df$padj) & res_df$padj < 0.2 & res_df$log2FoldChange < 0] <- "Down_in_Sm_Hf"

sig_df <- subset(res_df, !is.na(padj) & padj < 0.2)
sig_lfc_df <- subset(sig_df, abs(log2FoldChange) >= 1)

# -------------------------------------
# 7. Save result tables
# -------------------------------------
write.csv(
  res_df,
  file = file.path(output_dir, "deseq2_full_results_Sm_Hf_vs_nc.csv"),
  row.names = FALSE
)

write.csv(
  sig_df,
  file = file.path(output_dir, "deseq2_significant_padj_lt_0.2_Sm_Hf_vs_nc.csv"),
  row.names = FALSE
)

write.csv(
  sig_lfc_df,
  file = file.path(output_dir, "deseq2_significant_padj_lt_0.2_absLFC_ge_1_Sm_Hf_vs_nc.csv"),
  row.names = FALSE
)

# -------------------------------------
# 8. PCA plot
# -------------------------------------
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA of variance-stabilized counts") +
  theme_bw()

ggsave(
  filename = file.path(output_dir, "pca_Sm_Hf_vs_nc.png"),
  plot = p_pca,
  width = 7,
  height = 5,
  dpi = 300
)

# -------------------------------------
# 9. Dispersion plot
# -------------------------------------
png(
  filename = file.path(output_dir, "dispersion_estimates.png"),
  width = 1800,
  height = 1400,
  res = 220
)
plotDispEsts(dds)
dev.off()

# -------------------------------------
# 10. MA plot
# -------------------------------------
png(
  filename = file.path(output_dir, "ma_plot_Sm_Hf_vs_nc.png"),
  width = 1800,
  height = 1400,
  res = 220
)
plotMA(res, ylim = c(-5, 5), main = "MA plot: Sm_Hf vs n.c")
dev.off()

# -------------------------------------
# 11. Volcano plot
# -------------------------------------
volcano_df <- res_df

# Avoid Inf values when pvalue is exactly 0
volcano_df$pvalue_for_plot <- volcano_df$pvalue
volcano_df$pvalue_for_plot[is.na(volcano_df$pvalue_for_plot)] <- NA
volcano_df$pvalue_for_plot[volcano_df$pvalue_for_plot == 0] <- .Machine$double.xmin
volcano_df$negLog10P <- -log10(volcano_df$pvalue_for_plot)

volcano_df$category <- "Not significant"
volcano_df$category[!is.na(volcano_df$padj) & volcano_df$padj < 0.2] <- "padj < 0.2"
volcano_df$category[
  !is.na(volcano_df$padj) &
    volcano_df$padj < 0.2 &
    abs(volcano_df$log2FoldChange) > 1
] <- "padj < 0.2 & |log2FC| > 1"

volcano_df$category <- factor(
  volcano_df$category,
  levels = c("Not significant", "padj < 0.2", "padj < 0.2 & |log2FC| > 1")
)

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = negLog10P)) +
  geom_point(aes(color = category), alpha = 0.7, size = 1.8, na.rm = TRUE) +
  scale_color_manual(
    values = c(
      "Not significant" = "grey70",
      "padj < 0.2" = "blue",
      "padj < 0.2 & |log2FC| > 1" = "red"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Volcano plot: Sm_Hf vs n.c",
    x = "log2 fold change",
    y = "-log10(p-value)",
    color = NULL
  ) +
  theme_bw()

ggsave(
  filename = file.path(output_dir, "volcano_Sm_Hf_vs_nc.png"),
  plot = p_volcano,
  width = 7,
  height = 5,
  dpi = 300
)

# -------------------------------------
# 12. Heatmap
# -------------------------------------
# Use top significant genes if available; otherwise top variable genes
vsd_mat <- assay(vsd)

if (nrow(sig_df) >= 2) {
  top_genes <- head(sig_df$gene, 50)
  top_genes <- top_genes[top_genes %in% rownames(vsd_mat)]
  heatmap_mat <- vsd_mat[top_genes, , drop = FALSE]
} else {
  gene_vars <- apply(vsd_mat, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(50, length(gene_vars))]
  heatmap_mat <- vsd_mat[top_genes, , drop = FALSE]
}

if (nrow(heatmap_mat) < 2) {
  warning("Heatmap skipped because fewer than 2 genes are available.")
} else {
  # Row scaling for visualization
  heatmap_mat <- t(scale(t(heatmap_mat)))
  heatmap_mat[is.na(heatmap_mat)] <- 0

  annotation_col <- data.frame(condition = droplevels(coldata$condition))
  rownames(annotation_col) <- rownames(coldata)

  # Explicit annotation colors prevents common pheatmap annotation/fill errors.
  annotation_colors <- list(
    condition = c(
      "n.c"   = "grey50",
      "Hf"    = "darkorange",
      "Sm"    = "steelblue",
      "Sm_Hf" = "firebrick"
    )
  )

  png(
    filename = file.path(output_dir, "heatmap_top_genes_Sm_Hf_vs_nc.png"),
    width = 2200,
    height = 1800,
    res = 220
  )

  pheatmap(
    heatmap_mat,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 7,
    main = "Heatmap of top genes"
  )

  dev.off()
}

# -------------------------------------
# 13. Save normalized counts
# -------------------------------------
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df$gene <- rownames(norm_counts_df)
norm_counts_df <- norm_counts_df[, c("gene", setdiff(colnames(norm_counts_df), "gene"))]

write.csv(
  norm_counts_df,
  file = file.path(output_dir, "normalized_counts.csv"),
  row.names = FALSE
)

# ============================================================
# 13B. Annotation module: Ensembl Gene ID -> reference gene
#      This module creates additional annotated files.
#      It does not overwrite the existing DESeq2 output files.
# ============================================================

# ------------------------------------------------------------
# 13B.1. Prepare Ensembl IDs
# ------------------------------------------------------------
res_df_annot <- res_df

# Remove possible Ensembl version suffixes such as ENSG00000123456.5
res_df_annot$gene_ensembl <- sub("\\..*$", "", res_df_annot$gene)
unique_ensg <- unique(res_df_annot$gene_ensembl)

# ------------------------------------------------------------
# 13B.2. Map Ensembl IDs to gene symbols, gene names and Entrez IDs
# ------------------------------------------------------------
gene_symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique_ensg,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_name_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique_ensg,
  column = "GENENAME",
  keytype = "ENSEMBL",
  multiVals = "first"
)

entrez_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique_ensg,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res_df_annot$gene_symbol <- unname(gene_symbol_map[res_df_annot$gene_ensembl])
res_df_annot$gene_name   <- unname(gene_name_map[res_df_annot$gene_ensembl])
res_df_annot$entrez_id   <- unname(entrez_map[res_df_annot$gene_ensembl])

# Label used in plots: gene symbol if available, otherwise Ensembl Gene ID
res_df_annot$plot_label <- ifelse(
  !is.na(res_df_annot$gene_symbol) & res_df_annot$gene_symbol != "",
  res_df_annot$gene_symbol,
  res_df_annot$gene_ensembl
)

# Reorder columns for readability
res_df_annot <- res_df_annot[, c(
  "gene",
  "gene_ensembl",
  "gene_symbol",
  "gene_name",
  "entrez_id",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj",
  "direction",
  "plot_label"
)]

# ------------------------------------------------------------
# 13B.3. Save annotated DESeq2 result tables
# ------------------------------------------------------------
write.csv(
  res_df_annot,
  file = file.path(output_dir, "deseq2_full_results_Sm_Hf_vs_nc_annotated_gene_symbols.csv"),
  row.names = FALSE
)

sig_df_annot <- subset(res_df_annot, !is.na(padj) & padj < 0.2)

write.csv(
  sig_df_annot,
  file = file.path(output_dir, "deseq2_significant_padj_lt_0.2_Sm_Hf_vs_nc_annotated_gene_symbols.csv"),
  row.names = FALSE
)

sig_lfc_df_annot <- subset(
  sig_df_annot,
  abs(log2FoldChange) >= 1
)

write.csv(
  sig_lfc_df_annot,
  file = file.path(output_dir, "deseq2_significant_padj_lt_0.2_absLFC_ge_1_Sm_Hf_vs_nc_annotated_gene_symbols.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13B.4. Save annotated normalized counts
# ------------------------------------------------------------
norm_counts_df_annot <- norm_counts_df
norm_counts_df_annot$gene_ensembl <- sub("\\..*$", "", norm_counts_df_annot$gene)

norm_counts_df_annot$gene_symbol <- unname(
  gene_symbol_map[norm_counts_df_annot$gene_ensembl]
)

norm_counts_df_annot$gene_name <- unname(
  gene_name_map[norm_counts_df_annot$gene_ensembl]
)

norm_counts_df_annot <- norm_counts_df_annot[, c(
  "gene",
  "gene_ensembl",
  "gene_symbol",
  "gene_name",
  setdiff(colnames(norm_counts_df_annot), c("gene", "gene_ensembl", "gene_symbol", "gene_name"))
)]

write.csv(
  norm_counts_df_annot,
  file = file.path(output_dir, "normalized_counts_annotated_gene_symbols.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13B.5. Annotated volcano plot with gene symbols
# ------------------------------------------------------------
volcano_annot_df <- res_df_annot

volcano_annot_df$pvalue_for_plot <- volcano_annot_df$pvalue
volcano_annot_df$pvalue_for_plot[is.na(volcano_annot_df$pvalue_for_plot)] <- NA
volcano_annot_df$pvalue_for_plot[volcano_annot_df$pvalue_for_plot == 0] <- .Machine$double.xmin
volcano_annot_df$negLog10P <- -log10(volcano_annot_df$pvalue_for_plot)

volcano_annot_df$category <- "Not significant"
volcano_annot_df$category[
  !is.na(volcano_annot_df$padj) &
    volcano_annot_df$padj < 0.2
] <- "padj < 0.2"

volcano_annot_df$category[
  !is.na(volcano_annot_df$padj) &
    volcano_annot_df$padj < 0.2 &
    abs(volcano_annot_df$log2FoldChange) > 1
] <- "padj < 0.2 & |log2FC| > 1"

volcano_annot_df$category <- factor(
  volcano_annot_df$category,
  levels = c("Not significant", "padj < 0.2", "padj < 0.2 & |log2FC| > 1")
)

# Genes of interest to force-label on volcano and MA plots.
# Add or remove genes here if needed.
candidate_genes_to_label <- c(
  "ENSG00000132849",
  "ENSG00000176148"
)

# Label top significant genes by adjusted p-value and always include candidates.
label_df <- subset(
  volcano_annot_df,
  !is.na(padj) & padj < 0.2 & abs(log2FoldChange) >= 1
)

label_df <- label_df[order(label_df$padj), , drop = FALSE]
label_df <- head(label_df, 20)

candidate_label_df <- subset(
  volcano_annot_df,
  gene_ensembl %in% candidate_genes_to_label
)

label_df <- unique(rbind(label_df, candidate_label_df))

p_volcano_annot <- ggplot(
  volcano_annot_df,
  aes(x = log2FoldChange, y = negLog10P)
) +
  geom_point(aes(color = category), alpha = 0.7, size = 1.8, na.rm = TRUE) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = plot_label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = c(
      "Not significant" = "grey70",
      "padj < 0.2" = "blue",
      "padj < 0.2 & |log2FC| > 1" = "red"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Annotated volcano plot: Sm_Hf vs n.c",
    x = "log2 fold change",
    y = "-log10(p-value)",
    color = NULL
  ) +
  theme_bw()

ggsave(
  filename = file.path(output_dir, "volcano_annotated_gene_symbols_Sm_Hf_vs_nc.png"),
  plot = p_volcano_annot,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 13B.6. Annotated MA plot with gene symbols
# ------------------------------------------------------------
ma_annot_df <- volcano_annot_df
ma_annot_df$log10_baseMean <- log10(ma_annot_df$baseMean + 1)

ma_label_df <- subset(
  ma_annot_df,
  gene_ensembl %in% candidate_genes_to_label |
    (!is.na(padj) & padj < 0.2 & abs(log2FoldChange) >= 1)
)

ma_label_df <- ma_label_df[order(ma_label_df$padj), , drop = FALSE]
ma_label_df <- head(ma_label_df, 20)

p_ma_annot <- ggplot(
  ma_annot_df,
  aes(x = log10_baseMean, y = log2FoldChange)
) +
  geom_point(aes(color = category), alpha = 0.7, size = 1.5, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = ma_label_df,
    aes(label = plot_label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = c(
      "Not significant" = "grey70",
      "padj < 0.2" = "blue",
      "padj < 0.2 & |log2FC| > 1" = "red"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Annotated MA plot: Sm_Hf vs n.c",
    x = "log10(baseMean + 1)",
    y = "log2 fold change",
    color = NULL
  ) +
  theme_bw()

ggsave(
  filename = file.path(output_dir, "ma_plot_annotated_gene_symbols_Sm_Hf_vs_nc.png"),
  plot = p_ma_annot,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 13B.7. Annotated heatmap with gene symbols
# ------------------------------------------------------------
vsd_mat_annot <- assay(vsd)

if (nrow(sig_df_annot) >= 2) {
  heatmap_genes_annot <- head(sig_df_annot$gene, 50)
  heatmap_genes_annot <- heatmap_genes_annot[
    heatmap_genes_annot %in% rownames(vsd_mat_annot)
  ]
} else {
  gene_vars_annot <- apply(vsd_mat_annot, 1, var)
  heatmap_genes_annot <- names(sort(gene_vars_annot, decreasing = TRUE))[
    1:min(50, length(gene_vars_annot))
  ]
}

heatmap_mat_annot <- vsd_mat_annot[heatmap_genes_annot, , drop = FALSE]

if (nrow(heatmap_mat_annot) < 2) {
  warning("Annotated heatmap skipped because fewer than 2 genes are available.")
} else {
  heatmap_mat_annot <- t(scale(t(heatmap_mat_annot)))
  heatmap_mat_annot[is.na(heatmap_mat_annot)] <- 0

  heatmap_ensg <- sub("\\..*$", "", rownames(heatmap_mat_annot))
  heatmap_symbols <- unname(gene_symbol_map[heatmap_ensg])

  heatmap_row_labels <- ifelse(
    !is.na(heatmap_symbols) & heatmap_symbols != "",
    paste0(heatmap_symbols, " (", heatmap_ensg, ")"),
    heatmap_ensg
  )

  annotation_col <- data.frame(condition = droplevels(coldata$condition))
  rownames(annotation_col) <- rownames(coldata)

  annotation_colors <- list(
    condition = c(
      "n.c"   = "grey50",
      "Hf"    = "darkorange",
      "Sm"    = "steelblue",
      "Sm_Hf" = "firebrick"
    )
  )

  png(
    filename = file.path(output_dir, "heatmap_annotated_gene_symbols_Sm_Hf_vs_nc.png"),
    width = 2600,
    height = 2000,
    res = 220
  )

  pheatmap(
    heatmap_mat_annot,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    labels_row = heatmap_row_labels,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 7,
    main = "Annotated heatmap of top genes"
  )

  dev.off()
}

# ------------------------------------------------------------
# 13B.8. Candidate-gene annotated table for key transcripts
# ------------------------------------------------------------
candidate_results_annot <- subset(
  res_df_annot,
  gene_ensembl %in% candidate_genes_to_label
)

write.csv(
  candidate_results_annot,
  file = file.path(output_dir, "candidate_transcripts_annotated_gene_symbols.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13B.9. Annotation summary
# ------------------------------------------------------------
annotation_summary <- data.frame(
  total_genes_tested = nrow(res_df_annot),
  genes_with_symbol = sum(!is.na(res_df_annot$gene_symbol) & res_df_annot$gene_symbol != ""),
  genes_without_symbol = sum(is.na(res_df_annot$gene_symbol) | res_df_annot$gene_symbol == ""),
  significant_padj_lt_0.2_with_symbol = sum(!is.na(sig_df_annot$gene_symbol) & sig_df_annot$gene_symbol != "")
)

write.csv(
  annotation_summary,
  file = file.path(output_dir, "annotation_summary.csv"),
  row.names = FALSE
)

# -------------------------------------
# 14. Save sample metadata and QC summary
# -------------------------------------
write.csv(
  coldata,
  file = file.path(output_dir, "sample_metadata.csv"),
  row.names = TRUE
)

qc_summary <- data.frame(
  sample = colnames(counts),
  condition = as.character(condition),
  library_size_raw = colSums(counts),
  library_size_filtered = colSums(counts_filtered)
)

write.csv(
  qc_summary,
  file = file.path(output_dir, "qc_library_sizes.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

# -------------------------------------
# 15. Console summary
# -------------------------------------
cat("\nDESeq2 analysis completed successfully.\n")
cat("Input file: ", input_counts, "\n", sep = "")
cat("Output directory: ", output_dir, "\n", sep = "")
cat("Samples detected: ", ncol(counts), "\n", sep = "")
cat("Groups detected: ", paste(names(group_table), as.integer(group_table), sep = "=", collapse = ", "), "\n", sep = "")
cat("Genes before filtering: ", nrow(counts), "\n", sep = "")
cat("Genes after filtering: ", nrow(counts_filtered), "\n", sep = "")
cat("Total DE results: ", nrow(res_df), "\n", sep = "")
cat("Significant genes padj < 0.2: ", nrow(sig_df), "\n", sep = "")
cat("Significant genes padj < 0.2 and |log2FC| >= 1: ", nrow(sig_lfc_df), "\n", sep = "")
cat("Up in Sm_Hf, padj < 0.2: ", sum(sig_df$direction == "Up_in_Sm_Hf"), "\n", sep = "")
cat("Down in Sm_Hf, padj < 0.2: ", sum(sig_df$direction == "Down_in_Sm_Hf"), "\n", sep = "")

cat("\nGene annotation completed.\n")
cat("Annotated genes with symbol: ", annotation_summary$genes_with_symbol, " / ", annotation_summary$total_genes_tested, "\n", sep = "")
cat("Annotated full table: deseq2_full_results_Sm_Hf_vs_nc_annotated_gene_symbols.csv\n")
cat("Annotated volcano plot: volcano_annotated_gene_symbols_Sm_Hf_vs_nc.png\n")
cat("Annotated MA plot: ma_plot_annotated_gene_symbols_Sm_Hf_vs_nc.png\n")
cat("Annotated heatmap: heatmap_annotated_gene_symbols_Sm_Hf_vs_nc.png\n")
cat("Candidate annotated table: candidate_transcripts_annotated_gene_symbols.csv\n")
