###############################################################################
# EBSeq-based Differential Expression Analysis 
#
# Purpose:
#   Reanalyse RNA-seq dataset using EBSeq 
#   empirical Bayes framework to compare with the primary edgeR analysis.
#
# Notes:
#   - Input data are gene-level raw counts.
#   - A MAD-based feature filter is applied separately for each comparison.
#   - EBSeq median normalization factors are estimated separately for each
#     filtered subset.
#   - Significant results are exported as CSV tables.
###############################################################################

############################
# 1. SETUP
############################

# Working directory containing input RDS files and desired outputs.
setwd("xxxxxxxxx")
# Required packages.
library(EBSeq)#version '1.28.0'
library(pheatmap)#version '1.0.12'
library(edgeR)#version '4.0.16'

############################
# 2. IMPORT INPUT DATA
############################

# # Load count matrix
# 
Rcounts <- readRDS("ribo_minus_counts.rds")
dim(Rcounts)#52114 36

gene <- readRDS("ensembl_gchr37_annotation.rds")#to reload

# Remove genes with zero counts across all samples
Rcounts <- Rcounts[rowSums(Rcounts[,-1]) > 0, ]
dim(Rcounts)#52088  41

###############################
# 3. LOAD METADATA
###############################

sev_meta <- read.csv2("Metadata.csv",
                      sep = ",",
                      row.names = 1)

sev_meta <- sev_meta[colnames(Rcounts)[-1],]

# Ensure sample alignment
stopifnot(all(rownames(sev_meta) %in% colnames(Rcounts)))

# Reorder counts to match metadata
Rcounts <- data.frame(Rcounts[,1],Rcounts[,rownames(sev_meta)])
colnames(Rcounts)[1] <- "gene"
identical(rownames(sev_meta),colnames(Rcounts)[-1])#TRUE

############################
# 3. HELPER FUNCTIONS
############################

# Subset metadata and count matrix for a two-group comparison.
subset_comparison <- function(meta, counts, group_a, group_b) {
  meta_sub <- meta[meta$Group %in% c(group_a, group_b), , drop = FALSE]
  meta_sub$Group <- factor(as.character(meta_sub$Group))
  
  counts_sub <- counts[, rownames(meta_sub), drop = FALSE]
  
  list(meta = meta_sub, counts = counts_sub)
}

# Apply a MAD filter, keeping genes above the 75th percentile of MAD values.
mad_filter_counts <- function(counts_mat) {
  mad_values <- apply(as.matrix(counts_mat), 1, mad)
  cutoff <- unname(quantile(mad_values,na.rm=T)[4])
  keep <- names(which(mad_values > cutoff))
  
  list(
    counts = as.matrix(counts_mat)[keep, , drop = FALSE],
    mad = mad_values,
    cutoff = cutoff
  )
}

# Run EBSeq for a single comparison and export significant results.
run_ebseq_comparison <- function(counts_mat,
                                 conditions,
                                 size_factors,
                                 annotation,
                                 out_prefix,
                                 fdr = 0.1,
                                 threshold_fc = 0.8,
                                 maxround = 10) {
  
  eb_out <- EBSeq::EBTest(
    Data = as.matrix(counts_mat),
    Conditions = conditions,
    sizeFactors = size_factors,
    maxround = maxround
  )
  
  de_res <- GetDEResults(
    eb_out,
    FDR = fdr,
    Threshold_FC = threshold_fc
  )
  
  n_de <- sum(de_res$Status == "DE")
  message(out_prefix, ": ", n_de, " DE genes")
  
  sig_table <- NULL
  gene_fc <- PostFC(eb_out)
  
  if (n_de > 0) {
    sig_ids <- which(de_res$Status == "DE")
    sig_anno <- annotation[names(sig_ids), , drop = FALSE]
    
    pp_mat <- de_res$PPMat[de_res$Status == "DE", , drop = FALSE]
    
    # For a single significant gene, preserve row structure.
    if (nrow(pp_mat) == 1) {
      pp_mat <- as.data.frame(pp_mat)
    } else {
      pp_mat <- as.data.frame(na.omit(pp_mat))
    }
    
    sig_table <- data.frame(pp_mat, sig_anno)
    sig_table$realFC <- gene_fc$RealFC[rownames(sig_table)]
    
    write.csv(
      sig_table,
      paste0(out_prefix, "_EBseq_sig_table_MAD_filtered.csv")
    )
  }
  
  list(
    EBOut = eb_out,
    EBDERes = de_res,
    GeneFC = gene_fc,
    sig_table = sig_table,
    size_factors = size_factors
  )
}



############################
# 4. APPLY MAD-BASED FEATURE FILTERING & CALCULATE SIZE FACTORS
############################
size_factors <- MedianNorm(Rcounts[,-1])#Calculate on full, unfiltered dataset

#MAD filter full dataset
flt_all <- mad_filter_counts(Rcounts[,-1])
dim(flt_all$counts)#10978


############################
# 5. DEFINE COMPARISON SUBSETS
############################
cmp_1_4 <- subset_comparison(sev_meta, flt_all$counts, "Cluster_1", "Cluster_4")
dim(cmp_1_4$counts)#10978    15


############################
# 6. PAIRWISE CLUSTER COMPARISONS
############################

fit_1_4 <- run_ebseq_comparison(
  counts_mat = cmp_1_4$counts,
  conditions = cmp_1_4$meta$Group,
  size_factors=size_factors[rownames(cmp_1_4$meta)],
  annotation = gene,
  out_prefix = "cluster_1_vs_4TEST",
  fdr = 0.1,
  threshold_fc = 0.8,
  maxround = 10
)


#MAD filter full dataset
all_mad <- apply(Rcounts[,-1],1,mad) 
quantile(all_mad)
# 0%         25%         50%         75%        100% 
# 0.0000      0.0000      0.0000      2.9652 117340.3770 
keep <- names(which(all_mad>quantile(all_mad)[4]))
all_F <- Rcounts[keep,-1]
dim(all_F)#10978

EBOut=EBTest(Data=as.matrix(all_F[,rownames(cmp_1_4$meta)]),Conditions = cmp_1_4$meta$Group,
             sizeFactors=size_factors[rownames(cmp_1_4$meta)],
             maxround=10)

#1.Calculate MAD across full dataset not subset and filter
#2.Calculate size factors on full dataset (unfiltered)




# EBOut=EBTest(Data=as.matrix(Rcounts_one_four),Conditions = one_four$Group,sizeFactors=size_one_four, maxround=10)
EBOut=EBTest(Data=as.matrix(one_four_F),Conditions = one_four$Group,sizeFactors=size_one_four, maxround=10)

clean_data <- cmp_1_4$counts[complete.cases(cmp_1_4$counts), ]
dim(clean_data)
fit_1_4 <- run_ebseq_comparison(
  counts_mat = clean_data,
  conditions = cmp_1_4$meta$Group,
  annotation = gene,
  out_prefix = "cluster_1_vs_4_noMAD",
  fdr = 0.1,
  threshold_fc = 0.8,
  maxround = 10
)

# Removing transcripts with 100 th quantile < = 0 
# 9033 transcripts will be tested


############################
# 11. HEATMAPS OF SIGNIFICANT GENES
############################
logCPM_all <- cpm(as.matrix(Rcounts[, -1]), prior.count = 2, log = TRUE)

# Cluster 1 vs 4 heatmaps.
if (!is.null(fit_1_4$sig_table) && nrow(fit_1_4$sig_table) > 0) {
  counts_sig <- logCPM_all[rownames(fit_1_4$sig_table), colnames(cmp_1_4$counts), drop = FALSE]
  rownames(counts_sig) <- fit_1_4$sig_table$external_gene_name
  t_track <- cmp_1_4$meta
  
  pdf("cluster_1_vs_4_heatmap_row_scaled.pdf")
  pheatmap(counts_sig, annotation_col = t_track, scale = "row")
  dev.off()
  
  pdf("cluster_1_vs_4_heatmap_manual_sort.pdf")
  pheatmap(counts_sig, annotation_col = t_track, Colv = NA)
  dev.off()
}

###############################################################################
# END OF SCRIPT
###############################################################################

