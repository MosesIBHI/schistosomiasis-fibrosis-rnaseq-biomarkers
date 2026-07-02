###############################################################################
# DESCRIPTION:
# This script performs differential gene expression (DGE) analysis using edgeR.
# No filtering of samples or genes is applied prior
# to modelling, to preserve biological variation.
###############################################################################

###############################
# 1. SETUP & LIBRARIES
###############################

setwd("xxxxxxxxxxx")
# Core packages
library(edgeR)# version '4.0.16'
library(dplyr)# version '1.1.4'
library(tibble)# version '3.2.1'
library(magrittr)# version '2.0.3'
library(pheatmap)# version '1.0.12'
# Plotting / utilities
library(gplots)# version '3.1.3.1'

###############################
# 2. LOAD DATA
###############################

# Load count matrix
Rcounts <- read.csv("merged_gene_counts_2.txt",
                    sep = "\t",
                    row.names = 1)
dim(Rcounts)#63677    41
# Clean sample names
colnames(Rcounts) <- sub(".fwAligned.sortedByCoord.out.bam", "", colnames(Rcounts))
colnames(Rcounts) <- sapply(strsplit(colnames(Rcounts), "_S"), `[`, 1)

# Remove genes with zero counts across all samples
Rcounts <- Rcounts[rowSums(Rcounts[,-1]) > 0, ]
dim(Rcounts)#53541    41
###############################
# 3. LOAD METADATA
###############################

sev_meta <- read.csv2("Metadata.csv",
                      sep = ",",
                      row.names = 1)
# Ensure sample alignment
stopifnot(all(rownames(sev_meta) %in% colnames(Rcounts)))

# Reorder counts to match metadata
Rcounts <- data.frame(Rcounts[,1],Rcounts[,rownames(sev_meta)])
colnames(Rcounts)[1] <- "gene"
identical(rownames(sev_meta),colnames(Rcounts)[-1])#TRUE


###############################
# 4. CREATE DGE OBJECT
###############################

dge <- DGEList(
  counts = as.matrix(Rcounts[,-1]),
  genes = Rcounts[,1],
  group = factor(sev_meta$Group)
)

###############################
# 5. EXPLORATORY QC
###############################

pdf("MDS_all_samples.pdf")
plotMDS(dge, col = as.numeric(sev_meta$Group))
legend("topright", legend = levels(sev_meta$Group), col = 1:4, pch = 16)
dev.off()

# Library size diagnostics
lib_sizes <- colSums(dge$counts)
boxplot(lib_sizes, main = "Library sizes")
median(lib_sizes)#

###############################
# 6. DESIGN MATRIX
###############################

design <- model.matrix(~0 + Group, data = sev_meta)
colnames(design) <- gsub("Group", "", colnames(design))

###############################
# 7. NORMALISATION & DISPERSION
###############################

dge <- calcNormFactors(dge, log=TRUE)

dge <- estimateDisp(dge, design, robust = TRUE)

# Visualise biological coefficient of variation
plotBCV(dge)

###############################
# 8. MODEL FITTING (QL FRAMEWORK)
###############################

fit <- glmQLFit(dge, design, robust = TRUE)

###############################
# 9. DIFFERENTIAL EXPRESSION TESTING
###############################
gene <- readRDS("../ensembl_gchr37_annotation.rds")#to reload

contrast <- makeContrasts(Cluster_1-Cluster_4, levels=design)
res <- glmQLFTest(fit, contrast=contrast)
TT_cluster_1_4 <- topTags(res, n = Inf)
save <- TT_cluster_1_4$table[TT_cluster_1_4$table$FDR <= 0.15,]
save$annotation <- gene[rownames(save) ,"description"]
write.csv(save, file = "Unfiltered_input_cluster1_vs_4_FDR_0_15.csv")

###############################
# 10. VISUALISATION (HEATMAP)
###############################

# LogCPM transformation
logCPM <- cpm(dge, prior.count = 2, log = TRUE)

# Extract significant genes (C1 vs C4 example)
sig_genes <- TT_cluster_1_4$table[TT_cluster_1_4$table$FDR <= 0.15, ]
samples <- rownames(sev_meta[sev_meta$Group %in% c("Cluster_1", "Cluster_4"), ])

heat_data <- logCPM[rownames(sig_genes), samples]

# Replace rownames with gene symbols
rownames(heat_data) <- gene[rownames(heat_data), "external_gene_name"]

# Annotation track
ann <- sev_meta[samples, ]

pdf("heatmap_C1_vs_C4_FDR_015.pdf", 15, 10)
pheatmap(heat_data, annotation_col = ann)
dev.off()

###############################################################################
# END OF SCRIPT
###############################################################################