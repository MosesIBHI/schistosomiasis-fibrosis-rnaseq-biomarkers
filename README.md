# RNA-seq Biomarker Discovery in Schistosoma mansoni-Associated Hepatic Fibrosis
This repository contains the scripts, workflows and supporting files used for the transcriptomic analyses performed in the study investigating biomarkers associated with Schistosoma mansoni-induced hepatic fibrosis in children from Cameroon.
The repository was developed to ensure full reproducibility of the differential expression analyses described in the manuscript.

# Study design
A total of 80 participants were recruited and classified into four epidemiological clusters:

- KK+US+ : Schistosoma-positive with ultrasound-detected hepatic fibrosis
- KK+US− : Schistosoma-positive without fibrosis
- KK−US+ : Schistosoma-negative with fibrosis
- KK−US− : Schistosoma-negative without fibrosis

To reduce sequencing costs and inter-individual variability, samples were pooled two-by-two within each cluster, resulting in 40 RNA-seq libraries.

The primary comparison investigated in this repository is:
KK+US+ versus KK−US−

# Sequencing and preprocessing

RNA libraries were prepared following ribosomal RNA depletion and sequenced using Illumina paired-end 150 bp chemistry.

Reads were processed using:

- FastQC for quality assessment
- STAR for alignment to the human reference genome (GRCh38)
- SAMtools for BAM processing and quality control
- RSEM for transcript quantification

Differential expression analyses were performed independently using DESeq2, edgeR and EBSeq.

# Differential expression workflows

## DESeq2

Institution:
IBHI Laboratory, Cameroon

Comparison:
KK+US+ (n=10) versus KK−US− (n=10)

Thresholds:
- FDR ≤ 0.20
- |log2FC| ≥ 1

Rationale:
Given the exploratory nature of the study and the relatively small sample size, a permissive FDR threshold was intentionally used to minimize false-negative discoveries during candidate biomarker screening.


## edgeR

Institution:
CBIO Unit, University of Cape Town, South Africa

Comparison:
KK+US+ (n=10) versus KK−US− (n=10)

Thresholds:
- FDR ≤ 0.05
- |log2FC| ≥ 1


## EBSeq

Institution:
CBIO Unit, University of Cape Town, South Africa

Comparison:
KK+US+ (n=10) versus KK−US− (n=5)

Thresholds:
- FDR ≤ 0.05
- |log2FC| ≥ 1

Five KK−US− control libraries were excluded following quality-control assessment because of globally low expression profiles that negatively affected model fitting and downstream inference.


# Cross-method comparison

Differentially expressed genes identified by DESeq2, edgeR and EBSeq were compared using InteractiveVenn.

Genes detected by all three methods were considered the most robust candidate biomarkers and were prioritized for downstream validation and biological interpretation.

# Candidate biomarkers
The recurrent candidate biomarkers identified across analyses included:

- INADL
- TCP11L1
- AMPD3
- SEMA7A
- DUS2








