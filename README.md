# RNA-Seq Differential Expression Analysis

This repository contains R scripts for analyzing transcriptomic data using `limma` and related Bioconductor packages.

## Dataset
* **Accession:** GSE313717
* **Platform:** RNA-Sequencing 

## Contents
* `GSE313717_Limma_Code.R`: Complete R script for data pre-processing, normalization, linear modeling, and differential expression evaluation.

## Prerequisites & Dependencies
To run the code, ensure you have R installed along with the following packages:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "limma",
  "edgeR",
  "RColorBrewer",
  "dplyr",
  "umap",
  "ggplot2",
  "tidyverse",
  "ggrepel",
  "pheatmap",
  "biomaRt",
  "VennDiagram",
  "matrixStats",
  "sparseMatrixStats",
  "DelayedMatrixStats",
  "readxl",
  "EnhancedVolcano",
  "openxlsx"
))
