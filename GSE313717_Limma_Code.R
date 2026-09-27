###############################################################################################################################################################
# Title: 	Gene expression signature of liver of mice with and without functional NNT

# Overall design: At the age of 11-13 weeks, half of the mouse cohort (female as well 
# as male mice) were switched from chow to a 45 kcal% Fat (ssniff E15744-347, according 
# to research diets D12451). After 11 weeks of HFD feeding, mice were sacrificed 2–4,5 h 
# after the lights went on, and serum and tissue samples were collected.


# Platform: GPL24247	Illumina NovaSeq 6000 (Mus musculus)


# Call required package
library(limma)
library(edgeR)
library(RColorBrewer)
library(dplyr)
library(umap)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(biomaRt)
library(VennDiagram)
library(matrixStats)
library(sparseMatrixStats)
library(DelayedMatrixStats)
library(readxl)
library(EnhancedVolcano)
library(openxlsx)

# Read in the raw counts matrix 
## ***GSE313717 datasets*** I got the 5 different datasets (depending on the diet comparison)
# Set file path
fp <- "C:/Users/dnnamani/OneDrive - University of Tennessee/UTHSC/Personal_Ideas/"

# get to the data folder
dfp <- paste0(fp, "Diabetic_Liver_Data/") 

# Set the NCD&HFD file path
HFDfp <- paste0(dfp, "GSE313717_count.xlsx")

#read in count data 
HFD = as.data.frame(read_xlsx(HFDfp)) 

# Connect to the official Ensembl Mouse dataset
ensembl <- useEnsembl(biomart = "genes",
                      dataset = "mmusculus_gene_ensembl",
                      mirror = "useast") 

# Query BioMart for the matching common gene symbols
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"), 
  filters    = "ensembl_gene_id",                
  values     = HFD$Geneid,                      
  mart       = ensembl
)

# Merge the new gene symbols back into the original dataframe
HFD_mapped <- merge(HFD, gene_map, 
                      by.x = "Geneid", 
                      by.y = "ensembl_gene_id", 
                      all.x = TRUE)

HFD <- HFD_mapped 

# Relocate the new column to the front so it's easy to see
HFD <- HFD %>% relocate(mgi_symbol, .after = Geneid)

# Remove the first column from the data
HFD$Geneid <- NULL

# Handle any duplicate gene names
HFD <- HFD %>%
  group_by(mgi_symbol) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  as.data.frame()

# Handle NA
HFD <- na.omit(HFD)

# Set row names to the first column
rownames(HFD) <- HFD$mgi_symbol

# Remove the symbol column column from the data
HFD$mgi_symbol <- NULL

# Confirm if the data was normalized: Massive difference between reads shows that it wasn't
colSums(HFD)

countsHFD <- HFD

# Set up metadata with group labels for each sample
# 5 replicates per group (In this case)
coldata <- data.frame(
  row.names = colnames(countsHFD),
  group = factor(c(
    rep("Male_Chow", 5),
    rep("Male_HFD", 5),
    rep("Female_Chow", 5),
    rep("Female_HFD", 5)
  ))
)

# Create DGEList object
dge <- DGEList(counts = countsHFD, group = coldata$group)
dge_raw     <- dge
lcpm_raw    <- cpm(dge_raw, log = TRUE)
L <- mean(dge_raw$samples$lib.size) * 1e-6
M <- median(dge_raw$samples$lib.size) * 1e-6
lcpm.cutoff <- log2(10/M + 2/L)

#  Filter low-count genes 
keep <- filterByExpr(dge, group = coldata$group)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalize counts
dge <- calcNormFactors(dge)
lcpm_filt <- cpm(dge, log = TRUE)

cat("Genes retained after filtering: ", nrow(dge), "\n", sep = "")

########### Design matrix##########
design <- model.matrix(~ 0 + group, data = coldata)
colnames(design) <- levels(coldata$group)

########### Define multiple contrasts##########
contrast.matrix <- makeContrasts(
  Male_HFD_vs_Male_Chow = Male_HFD - Male_Chow,
  Female_HFD_vs_Female_Chow = Female_HFD - Female_Chow,
  Male_HFD_vs_Female_HFD = Male_HFD - Female_HFD,
  Male_Chow_vs_Female_Chow = Male_Chow - Female_Chow,
  levels = design
)

## Estimate dispersion. Mean-variance modelling
v <- voom(dge, design, plot = TRUE)

# Fit model
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

plotSA(fit2, main = "Final model: Mean-variance trend")

############## density plots ################################################
nsamples    <- ncol(dge_raw)
samplenames <- colnames(dge_raw)
col <- colorRampPalette(brewer.pal(12, "Paired"))(nsamples)  

par(mfrow = c(1, 2))
plot(density(lcpm_raw[, 1]), col = col[1], lwd = 2, ylim = c(0, 0.26), las = 2,
     main = "", xlab = "")
title(main = "A. Raw data", xlab = "Log-cpm")
abline(v = lcpm.cutoff, lty = 3)
for (i in 2:nsamples) lines(density(lcpm_raw[, i])$x, density(lcpm_raw[, i])$y,
                            col = col[i], lwd = 2)
legend("topright", samplenames, text.col = col, bty = "n", cex = 0.6)

plot(density(lcpm_filt[, 1]), col = col[1], lwd = 2, ylim = c(0, 0.26), las = 2,
     main = "", xlab = "")
title(main = "B. Filtered data", xlab = "Log-cpm")
abline(v = lcpm.cutoff, lty = 3)
for (i in 2:nsamples) lines(density(lcpm_filt[, i])$x, density(lcpm_filt[, i])$y,
                            col = col[i], lwd = 2)
legend("topright", samplenames, text.col = col, bty = "n", cex = 0.6)
par(mfrow = c(1, 1))
###############################################################################

# Perform LRT and Extract Results for Each Comparison

for (cn in colnames(contrast.matrix)) {
  res <- topTable(fit2, coef = cn, number = Inf, sort.by = "P")
  res$SE <- res$logFC / res$t
  names(res)[names(res) == "P.Value"]   <- "PValue"
  names(res)[names(res) == "adj.P.Val"] <- "FDR"
  
  cat(cn, ": ", sum(res$FDR < 0.05 & abs(res$logFC) > 1), " DEGs\n", sep = "")
  write.csv(res, file = paste0(dfp, "limma_DEG_", cn, ".csv"))
}

# UMAP plot (instead of pca)
# Create UMAP plot from expression data
# Get log2 CPM (counts per million)
logcpm <- cpm(dge, log = TRUE)

# Transpose the logcpm matrix to samples x genes
umap_input <- t(logcpm)

# Run UMAP
set.seed(42)  # For reproducibility
umap_res <- umap(umap_input, n_neighbors = 10, min_dist = 0.3, metric = "euclidean")

# Format Coordinates for Plotting
# Turn result into data frame
umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# Pull sample IDs from the layout rows so we can map metadata
umap_df$SampleID <- rownames(umap_input)

# Add metadata
umap_df$group <- coldata$group

# Plot
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = group, label= SampleID)) +
  geom_point(size = 4) +
  #geom_text(vjust = 1, size = 3) +
  #stat_ellipse(type = "norm", linetype = 2, size = 1, alpha = .2) + 
  theme_minimal() +
  scale_color_manual(values = c("Male_Chow" = "cornflowerblue",
                                "Male_HFD" = "indianred1",
                                "Female_Chow" = "lightgreen",
                                "Female_HFD" = "violet")) +
  labs(title = "UMAP of Normalized Expression Data : GSE313717", x= "UMAP1", y="UMAP2") +
  theme(text = element_text(size = 14))


ggsave(file = paste0(dfp, "UMAP_GSE313717.png"), width = 8, height = 6, dpi = 3000)


#############################
###############################################################################
logcpm <- cpm(dge, log = TRUE)
pca <- prcomp(t(logcpm), scale. = FALSE)
pv  <- round(100 * summary(pca)$importance[2, 1:2], 1)

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     SampleID = colnames(logcpm),
                     group = coldata$group,
                     lib_size = dge$samples$lib.size / 1e6)

ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(aes(label = SampleID), size = 2.5, show.legend = FALSE) +
  labs(title = "PCA: GSE313717",
       x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)")) +
  theme_minimal(base_size = 14)
##########################################################

# I have already log transform (when I did the umap) so
# Plot heatmap of top 100 genes (all treatment)
# Ensure filtered_data is treated as a matrix
cn  <- "Male_HFD_vs_Male_Chow"
res <- topTable(fit2, coef = cn, number = Inf, sort.by = "P")
top_degs <- head(rownames(res)[res$adj.P.Val < 0.05 & abs(res$logFC) > 1], 100)

# all 12 samples, ordered by group 
grp_order <- order(factor(coldata$group,
                          levels = c("Male_Chow","Male_HFD", "Female_Chow", "Female_HFD")))
mat <- as.matrix(logcpm)[top_degs, grp_order, drop = FALSE]

ann <- data.frame(
  group    = coldata$group[grp_order],
  row.names = colnames(mat)
)

pheatmap(mat,
         annotation_col = ann,
         scale = "row",
         cluster_rows = TRUE, 
         cluster_cols = FALSE,   # FALSE: keep group order
         treeheight_row = 0,
         show_rownames = TRUE, show_colnames = FALSE,
         color = colorRampPalette(c("blue","white","red"))(100),
         fontsize_row = 10,
         main = paste0("Top 100 DEGs across all groups: GSE313717"))

################################
contrast_names <- colnames(contrast.matrix)
titles <- c("Effect of High fat Diet (Male_HFD_vs_Male_Chow)",
            "Female_HFD_vs_Female_Chow",
            "Effect of HFD Diet (Male_vs_Female)",
            "Male_Chow_vs_Female_Chow")

deg_list <- list()
for (i in seq_along(contrast_names)) {
  cn  <- contrast_names[i]
  res <- topTable(fit2, coef = cn, number = Inf, sort.by = "none")
  names(res)[names(res) == "P.Value"]   <- "PValue"
  names(res)[names(res) == "adj.P.Val"] <- "FDR"
  
  print(EnhancedVolcano(res, lab = rownames(res),
                        x = 'logFC', y = 'FDR',
                        pCutoff = 0.05, FCcutoff = 1.0,
                        title = paste0("Volcano: ", titles[i]),
                        subtitle = NULL, pointSize = 3, labSize = 3, colAlpha = 0.8))
  
  deg <- res[res$FDR < 0.05 & abs(res$logFC) >= 1, ]
  deg_list[[cn]] <- deg
  write.csv(deg, file = paste0(dfp, "significant_limma_DEGs_", cn, ".csv"))
  cat(cn, ": ", nrow(deg), " DEGs\n", sep = "")
}


