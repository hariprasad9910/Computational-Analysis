# Gene Expression Analysis Pipeline

This repository contains an R script for analyzing gene expression data from the Gene Expression Omnibus (GEO) database, specifically focused on differential expression analysis between normal and tumor samples.

## Overview

The pipeline performs the following steps:
1. Data retrieval from GEO
2. Data preprocessing and normalization
3. Quality control and exploratory data analysis 
4. Differential expression analysis
5. Visualization of results
6. Gene annotation and pathway analysis

## Prerequisites

The following R packages are required:

```r
# Core analysis packages
install.packages(c("ggplot2", "dplyr", "readr", "ggrepel", "pheatmap"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "org.Hs.eg.db", "impute"))
```

For MetaDE package installation:
```r
# Option 1: From archive
# download.packages('MetaDE_1.0.5.tar.gz', destdir=".", repos=NULL)
# install.packages('MetaDE_1.0.5.tar.gz', repos=NULL, type='source')

# Option 2: From GitHub
# install.packages("devtools")
# devtools::install_github("MetaOmics/MetaDE")
```

## Script

```r
# Gene Expression Analysis Pipeline
# This script analyzes GEO dataset GSE33127 to identify differentially expressed genes
# between normal and tumor samples

# ---------- Load Libraries ----------
library(GEOquery)
library(MetaDE)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(readr)
library(limma)
library(ggrepel)
library(org.Hs.eg.db)

# ---------- Data Retrieval ----------
# Set working directory if needed
# setwd("your/project/directory")

# Download GEO dataset (uses cache if already downloaded)
data <- getGEO("GSE33127", destdir = "data/")
gse <- data[[1]]  # Extract the ExpressionSet object

# ---------- Data Exploration ----------
# Inspect the dataset
dim(gse)  # Number of genes and samples
# Print sample information
sampleInfo <- pData(gse)
# Extract relevant sample information
sampleInfo <- select(sampleInfo, source_name_ch1, characteristics_ch1.1)
# Rename columns for clarity
sampleInfo <- rename(sampleInfo, group = source_name_ch1, patient = characteristics_ch1.1)

# ---------- Data Preprocessing ----------
# Log2 transform expression values if they're not already
if(max(exprs(gse), na.rm=TRUE) > 100) {
  message("Applying log2 transformation to expression data")
  exprs(gse) <- log2(exprs(gse))
}

# ---------- Quality Control ----------
# Create boxplot of expression values
pdf("figures/expression_boxplot.pdf")
boxplot(exprs(gse), outline=FALSE, main="Expression Distribution", 
        xlab="Samples", ylab="log2 Expression")
dev.off()

# Generate correlation heatmap
corMatrix <- cor(exprs(gse), use="complete.obs")
# Make sure rownames match column names for annotations
rownames(sampleInfo) <- colnames(corMatrix)

pdf("figures/correlation_heatmap.pdf", width=10, height=8)
pheatmap(corMatrix, annotation_col=sampleInfo)
dev.off()

# Principal Component Analysis
pca <- prcomp(t(exprs(gse)))

pdf("figures/pca_plot.pdf")
# Create PCA plot with sample labels
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group, label=paste("Patient", patient))) + 
  geom_point(size=3) +
  geom_text_repel() +
  theme_bw() +
  labs(title="PCA of Gene Expression Data",
       x=paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y=paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"))
dev.off()

# ---------- Differential Expression Analysis ----------
# Create design matrix for the linear model
design <- model.matrix(~0+sampleInfo$group)
colnames(design) <- c("Normal", "Tumour")

# Filter low-expressed genes
# Calculate median expression level
cutoff <- median(exprs(gse))
is_expressed <- exprs(gse) > cutoff
# Keep genes expressed in more than 2 samples
keep <- rowSums(is_expressed) > 2
# Subset to expressed genes
gse_filtered <- gse[keep,]

# Calculate array weights to account for variable quality
aw <- arrayWeights(exprs(gse_filtered), design)

# Fit linear model with weights
fit <- lmFit(exprs(gse_filtered), design, weights = aw)

# Define contrasts
contrasts <- makeContrasts(Tumour - Normal, levels=design)
fit2 <- contrasts.fit(fit, contrasts)

# Apply empirical Bayes method
fit2 <- eBayes(fit2)

# Prepare gene annotations
anno <- fData(gse_filtered)
anno <- select(anno, Symbol, Entrez_Gene_ID, Chromosome, Cytoband)
fit2$genes <- anno

# Extract full results table
full_results <- topTable(fit2, number=Inf)
full_results <- tibble::rownames_to_column(full_results, "ID")

# Save complete results
write_csv(full_results, "results/full_differential_expression_results.csv")

# ---------- Results Visualization ----------
# Set cutoffs for significant genes
p_cutoff <- 0.05
fc_cutoff <- 1  # log2 fold change (corresponds to 2-fold)

# Filter significant genes
significant_genes <- filter(full_results, adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff)
write_csv(significant_genes, "results/significant_differentially_expressed_genes.csv")

# Create volcano plot
pdf("figures/volcano_plot.pdf")
topN <- 20  # Number of top genes to label

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff & abs(logFC) > fc_cutoff) %>% 
  mutate(Rank = rank(adj.P.Val), Label = ifelse(Rank <= topN, Symbol, "")) %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), col=Significant, label=Label)) + 
  geom_point(alpha=0.7) + 
  geom_text_repel(col="black", max.overlaps=20) +
  geom_vline(xintercept=c(-fc_cutoff, fc_cutoff), linetype="dashed", color="gray") +
  geom_hline(yintercept=-log10(p_cutoff), linetype="dashed", color="gray") +
  scale_color_manual(values=c("FALSE"="gray", "TRUE"="red")) +
  theme_bw() +
  labs(title="Volcano Plot of Differential Expression",
       x="log2 Fold Change (Tumor vs Normal)",
       y="-log10 Adjusted P-value")
dev.off()

# Create heatmap of top differentially expressed genes
topN <- 20
ids_of_interest <- mutate(full_results, Rank = rank(adj.P.Val)) %>% 
  filter(Rank <= topN) %>% 
  pull(ID)
gene_names <- full_results %>% 
  filter(ID %in% ids_of_interest) %>% 
  pull(Symbol)

gene_matrix <- exprs(gse_filtered)[ids_of_interest,]

pdf("figures/top_genes_heatmap.pdf", width=10, height=8)
pheatmap(gene_matrix,
         labels_row = gene_names,
         scale="row",
         annotation_col = sampleInfo,
         main="Top Differentially Expressed Genes",
         show_colnames = FALSE)
dev.off()

# ---------- Gene Set Analysis ----------
# Example: Analyze specific genes of interest
my_genes <- c("KIAA1199", "CDH3", "CLDN1", "GUCA2A")

# Get gene IDs and extract expression data
ids_of_interest <- filter(full_results, Symbol %in% my_genes) %>% pull(ID)
gene_names <- filter(full_results, Symbol %in% my_genes) %>% pull(Symbol)

if(length(ids_of_interest) > 0) {
  gene_matrix <- exprs(gse_filtered)[ids_of_interest,]
  
  pdf("figures/genes_of_interest_heatmap.pdf", width=8, height=6)
  pheatmap(gene_matrix,
           labels_row = gene_names,
           scale="row",
           annotation_col = sampleInfo,
           main="Expression of Genes of Interest")
  dev.off()
}

# Example: Analyze genes in a specific GO term (chromatin remodeling)
go_genes <- AnnotationDbi::select(org.Hs.eg.db,
                                 columns="SYMBOL",
                                 keys="GO:0006338",  # Chromatin remodeling
                                 keytype="GO")

if(nrow(go_genes) > 0) {
  go_ids <- filter(full_results, Symbol %in% go_genes$SYMBOL) %>% pull(ID)
  go_gene_names <- filter(full_results, Symbol %in% go_genes$SYMBOL) %>% pull(Symbol)
  
  if(length(go_ids) > 0) {
    go_matrix <- exprs(gse_filtered)[go_ids,]
    
    pdf("figures/go_term_genes_heatmap.pdf", width=10, height=8)
    pheatmap(go_matrix,
             labels_row = go_gene_names,
             scale="row",
             annotation_col = sampleInfo,
             main="Expression of Chromatin Remodeling Genes (GO:0006338)")
    dev.off()
  }
}
```

## Usage

1. Clone this repository
2. Install the required R packages
3. Run the script in R or RStudio
4. Results will be saved in the 'results' directory and figures in the 'figures' directory

## Script Workflow

### Data Retrieval
The script downloads data from GEO dataset GSE33127, which contains gene expression data comparing tumor and normal samples.

### Preprocessing
- Log2 transformation of expression values
- Quality control using boxplots, correlation heatmaps, and PCA
- Filtering of lowly expressed genes

### Differential Expression Analysis
- Uses limma package to identify differentially expressed genes
- Applies empirical Bayes method for robust statistics
- Accounts for array quality using weights

### Results Visualization
- Volcano plot showing significantly up/down-regulated genes
- Heatmaps of top differentially expressed genes
- Visualization of specific genes of interest

### Gene Set Analysis
- Analysis of user-defined genes of interest
- GO term-based gene set analysis

## Output Files

The script generates the following outputs:
- `full_differential_expression_results.csv`: Complete results for all genes
- `significant_differentially_expressed_genes.csv`: Filtered list of significant genes
- Various visualizations in the figures directory


## Acknowledgments

This analysis pipeline is based on techniques from various bioinformatics resources and utilizes the publicly available GEO dataset GSE33127.