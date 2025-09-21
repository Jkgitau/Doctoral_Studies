

# Set working directory ---------------------------------------------------

# setwd("~/Desktop/trial/clean_analyses/Final_6D")


# Load required libraries
library(ggplot2)
library(dplyr)
library(ggsignif)
library(ggpubr)

# Read your files
dREG <- read.table("filtered_gene_expression_in_Pore_C 2.tsv", header = TRUE, sep = "\t")
ATAC <- read.table("ATAC_in_gene-no_gene.tsv", header = TRUE, sep = "\t")
NONE <- read.table('matched_genes.txt', header = TRUE, sep = '\t')

# --- Step 1: Prepare input data
dREG$group <- "dREG"
ATAC$group <- "ATAC"
NONE$group <- "NONE"

# --- Combine into one dataframe
combined <- bind_rows(
  dREG %>% select(gene.log2FoldChange, group),
  ATAC %>% select(gene.log2FoldChange, group),
  NONE %>% select(gene.log2FoldChange, group)
)

# --- Remove rows with NA if any
combined <- combined %>% filter(!is.na(gene.log2FoldChange))

# --- Define all pairwise comparisons
comparisons <- list(
  c("dREG", "ATAC"),
  c("dREG", "NONE"),
  c("ATAC", "NONE")
)

# --- Boxplot with Bonferroni-corrected Mann–Whitney U tests
ggboxplot(combined, x = "group", y = "gene.log2FoldChange", 
          color = "group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          ylab = "gene.log2FoldChange", xlab = "Groups") +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",           # Mann–Whitney U test
    p.adjust.method = "bonferroni",  # Bonferroni correction
    label = "p.format"                # Show the actual p-value
  )
