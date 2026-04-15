#Loading required libraries for differential expression and plotting
library(DESeq2)
library(tibble)
library(ggbreak)
library(EnhancedVolcano)

#Setting working directory
setwd(" ")

#Reading the phenodata
Phenodata<- read.csv("") # 47;29

#Reading the raw count matrix
Data<- read.delim(" ", check.names = F, sep = "")

#Identifying common samples between count data and phenodata
common_samples <- intersect(colnames(Data), Phenodata$Sample.ID)

#Subsetting count matrix and phenodata to common samples
Data <- Data[, common_samples] # 59050;47
Phenodata <- Phenodata[Phenodata$Sample.ID %in% common_samples, ]

#Reordering count matrix columns to match phenodata sample order
Data <- Data[, match(Phenodata$Sample.ID, colnames(Data))]

#Ensuring sample order is identical between both datasets
stopifnot(all(colnames(Data) == Phenodata$Sample.ID))

#Converting diagnosis column to factor for DESeq2 design
Phenodata$Diagnosis <- factor(Phenodata$Diagnosis)

#Creating DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(Data,
                              Phenodata,
                              design = ~ Diagnosis)

# Create DESeq
#Running differential expression analysis pipeline
dds2 <- DESeq(dds)

#filtering the genes
#Filtering low-count genes based on minimum count threshold across samples
smallestGroupSize <- 22
keep <- rowSums(counts(dds2) >= 30) >= smallestGroupSize
dds2<- dds2[keep,]

#Variance stabilizing transformation for downstream visualization
vsd <- vst(dds2, blind=TRUE)
mat <- assay(vsd)

#PCA plot to visualize sample clustering by diagnosis
plotPCA(vsd, intgroup=c("Diagnosis"))

# Differential expression analysis
#Checking available result names for contrasts
resultsNames(dds2)

#Extracting results for SCZ vs Control comparison
res_group_CTR_vs_SCZ <- results(object = dds2, name="Diagnosis_SCZ_vs_CNT", alpha = 0.05)

# Summary of DE analysis
#Overview of differentially expressed genes
summary(res_group_CTR_vs_SCZ)

# MA plot
#Visualizing log fold changes vs mean expression
#plotMA(res_group_CTR_vs_SCZ)
DESeq2::plotMA(res_group_CTR_vs_SCZ)

#Volcano plot for visualizing significance vs fold change
EnhancedVolcano(res_group_CTR_vs_SCZ,
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.0)  +  ggplot2::coord_cartesian(ylim=c(0, 10))

# Define log2 fold change threshold
log2FC_threshold <- 0  # Use 0 to identify any changes above or below

# Separate upregulated and downregulated genes based on log2 fold change
upregulated <- res_group_CTR_vs_SCZ[res_group_CTR_vs_SCZ$log2FoldChange > log2FC_threshold, ]
downregulated <- res_group_CTR_vs_SCZ[res_group_CTR_vs_SCZ$log2FoldChange < -log2FC_threshold, ]

# Sort by log2 fold change
#Sorting genes based on magnitude of differential expression
upregulated <- upregulated[order(upregulated$log2FoldChange, decreasing = TRUE), ]
upregulated<- as.data.frame(upregulated) # 6725;6
upregulated<- rownames_to_column(upregulated, var = "Genes")
rownames(upregulated)<- NULL

downregulated <- downregulated[order(downregulated$log2FoldChange, decreasing = FALSE), ]
downregulated<- as.data.frame(downregulated) # 7978;6
downregulated<- rownames_to_column(downregulated, var = "Genes")
rownames(downregulated)<- NULL

# View the top 10 upregulated genes
top_upregulated <- head(upregulated, 10)
top_upregulated

# View the top 10 downregulated genes
top_downregulated <- head(downregulated, 10)
top_downregulated

# Save top upregulated genes
#Converting gene names to character for export
upregulated$Genes<- as.character(upregulated$Genes)
#write.csv(upregulated, file = "upregulated_genes_16_03_2026.csv", row.names = F, quote = T)

#write.table(
#  as.character(upregulated$Genes),
#  file = "upregulated_genes_16_03_2026.csv",
#  row.names = FALSE,
#  col.names = FALSE,
#  quote = FALSE
#)

# Save top downregulated genes
#Converting gene names to character for export
downregulated$Genes<- as.character(downregulated$Genes)
#write.csv(downregulated, file = "downregulated_genes_16_03_2026.csv", row.names = F, quote = T)

#write.table(
#  as.character(downregulated$Genes),
#  file = "downregulated_genes_16_03_2026.csv",
#  row.names = FALSE,
#  col.names = FALSE,
#  quote = FALSE
#)

#Plot normalized counts for a specific gene across groups
plotCounts(dds2, gene="KRT5", intgroup="Diagnosis")

#Extract normalized counts from DESeq2 object
norm_counts <- counts(dds2, normalized=TRUE)

# Subset gene of interest
#Extract expression values for selected gene
gene_expr <- norm_counts["TADA2A", , drop = FALSE]

#Create dataframe for plotting
gene_df <- data.frame(
  Sample = colnames(gene_expr),
  Expression = as.numeric(gene_expr[1, ])
)

#Ensure sample IDs are in character format for merging
Phenodata$Sample.ID<- as.character(Phenodata$Sample.ID)

# Merge phenotype
#Merge gene expression data with phenodata
gene_df <- gene_df %>%
  dplyr::left_join(
    Phenodata,
    by = c("Sample" = "Sample.ID")   
  )

# Ensure factor order
#Setting order of diagnosis groups for plotting
gene_df$Diagnosis <- factor(gene_df$Diagnosis, levels = c("CNT", "SCZ"))

#Compute mean expression per group
gene_summary <- gene_df %>%
  dplyr::group_by(Diagnosis) %>%
  dplyr::summarise(
    mean_expr = mean(Expression, na.rm = TRUE)
  )

# ---- ADD SD HERE ----
#Compute mean and standard deviation for each group
gene_summary <- gene_df %>%
  dplyr::group_by(Diagnosis) %>%
  dplyr::summarise(
    mean_expr = mean(Expression, na.rm = TRUE),
    sd_expr   = sd(Expression, na.rm = TRUE)
  )

# ---- BARPLOT WITH SD ----
#Barplot showing mean expression with standard deviation
ggplot(gene_summary, aes(x = Diagnosis , y = mean_expr, fill = Diagnosis)) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(
      ymin = mean_expr - sd_expr,
      ymax = mean_expr + sd_expr
    ),
    width = 0.2
  ) +
  labs(
    title = "GPR18 mean expression",
    x = "Group",
    y = "Mean expression "
  ) +
  theme_minimal()

#If break is needed in the barplot
#Barplot with y-axis break for better visualization of large differences
ggplot(gene_summary, aes(x = Diagnosis , y = mean_expr, fill = Diagnosis)) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(
      ymin = mean_expr - sd_expr,
      ymax = mean_expr + sd_expr
    ),
    width = 0.1
  ) +
  scale_y_break(c(100,1200)) +   
  labs(
    title = "MAPKAP1 mean expression",
    x = "Group",
    y = "Mean expression"
  ) +
  theme_minimal()
