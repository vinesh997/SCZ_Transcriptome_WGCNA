library(DESeq2)
library(tibble)
setwd(" ")

#Reading the phenodata
Phenodata<- read.csv("SCZ_CNT_Sample.csv")

#Filtering samples to remove the eight samples with low sequencing depth
Data<- read.delim("raw_counts_2.txt", check.names = F, sep = "") #dim 59050;38
Data<- t(Data)
Data<- as.data.frame(Data)
Data<- rownames_to_column(Data, var = "Sample")
Data<- Data[Data$Sample %in% Phenodata$Sample,]
Data<- Data %>% remove_rownames() %>% column_to_rownames(var = "Sample")
Data<- t(Data) #59050;29

write.table(Data, "Raw_counts.txt") ##use this file for further analysis as well as for WGCNA

#filtering the genes
smallestGroupSize <- 13
keep <- rowSums(Data >= 10) >= smallestGroupSize
data_filt<- Data[keep,]


###Adjusting all the covariates 

# Convert numeric variables to factors if they represent categories
Phenodata$Batch <- factor(Phenodata$Batch)
Phenodata$Sex <- factor(Phenodata$Sex)
Phenodata$Diagnosis <- factor(Phenodata$Diagnosis)

# List of variables to divide by median
variables <- c("Age", "NEUTROPHILS", "LYMPHOCYTES", "MONOCYTES", "EOSINOPHILS", "BASOPHILS")

# Loop over each variable and create a new column based on median grouping
for (var in variables) {
  median_value <- median(Phenodata[[var]], na.rm = TRUE)
  group_name <- paste0(var, "_group")
  
  # Create new binary group column for each variable
  Phenodata[[group_name]] <- ifelse(Phenodata[[var]] > median_value, "A", "B")
}

# View updated Phenodata with new group columns
head(Phenodata)


# Convert columns to factors
Phenodata[c("Batch", "Sex", "Diagnosis", 
            "Age_group", 
            "NEUTROPHILS_group", 
            "LYMPHOCYTES_group", 
            "MONOCYTES_group" ,
            "EOSINOPHILS_group",
            "BASOPHILS_group")] <- lapply(Phenodata[c("Batch", "Sex", "Diagnosis", 
                                                      "Age_group", 
                                                      "NEUTROPHILS_group", 
                                                      "LYMPHOCYTES_group", 
                                                      "MONOCYTES_group" ,
                                                      "EOSINOPHILS_group",
                                                      "BASOPHILS_group")], factor)

# Check to confirm they have been converted
str(Phenodata)


# Creating DESeq dataset
dds <- DESeqDataSetFromMatrix(data_filt,
                              Phenodata,
                              design = ~ NEUTROPHILS_group + LYMPHOCYTES_group + MONOCYTES_group 
                              + BASOPHILS_group + EOSINOPHILS_group + Age_group + Sex + Batch + 
                                + Diagnosis)



dds$condition <- factor(dds$Diagnosis, levels = c("CNT","SCZ"))



# Differential gene expression
dds <- DESeq(dds)

########## PCA and MA plot for visualization #####################
vsdata <- varianceStabilizingTransformation(dds, blind=FALSE)

plotPCA(vsdata, intgroup="Diagnosis") 
plotPCA(vsdata, intgroup="Batch")

res <- results(dds)
res

res <- results(dds, name="Diagnosis_SCZ_vs_CNT")
res <- results(dds, contrast=c("Diagnosis","SCZ","CNT"))

resultsNames(dds)

DESeq2::plotMA(res)

resLFC <- lfcShrink(dds, coef="Diagnosis_SCZ_vs_CNT", type="apeglm")
resLFC

resOrdered <- res[order(res$pvalue),]
summary(res)

res05 <- results(dds, alpha=0.05)
summary(res05)

sum(res05$padj < 0.05, na.rm=TRUE)


# Volcano plot
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.0)  +  ggplot2::coord_cartesian(ylim=c(0, 6))

# Define log2 fold change threshold
log2FC_threshold <- 0  # Use 0 to identify any changes above or below

# Separate upregulated and downregulated genes based on log2 fold change
upregulated <- res[res$log2FoldChange > log2FC_threshold, ]
downregulated <- res[res$log2FoldChange < -log2FC_threshold, ]

# Sort by log2 fold change
upregulated <- upregulated[order(upregulated$log2FoldChange, decreasing = TRUE), ]
upregulated<- as.data.frame(upregulated)
upregulated<- rownames_to_column(upregulated, var = "Genes")
rownames(upregulated)<- NULL

downregulated <- downregulated[order(downregulated$log2FoldChange, decreasing = FALSE), ]
downregulated<- as.data.frame(downregulated)
downregulated<- rownames_to_column(downregulated, var = "Genes")
rownames(downregulated)<- NULL

# View the top 10 upregulated genes
top_upregulated <- head(upregulated, 10)
top_upregulated

# View the top 10 downregulated genes
top_downregulated <- head(downregulated, 10)
top_downregulated

# Save top upregulated genes
upregulated$Genes<- as.character(upregulated$Genes)
write.csv(upregulated, file = "upregulated_genes_18_03_2025.csv", row.names = F, quote = T)

# Save top downregulated genes
downregulated$Genes<- as.character(downregulated$Genes)
write.csv(downregulated, file = "downregulated_genes_18_03_2025.csv", row.names = F, quote = T)

# Visualize gene expression between the groups
plotCounts(dds, gene=" ", intgroup="Diagnosis")
