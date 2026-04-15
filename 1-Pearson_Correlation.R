#Loading the libraries
library(tibble)
library(DESeq2)
library(ggplot2)
library(ggfortify)
library(Hmisc)
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(grid)

#Setting working directory to correlation analysis folder
setwd(" ")

#Reading the count matrix
Data<- read.delim("", check.names = F, sep = "") # 59050;47

#filtering the genes
smallestGroupSize <- 22
keep <- rowSums(Data >= 30) >= smallestGroupSize
data_filt<- Data[keep,] # 14703;47

#reading the phenodata
Pheno<- read.csv("") # 47;29

#filtering SCZ samples from phenodata
#subseting based on diagnosis column in phenodata
Pheno_SCZ<- Pheno[Pheno$Diagnosis == "SCZ", ] # 25;29
#table(Pheno_SCZ$Batch) #checking the number of batches

#filtering only SCZ samples
#Transpose and reformat count data for sample-wise filtering
data_filt<- t(data_filt)
data_filt<- as.data.frame(data_filt)
data_filt<- rownames_to_column(data_filt, var = "Sample")

#Retain only SCZ samples present in phenodata
dat_sub<- data_filt[data_filt$Sample %in% Pheno_SCZ$Sample.ID,] # 25;14704

#Restore gene × sample format
dat_sub<- dat_sub %>% remove_rownames() %>% column_to_rownames(var = "Sample")
dat_sub<-t(dat_sub)

#Ensuring sample matching between phenodata and the data
Pheno_SCZ <- Pheno_SCZ[match(colnames(dat_sub), Pheno_SCZ$Sample.ID), ]
stopifnot(all(colnames(dat_sub) == Pheno_SCZ$Sample.ID))

#transforming the data using vst of DESeq2
vsd<- varianceStabilizingTransformation(dat_sub)

#visualizing for batch on PCA
Pheno_SCZ$Batch<- factor(Pheno_SCZ$Batch)
autoplot(prcomp(t(vsd), scale. = T), data = Pheno_SCZ, colour="Batch") + ggtitle("Before Batch Correction")

#Editing the dataframe for keeping samples as rows
expr <- t(vsd) %>% as.data.frame()  
expr <- rownames_to_column(expr, "Sample_ID")

# Match samples between expression data and phenodata
Pheno_SCZ <- Pheno_SCZ[Pheno_SCZ$Sample.ID %in% expr$Sample_ID,]
expr  <- expr[expr$Sample_ID %in% Pheno_SCZ$Sample.ID,]

#Reorder samples and convert back to gene × sample format
expr  <- expr %>% remove_rownames() %>% column_to_rownames(var = "Sample_ID")

#Editing phenodata for correlation analysis
#keeping only Age, Batch, Sex, DoD, Lymphocytes, Monocytes, P, N, G, T
Pheno_SUHRC_2 <- Pheno_SCZ %>% remove_rownames() %>% column_to_rownames("Sample.ID")
Pheno_SUHRC_2<- Pheno_SUHRC_2[,-c(1:2,4,10,13:14,18:23,25:28)] # 25;12

# final alignment
#Ensure both matrices contain the same samples in the same order
common <- intersect(rownames(expr), rownames(Pheno_SUHRC_2))
expr  <- expr[common, ]
Pheno_SUHRC_2 <- Pheno_SUHRC_2[common, ]

#Pearson Correlation: Gene~Traits
#Compute gene-wise correlations with all clinical traits
results_list <- list()

for (gene in colnames(expr)) {
  
  df <- cbind(GENE = expr[, gene], Pheno_SUHRC_2)
  rc <- rcorr(as.matrix(df), type = "pearson")
  
  results_list[[gene]] <- data.frame(
    Gene        = gene,
    PANSS       = colnames(Pheno_SUHRC_2),
    Correlation = rc$r["GENE", colnames(Pheno_SUHRC_2)],
    Pvalue      = rc$P["GENE", colnames(Pheno_SUHRC_2)]
  )
}

 #Combine all gene-wise correlation results
results <- bind_rows(results_list) # 176436;4

#Creating Correlation and P-value Matrices
#Reshape long-format results into wide matrices
corr_matrix_all <- dcast(results, Gene ~ PANSS, value.var = "Correlation") # 14703;13
p_matrix_all    <- dcast(results, Gene ~ PANSS, value.var = "Pvalue") # 14703;13

#Assign gene names as rownames
rownames(corr_matrix_all) <- corr_matrix_all$Gene
rownames(p_matrix_all)    <- p_matrix_all$Gene

#Remove redundant Gene column
corr_matrix_all <- corr_matrix_all[, -1]
p_matrix_all    <- p_matrix_all[, -1]

#identify genes significant for PANSS Negative syndrome

sig_neg <- results %>%
  dplyr::filter(
    PANSS == "N",
    Pvalue < 0.05
  ) %>%
  dplyr::pull(Gene) %>%
  unique()

#Subset correlation and p-value matrices (only genes significant for N)

filtered_corr <- corr_matrix_all[sig_neg, , drop = FALSE] # 2211;12
filtered_p    <- p_matrix_all[sig_neg, , drop = FALSE] # 2211;12

# Identify confounder/covariate-associated genes
#Traits considered as potential confounders

batch_trait <- c(
  "Lymphocytes",
  "Monocytes",
  "Batch_1",
  "Batch_2",
  "Batch_3", 
  "Age", 
  "Sex", 
  "DoD"
)

#Genes significantly associated with any confounder
sig_batch <- rownames(filtered_p)[
  rowSums(filtered_p[, batch_trait, drop = FALSE] < 0.05) > 0
]

# List of significant genes per PANSS negative trait
#    (PANSS traits only — confounder/covariate-associated excluded)

panss_traits <- c("N", "P", "G", "T")

sig_list <- lapply(panss_traits, function(tr) {
  rownames(filtered_p)[filtered_p[, tr] < 0.05]
})
names(sig_list) <- panss_traits

# Keep Negative genes but REMOVE confounder/covariate-associated genes

genes_negative <- setdiff(sig_list$N, sig_batch) # 141 genes

#cutoff of the correlated genes
#Apply correlation threshold to retain strongly associated genes
corr_cutoff <- 0.4

genes_negative <- genes_negative[
  abs(corr_matrix_all[genes_negative, "N"]) >= corr_cutoff
] # 133 genes

#Sort genes by strength of correlation with PANSS-N
neg_unique_values <- results %>%
  dplyr::filter(
    Gene %in% genes_negative,
    PANSS == "N"
  ) %>%
  dplyr::arrange(desc(abs(Correlation)))

# Sample size
n_samples <- 25

#power calculation
#Estimate statistical power for each gene-trait correlation
neg_unique_values$Power <- sapply(abs(neg_unique_values$Correlation), function(r) {
  pwr::pwr.r.test(
    r = r,
    n = n_samples,
    sig.level = 0.05,
    alternative = "two.sided"
  )$power
})

#calculating FDR
#Adjust p-values for multiple testing using Benjamini-Hochberg method
p_matrix_all$FDRneg<- p.adjust(p_matrix_all$N, method = "BH")
p_matrix_all<- rownames_to_column(p_matrix_all, var = "Gene")

#Merge FDR values with gene-level results
neg_unique_values <- neg_unique_values %>%
  dplyr::left_join(
    p_matrix_all %>% dplyr::select(Gene, FDRneg),
    by = "Gene"
  )

#filtering genes based on power
#Retain genes with sufficient statistical power
neg_power_filtered <- neg_unique_values %>%
  dplyr::filter(Power >= 0.8) # 10 genes

#Final list of high-confidence PANSS-N associated genes
genes_high_power_N <- neg_power_filtered$Gene

#Define PANSS domains for multi-trait comparison
traits_main <- c("N", "P", "G", "T")

#Extract correlations of selected genes across all PANSS domains
multi_trait_df <- results %>%
  dplyr::filter(
    Gene %in% genes_high_power_N,
    PANSS %in% traits_main
  )

#Compute power for each gene across all PANSS domains
multi_trait_df <- multi_trait_df %>%
  dplyr::mutate(
    Power_all = sapply(abs(Correlation), function(r) {
      pwr::pwr.r.test(
        r = r,
        n = n_samples,
        sig.level = 0.05,
        alternative = "two.sided"
      )$power
    })
  )

#Create a summary table with correlation, p-value, and power
summary_table <- multi_trait_df %>%
  dplyr::mutate(
    r_val = round(Correlation, 2),
    p_val = signif(Pvalue, 2),
    pw_val = round(Power_all, 2)
  ) %>%
  dplyr::select(Gene, PANSS, r_val, p_val, pw_val) %>%
  tidyr::pivot_wider(
    names_from = PANSS,
    values_from = c(r_val, p_val, pw_val)
  )

#Prepare matrices for heatmap visualization

#Correlation matrix
heat_corr_mat2 <- multi_trait_df %>%
  dplyr::select(Gene, PANSS, Correlation) %>%
  tidyr::pivot_wider(names_from = PANSS, values_from = Correlation) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

#Text annotations for each heatmap cell
heat_text_mat2 <- multi_trait_df %>%
  dplyr::mutate(
    cell_label = paste0(
      "r=", round(Correlation, 2),
      "\np=", signif(Pvalue, 2),
      "\npw=", round(Power_all, 2)
    )
  ) %>%
  dplyr::select(Gene, PANSS, cell_label) %>%
  tidyr::pivot_wider(names_from = PANSS, values_from = cell_label) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

#P-value matrix
heat_p_mat2 <- multi_trait_df %>%
  dplyr::select(Gene, PANSS, Pvalue) %>%
  tidyr::pivot_wider(names_from = PANSS, values_from = Pvalue) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

#Power matrix
heat_power_mat2 <- multi_trait_df %>%
  dplyr::select(Gene, PANSS, Power_all) %>%
  tidyr::pivot_wider(names_from = PANSS, values_from = Power_all) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

#Define color scale for correlation values
col_fun_nat <- colorRamp2(
  c(-1, 0, 1),
  c("#3B4CC0", "white", "#B40426")
)

#Generate heatmap of gene-PANSS correlations
ht_final <- Heatmap(
  heat_corr_mat2,   # <-- NOT transposed
  
  name = "Pearson r",
  col = col_fun_nat,
  
  cluster_rows = TRUE,       # cluster genes
  cluster_columns = FALSE,   # keep PANSS order N P G T
  
  row_names_side = "right",
  row_title_side = "left",
  row_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  
  column_title = "PANSS Groups",
  column_names_side = "bottom",
  row_title = "Genes",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  rect_gp = gpar(col = "grey85", lwd = 0.6),
  
  heatmap_legend_param = list(
    title = "r",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 11),
    legend_height = unit(4, "cm"),
    border = "black"
  ),
  
  #Annotate only significant and high-power cells
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    if (!is.na(heat_p_mat2[i, j]) &&
        heat_p_mat2[i, j] <= 0.05 &&
        heat_power_mat2[i, j] >= 0.8) {
      
      grid.text(
        paste0(
          sprintf("r=%.2f", heat_corr_mat2[i, j]), "\n",
          sprintf("p=%.3f", heat_p_mat2[i, j])
        ),
        x, y,
        gp = gpar(fontsize = 9, fontface = "bold")
      )
    }
    
  }
)

#Draw the final heatmap
draw(ht_final)
