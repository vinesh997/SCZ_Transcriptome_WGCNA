#loading the libraries
library(BiocManager)
library(WGCNA)
library(tidyr)
library(DESeq2)
library(dplyr)
library(tibble)
library(limma)
library(ggplot2)
library(ggfortify)

#setting the working directory
setwd("") #set the working directory where the files are stored.

#Reading the raw count matrix (raw_counts) 
Data<- read.delim("raw_counts.txt", check.names = F, sep = " ") #dim 59050;29

#filtering the genes 
smallestGroupSize <- 13
keep <- rowSums(Data >= 10) >= smallestGroupSize
data_filt<- Data[keep,] #dim 18840;29

#retaining only SCZ samples
data_filt<- t(data_filt)
data_filt<- as.data.frame(data_filt)
data_filt<- rownames_to_column(data_filt, var = "Sample")
Pheno<- read.csv("filename.csv") #dim 16;14, reading the phenodata (Phenodata 4) 
dat_sub<- data_filt[data_filt$Sample %in% Pheno$Sample,]
dat_sub<- dat_sub %>% remove_rownames() %>% column_to_rownames(var = "Sample")
dat_sub<-t(dat_sub) #dim 18840;16

#creating a DESeq DataSet object
dds <- DESeqDataSetFromMatrix(dat_sub,
                              Pheno,
                              design = ~ 1)

#countmatrix from DESeqDataSet
count_matrix <- counts(dds) #dim 18840;16

# Transforming the data using vst
vsd<- varianceStabilizingTransformation(dds)
mat<- assay(vsd)


#creating a model with all the covariates to be adjusted
design.Covariates <- model.matrix(~ Age + Sex + Duration + NEUTROPHILS + LYMPHOCYTES + 
                                    MONOCYTES + EOSINOPHILS + BASOPHILS, data=Pheno)

# Adjusting covariates using limma 
mat.adjusted <- limma::removeBatchEffect(mat, batch = Pheno$Batch, covariates = design.Covariates)

#visualizing for batch correction on PCA
autoplot(prcomp(t(mat), scale. = T), data = Pheno, colour="Batch") + ggtitle("Before Batch Correction")
autoplot(prcomp(t(mat.adjusted), scale. = T), data = Pheno, colour= "Batch") + ggtitle("After Batch Correction")

#transposing the data for WGCNA 
mat_adjusted_t<- t(mat.adjusted)



########################## WGCNA analysis ####################################


#Checking of missing values.
gsg <- goodSamplesGenes(mat_adjusted_t, verbose = 3);
gsg$allOK

#cluster the samples
sampleTree <- hclust(dist(mat_adjusted_t), method = "average");

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

#Reading and editing the trait data (Phenodata 5) to be used for correlation 
PS<- read.csv(" .csv", check.names = F)
PS<- column_to_rownames(PS, var = "Scale")
traitData <-as.data.frame(t(PS)) #dim 16;22
symptom_names <- rownames(PS)
PS<- t(mutate_all(PS, as.numeric))
colnames(PS) <- symptom_names
PS<- as.data.frame(PS) #dim 16;22

#checking the sequence of samples 
Samples <- rownames(mat_adjusted_t)
traitData<- rownames_to_column(traitData, var = "Samples")
traitRows<- match(Samples, traitData$Samples)
datTraits<- traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]
collectGarbage()


####### Re-cluster samples with clinical traits
sampleTree2<- hclust(dist(mat_adjusted_t), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors<- numbers2colors(datTraits, signed = TRUE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = colnames(datTraits),
                    main = "Sample dendrogram and trait heatmap")

# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 30, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(mat_adjusted_t,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

# Scale-free topology fit index as a function of the soft-thresholding power
sizeGrWindow(12,9)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=power,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off 
abline(h = 0.8, col="red") 

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=power, cex=cex1,col="red")



############# MODULE CONSTRUCTION & TRAIT RELATIONSHIPS ########################

# Constructing signed network in single block using sft power
cor <- WGCNA::cor
#selecting power as 18 according to WGCNA FAQ
net<- blockwiseModules(mat_adjusted_t, power = 18, maxBlockSize = 20000,  
                       deepSplit = 0, networkType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = FALSE,
                       saveTOMs = FALSE, verbose = 5)

# get number of genes for each module and plot dendrogram
table(net$colors)
sizeGrWindow(12, 9)
mergedColors<- labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], net$colors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(net$dendrograms[[1]], cbind(net$unmergedColors, net$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)



## Relating modules with external traits
moduleColors<- net$colors
MEs<- net$MEs
# Define number of genes and samples
nGenes<- ncol(mat_adjusted_t)
nSamples<- nrow(mat_adjusted_t)


# Recalculate MEs with color labels
MEs0<- moduleEigengenes(mat_adjusted_t, moduleColors)$eigengenes
MEs<- orderMEs(MEs0)
moduleTraitCor<- cor(MEs, PS, use = "p")
moduleTraitPvalue<- corPvalueStudent(moduleTraitCor, nSamples)


# Will display correlations and their p-values
textMatrix<- paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

## Display the correlation values within a heatmap plot
sizeGrWindow(10,6)
par(mar = c(10, 9, 1, 1))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(PS),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#dev.off()


############################## INTRAMODULAR ANALYSIS ###########################

#Defining the variable 
Trait<- as.data.frame(PS$Activation)
names(Trait)<- "Trait"

# names (colors) of the modules
modNames<- substring(names(MEs), 3)

# calculating module membership and p-value by using Pearson correlation between expression data and module eigengens
geneModuleMembership<- as.data.frame(cor(mat_adjusted_t, MEs, use = "p"))
MMPvalue<- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))


names(geneModuleMembership)<- paste("MM", modNames, sep="")
names(MMPvalue)<- paste("p.MM", modNames, sep="")

# calculating gene significance and p-value
geneTraitSignificance<- as.data.frame(cor(mat_adjusted_t, Trait, use = "p")) #pearson method
GSPvalue<- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance)<- paste("GS.", names(Trait), sep="")
names(GSPvalue)<- paste("p.GS.", names(Trait), sep="")

module<- "royalblue"   ######################### putting the color below the plot
column<- match(module, modNames)
moduleGenes<- moduleColors == module
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Guilt Feelings",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, 
                   col = module, abline = TRUE)

# Add vertical and horizontal cutoff lines
abline(v = 0.8, col = "red", lty = 2)  # Vertical line at X = 0.8
abline(h = 0.2, col = "blue", lty = 2) # Horizontal line at Y = 0.2 

#Identifying most important genes for one determined characteristic inside of the cluster
geneInfo0<- data.frame(EST = colnames(mat_adjusted_t),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
modOrder<- order(-abs(cor(MEs, Trait, use = "p")))
modOrder<- modOrder[modOrder > 0 & modOrder <= ncol(geneModuleMembership)]
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder<- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Trait))
geneInfo<- geneInfo0[geneOrder,]
#write.csv(geneInfo, file = " ")


# Subset genes belonging to a module
ModuleGenes<- geneInfo[geneInfo$moduleColor == "cyan", ]
#write.csv(ModuleGenes, " ")

# Filter the GS-MM file
filtered_data<- ModuleGenes[ModuleGenes$GS.Trait > 0.2 & ModuleGenes$MM.cyan >  0.8, ]

# View the filtered data
print(filtered_data)

# Save the filtered data to a new file
#write.csv(filtered_data, " ", row.names = FALSE)
