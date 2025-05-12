# set working directory
setwd("")

# sample sheet should have "bam file" location, "sample name" and "group" columns
bam_files <- read.delim("./bam_files.txt", header = TRUE, sep = "\t")
bam_files$bam_files <- as.character(bam_files$bam_files)


library(Rsubread)
# Run featureCounts for transcripts 
Counts_fc <- featureCounts(bam_files$bam_files, 
                           annot.ext = "/data/genomes/gencode/GRCh38_p12/gencode.v31.annotation.gtf", 
                           isGTFAnnotationFile = TRUE, GTF.attrType = "gene_name", isPairedEnd=TRUE, nthreads = 30)

library(DESeq2)
# make simple sample names
colnames(Counts_fc$counts) <- bam_files$Sample
# Create ColData object
samples<-subset(bam_files, select=c("Batch","Age","Sex","Diagnosis"))
# assign sample name as rowname
rownames(samples) <- bam_files[,2]
# Create mycols object
mycols = data.frame(row.names = factor(bam_files$Sample))
# Check the column names of count matrix and row names of "mycols" are same
all(mycols %in% colnames(Counts_fc$counts))

# Create DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = Counts_fc$counts, colData = samples, 
                              design = ~ Batch + Sex + Diagnosis)

# Create DESeq
dds2 <- DESeq(dds)

#raw counts all
raw_counts <- counts(dds2, normalized=FALSE)
write.table(raw_counts, file = "raw_counts_2.txt", quote = FALSE, sep = "\t", 
            row.names = TRUE)
