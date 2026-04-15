**Title:** Gene Expression Changes in Peripheral Blood as Markers for Negative symptoms in Schizophrenia: Exploratory Findings in the Indian Population

**Data type:** RNA-seq data of SCZ and healthy control (HC) participants

**About the script:** 
The analysis was performed in two sections: 
1- Pearson's Correlation analysis for association of Transcriptome profile with the symptom severity
2- Transcriptome analysis for differential gene expression analysis (DGEA). 
The phenodata information can be obtained on request.  

**Methods**
For differential gene expression analysis and correlation analysis, raw RNA-seq files were processed using the scripts (https://github.com/satyajeetkhare/RNASeq_HISAT2_Gencode_hg38). The fastq files were obtained from INDA-CA (INCARP000309). The raw count matrix (raw_counts) generated was used for correlation analysis and DGEA.

1- Pearson_Correlation.R: The correlation of Transcriptome profile with the external trait (PANSS scores) is performed using this script. The raw count matrix (raw_counts) was used in this step and filtered further to retain only SCZ samples. The Phenodata for this step consists of the following columns: Sample ID, Batch, Sex, Age, Diagnosis, Duration, Lymphocytes, Monocytes, P, N, G, and T. The transcriptome profile was correlated with the symptom severity (P, N, G, and T) along with covariates and confounders (Batch, Sex, Age, Duration, Lymphocytes, Monocytes).

2- Differential_Gene_Expression.R: The expression pattern was validated using this script by first performing differential gene expression. The Phenodata consisted of columns: Sample_ID, Age, Sex, Batch, Diagnosis, Neutrophils, Basophils, Eosinophils, Lymphocytes, and Monocytes.


Contact: vineshkamble23@gmail.com
