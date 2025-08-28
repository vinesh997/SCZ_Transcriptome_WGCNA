**Title:** Gene Expression Changes in Peripheral Blood as Markers for Anergia in Schizophrenia: Exploratory Findings in the Indian Population

**Data type:** RNA-seq data of SCZ and healthy control (HC) participants

**About the script:** 
The analysis was performed in two sections: 
1- Transcriptome analysis for differential gene expression analysis (DGEA) 
2- Weighted gene coexpression analysis (WGCNA) for association with the symptom severity. 
The phenodata information can be obtained on request.  

For differential gene expression analysis and WGCNA, raw RNA-seq files were processed using the scripts(https://github.com/satyajeetkhare/RNASeq_HISAT2_Gencode_hg38). The fastq files were obtained from INDA-CA (INCARP000309). The raw count matrix (raw_counts) generated was used for DGEA and WGCNA.

1- Differential_Gene_Expression.R: The DEGs used for filtering the module genes were obtained through this script. The Phenodata 1 consisted of columns: Sample_ID, Age, Sex, Batch, Diagnosis, Neutrophils, Basophils, Eosinophils, Lymphocytes, and Monocytes.

2- WGCNA_Transcriptome.R: The correlation of Transcriptome profile with the external trait (PANSS scores) is performed using this script. The raw count matrix (raw_counts) was used in this step and filtered further to retain only SCZ samples. The Phenodata 4 for this step consists of the following columns: Sample ID, Batch, Sex, Age, Diagnosis, Duration, Neutrophils, Basophils, Eosinophils, Lymphocytes, Monocytes, Tobacco, Alcohol, and Cannabis. The adjustment of the raw count matrix for batch, age, sex, duration of disorder, and hemogram was performed. The external trait data (Phenodata 5) consisted of Sample_ID, Positive, Negative, General, Anergia, Thought Disturbance, Activation, Paranoid, Depression, Composite, Total, Neutrophils, Basophils, Eosinophils, Lymphocytes, Monocytes, Tobacco, Alcohol, Cannabis, Sex, Age, Batch, and Duration of Disorder. 

Contact: vineshkamble23@gmail.com
