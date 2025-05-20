**Title:** Cellular and Gene Expression Changes in Peripheral Blood as Markers for Schizophrenia Symptoms: Exploratory Findings from a Pilot Study in the Indian Population

**Data type:** Hemogram data and RNA-seq data of SCZ and healthy control (HC) participants

**About the script:** The analysis is divided into two sections. 1- Hemogram analysis for partial correlation of hemogram data with the PANSS scores (symptom severity) and 2- Transcriptome analysis for differential gene expression analysis (DGEA) and Weighted gene coexpression analysis (WGCNA) for association with the symptom severity. The phenodata information for partial correlation and WGCNA can be obtained on request.  

**Hemogram analysis** 

1-Partial_Correlation.R- The correlation analysis between hemogram and the PANSS scores was performed using this script. For the initial analysis, a single file (Phenodata 2) is generated using data from only SCZ participants. The columns for Phenodata 2 are as follows: Sample_ID, Positive, Negative, General, Anergia, Thought Disturbance, Activation, Paranoid, Depression, Composite, Total, Neutrophils, Basophils, Eosinophils, Lymphocytes, Monocytes, Hemoglobin, RBC count, Platelet, NLR, PLR, MLR, Tobacco, Alcohol, Cannabis, Sex, Age, and Duration of Disorder. The next file (Phenodata 3) for filtration is generated using data from both SCZ and HC samples. The columns for Phenodata 3 are as follows: Sample_ID, Groups, Neutrophils, Basophils, Eosinophils, Lymphocytes, Monocytes, Hemoglobin, RBC count, Platelet, NLR, PLR, MLR.

**Transcriptome analysis** 
For differential gene expression analysis and WGCNA, raw RNA-seq files were processed using the scripts (https://github.com/satyajeetkhare/RNASeq_HISAT2_Gencode_hg38). The fastq files were obtained from INDA-CA (INCARP000275 and  INCARP000309). The raw count matrix generated was used for DGEA and WGCNA.

2A- Differential_Gene_Expression.R: The DEGs used for filtering the module genes were obtained through this script. The Phenodata 1 consisted of columns: Sample_ID, Age, Sex, Batch, Diagnosis, Neutrophils, Basophils, Eosinophils, Lymphocytes, and Monocytes. The raw count matrix (Raw_count1) was pre-processed to obtain filtered Raw_count2, which was used for DGEA as well as WGCNA.

2B- WGCNA_Transcriptome.R: The correlation of Transcriptome profile with the external trait (PANSS scores) is performed using this script. The filtered raw count matrix (raw_count2) was used in this step. The Phenodata 4 for this step consists of the following columns: Sample ID, Batch, Sex, Age, Diagnosis, Duration, Neutrophils, Basophils, Eosinophils, Lymphocytes, Monocytes, Tobacco, Alcohol, and Cannabis. The adjustment of the raw count matrix for batch, age, sex, duration of disorder, and hemogram was performed. The external trait data (Phenodata 5) consisted of Sample_ID, Positive, Negative, General, Anergia, Thought Disturbance, Activation, Paranoid, Depression, Composite, Total, Neutrophils, Basophils, Eosinophils, Lymphocytes, Monocytes, Tobacco, Alcohol, Cannabis, Sex, Age, Batch, and Duration of Disorder. 

Contact: vineshkamble23@gmail.com
