**Title:** Blood Cell Count and Transcriptome Profile Associates with Schizophrenia (SCZ) Symptom Severity: Insights from a Pilot Study of Indian Population

**Aim:** Association of peripheral blood parameters with the symptom severity in SCZ in the Indian population

**Data type:** RNA-seq data of SCZ and healthy control (HC) participants

**About the script:** The raw RNA-seq files were processed using the scripts (link). The fastq files are obtained from INDA-CA (INCARP000275 and  ). The raw count matrix is used for differential gene expression analysis and Weighted gene coexpression analysis (WGCNA). The scripts are as provided above, and the details of each script are mentioned below. 

1- Differential gene expression: The DEGs used for filtering the hub genes were obtained through this script. The pre-processing of the raw count matrix is performed in this stage, which is used for DEGs as well as WGCNA. 

2- Partial Correlation: The correlation between hemogram and the PANSS scores was performed. The filtering is based on the p-value as well as the mean hemogram count

3- WGCNA: The adjustment of the raw count matrix for batch, age, sex, duration of disorder, and hemogram is performed. Further, the correlation of the transcriptome (raw counts) with the symptom severity is also performed.
