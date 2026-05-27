**Title:** Cellular and Molecular Markers for Negative symptoms in Schizophrenia: Exploratory Findings in the Indian Population

**Data type:** RNA-seq data of SCZ and healthy control (HC) participants

**About the script:** 
The analysis was performed in three sections: 
1- Partial correlation of cellular proportions obtained from deconvolution of RNA-seq data with Negative syndrome severity. 
2- Weighted Gene Coexpression Network Analysis (WGCNA) for identification of gene clusters associated with Negative syndrome. 
3- Differential gene expression analysis (DGEA) for studying the patterns in SCZ and HC. 
The phenodata information can be obtained on request.  

**Methods**
For the analysis, raw RNA-seq files were processed using the scripts (https://github.com/satyajeetkhare/RNASeq_HISAT2_Gencode_hg38). The fastq files were obtained from INDA-CA (INCARP000309). The raw count matrix (raw_counts) generated was used for further analysis.

1- Partial_Correlation: The cellular proportions obtained using CIBERSORTx were transformed to obtain absolute cell counts. These cell proportions (18 subtypes- PBMC components) were correlated with external trait (PANSS scores) using this script.

2- Weighted gene coexpression network analysis- The raw count matrix (raw_counts) was used in this step and filtered further to retain only SCZ samples. The Phenodata for this step consists of the following columns: Sample ID, Batch, Sex, Age, Diagnosis, DoD, P, N, G, T, and T cells CD4 memory activated. The transcriptome profile was adjusted for Age, Sex, Batch, and T cells CD4 memory activated. The adjusted matrix was subjected to WGCNA and the modules were correlated with the symptom severity (P, N, G, and T) in this script. 

2- Differential_Gene_Expression.R: The expression pattern was validated using this script by first performing differential gene expression. The Phenodata consisted of columns: Sample_ID, Age, Sex, Batch, Diagnosis, and T cells CD4 memory activated.


Contact: vineshkamble23@gmail.com
