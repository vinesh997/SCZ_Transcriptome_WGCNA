#!/bin/bash
folder="Directory of the raw fastq files"
hisat2="Directory for hisat2-2.1.0"

mkdir "${folder}/hisat2"

# Human genome index
genome="Directory of human reference genome "

for i in $(ls ${folder}/cutadapt2/*.gz | xargs -n 1 basename | sed 's/\(.*\)_.*/\1/' | sed 's/\(.*\)_.*/\1/' | sort -u)

do
echo ${i}

echo -e "\n""Sample ${i} alignment">>${folder}/hisat2/hisat2.log

# for human
(${hisat2}/hisat2 -p 20 -5 15 -3 5 --dta --known-splicesite-infile ${genome}/gencode.v31.annotation_splicesites.txt -x ${genome}/GRCh38.p12.genome -1 ${folder}/cutadapt/${i}_R1.fastq.gz,${folder}/cutadapt2/${i}_R1.fastq.gz -2 ${folder}/cutadapt/${i}_R2.fastq.gz,${folder}/cutadapt2/${i}_R2.fastq.gz | samtools view -@ 20 -bS - | samtools sort -@ 20 - -o ${folder}/hisat2/${i}_sorted.bam) 2>> ${folder}/hisat2/hisat2.log

# Index sorted bam file for browser
samtools index -@ 20 ${folder}/hisat2/${i}_sorted.bam

# For count matrix generation
echo "${folder}/hisat2/${i}_sorted.bam">>${folder}/nutshell/bam_files.txt

done
