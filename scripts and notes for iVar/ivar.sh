#!/bin/bash
output="/home/lbashor/cat_SARS2_manuscript_June2021/iVar/iVar_11_12_21/output/"

mkdir -p ${output}

cat InputFileList.txt | while read -r line
do

echo ${line}

bwa mem -t 32 MN985325.fasta ${line}_R1.fastq.gz ${line}_R2.fastq.gz | samtools view -b -F 4 -F 2048 | samtools sort -o ${output}${line}.sorted.bam   

ivar trim -b MN985325.primer.bed -p ${output}${line}.trimmed -i ${output}${line}.sorted.bam

samtools sort -o ${output}${line}.trimmed.sorted.bam ${output}${line}.trimmed.bam
samtools index ${output}${line}.trimmed.sorted.bam

samtools mpileup -A -d 0 -B ${output}${line}.trimmed.sorted.bam | ivar variants -p ${output}${line}_ivar -q 30 -t 0.03 -r MN985325.fasta -g MN985325.gff3

done
