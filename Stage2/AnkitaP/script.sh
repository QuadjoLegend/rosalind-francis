#!/bin/bash

#To create the folders
mkdir Datasets ref qc trim trimmedqc alignment VCF 

#To download the datasets
wget -P ./Datasets https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz
wget -P ./Datasets https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz
wget -P ./Datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz
wget -P ./Datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz
wget -P ./Datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz
wget -P ./Datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz
wget -P ./Datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz
wget -P ./Datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz
wget -P ./Datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz
wget -P ./Datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz
wget -P ./Datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz
wget -P ./Datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz
 
#To download the reference files
wget -P ./ref https://zenodo.org/records/10886725/files/Reference.fasta
wget -P ./ref https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta


#To run the fastqc for quality control
fastqc ./datasets/*.fastq.gz --outdir qc
echo "QC completed,outputs were saved in the qc folder."


samples=("ACBarrie_R" "Alsen_R" "Baxter_R" "Chara_R" "Drysdale_R" "ERR8774458_")

for smp in "${samples[@]}"; do
  fastp \
    -i "Datasets/${smp}1.fastq.gz" \
    -I "Datasets/${smp}2.fastq.gz" \
    -o "trim/${smp}1_trim.fastp.gz" \
    -O "trim/${smp}2_trim.fastp.gz" \
    --html "qc/${smp}_fastp.html" \
    --json "qc/${smp}_fastp.json"
done
 
#To run the fastqc again on the trimmed reads to vie the quailty of the reads 
fastqc trim/*.fastp.gz --outdir trimmedqc

#To Index reference files
bwa index ref/reference.fasta
samtools faidx ref/reference.fasta
bwa index ref/Reference.fasta
samtools faidx ref/Reference.fasta

echo "Now performing the bwa(alignment) on trimmed reads"
bwa mem \
  "ref/Reference.fasta" \
  "trim/ERR8774458_1_trim.fastp.gz" \
  "trim/ERR8774458_2_trim.fastp.gz" \
  > "alignment/ERR8774458_.sam"

for smp in "${samples[@]}"; do
  if [ "$smp" != "ERR8774458_" ]; then
    bwa mem \
      "ref/reference.fasta" \
      "trim/${smp}1_trim.fastp.gz" \
      "trim/${smp}2_trim.fastp.gz" \
      > "alignment/${smp}.sam"
  fi
done

#To Convert sam file into bam fileformat
for smp in "${samples[@]}"; do
    samtools view -b -S -o "alignment/${smp}.bam" "alignment/${smp}.sam"
done


#To Sort the bam files
for smp in "${samples[@]}"; do
    samtools sort \
        "alignment/${smp}.bam" \
        -o "alignment/${smp}.sorted.bam"
done


#To Index sorted bam file
samtools index -M alignment/*.sorted.bam


#To run the Variant calling only for ERR8774458_
bcftools mpileup -Ou -f "ref/Reference.fasta" "alignment/ERR8774458_.sorted.bam" | \
  bcftools call -Ov -mv > "VCF/ERR8774458_.vcf"

#To run the Variant calling for other samples
for smp in "${samples[@]}"; do
  if [ "$smp" != "ERR8774458_" ]; then
    bcftools mpileup -Ou -f "ref/reference.fasta" "alignment/${smp}.sorted.bam" | \
      bcftools call -Ov -mv > "VCF/${smp}.vcf"
  fi
done

