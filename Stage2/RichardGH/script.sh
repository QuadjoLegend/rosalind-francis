#!/bin/bash

#Creating needed directories...
mkdir datasets ref vcf mapping qc trim trimmedqc repaired

#Downloading the datasets...
wget -P ./datasets https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz
wget -P ./datasets https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz
wget -P ./datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz
wget -P ./datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz
wget -P ./datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz
wget -P ./datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz
wget -P ./datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz
wget -P ./datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz
wget -P ./datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz
wget -P ./datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz
wget -P ./datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz
wget -P ./datasets https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz

#Downloading the reference files...
wget -P ./ref https://zenodo.org/records/10886725/files/Reference.fasta
wget -P ./ref https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta

#Performing Quality Controls...
fastqc datasets/*.fastq.gz --outdir qc
echo "QC complete, output saved in qc directory."

SAMPLES=("ACBarrie_R" "Alsen_R" "Baxter_R" "Chara_R" "Drysdale_R" "ERR8774458_")

# Repair raw reads
for smp in "${SAMPLES[@]}"; do
    repair.sh \
      in1="datasets/${smp}1.fastq.gz" \
      in2="datasets/${smp}2.fastq.gz" \
      out1="repaired/${smp}1_rep.fastq.gz" \
      out2="repaired/${smp}2_rep.fastq.gz"
done

#Trim repaired reads
for smp in "${SAMPLES[@]}"; do
  fastp \
    -i "repaired/${smp}1_rep.fastq.gz" \
    -I "repaired/${smp}2_rep.fastq.gz" \
    -o "trim/${smp}1_trim_rep.fastq.gz" \
    -O "trim/${smp}2_trim_rep.fastq.gz" \
    --html "qc/${smp}_fastq_rep.html" \
    --json "qc/${smp}_fastq_rep.json"
done

#Performing qc on trimmed reads...
fastqc trim/*rep.fastq.gz --outdir trimmedqc

#Indexing reference files...
bwa index ref/reference.fasta
samtools faidx ref/reference.fasta
bwa index ref/Reference.fasta
samtools faidx ref/Reference.fasta

echo "Performing allignments on trimmed reads..."
bwa mem \
  "ref/Reference.fasta" \
  "trim/ERR8774458_1_trim_rep.fastq.gz" \
  "trim/ERR8774458_2_trim_rep.fastq.gz" \
  > "mapping/ERR8774458_.sam"

for smp in "${SAMPLES[@]}"; do
  if [ "$smp" != "ERR8774458_" ]; then
    bwa mem \
      "ref/reference.fasta" \
      "trim/${smp}1_trim_rep.fastq.gz" \
      "trim/${smp}2_trim_rep.fastq.gz" \
      > "mapping/${smp}.sam"
  fi
done

#Convert sam to bam
for smp in "${SAMPLES[@]}"; do
    samtools view -b -S -o "mapping/${smp}.bam" "mapping/${smp}.sam"
done

#Delete sam files


#Sort bam files
for smp in "${SAMPLES[@]}"; do
    samtools sort \
        "mapping/${smp}.bam" \
        -o "mapping/${smp}.sorted.bam"
done
#Index sorted bam
samtools index -M mapping/*.sorted.bam


# Variant calling for ERR8774458_
bcftools mpileup -Ou -f "ref/Reference.fasta" "mapping/ERR8774458_.sorted.bam" | \
  bcftools call -Ov -mv > "vcf/ERR8774458_.vcf"

# Variant calling for other samples
for smp in "${SAMPLES[@]}"; do
  if [ "$smp" != "ERR8774458_" ]; then
    bcftools mpileup -Ou -f "ref/reference.fasta" "mapping/${smp}.sorted.bam" | \
      bcftools call -Ov -mv > "vcf/${smp}.vcf"
  fi
done


echo "Pipeline created by Richard Agyekum"

