#!/bin/bash
mkdir Mycobac
cd Mycobac
wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz https://zenodo.org/records/10886725/files/Reference.fasta
fastqc ERR8774458_1.fastq.gz  ERR8774458_2.fastq.gz
fastp -i ERR8774458_1.fastq.gz -I ERR8774458_2.fastq.gz -o trim_1.fastq.gz -O trim_2.fastq.gz
fastqc trim_1.fastq trim_2.fastq.gz
bwa index Reference.fasta
bwa mem Reference.fasta trim_1.fastq.gz trim_2.fastq.gz > Alignment.sam
samtools view -bS alignment.sam | samtools sort -o alignment_sorted.bam
samtools index alignment_sorted.bam
bcftools mpileup -Ou -f Reference.fasta alignment_sorted.bam | bcftools call -mv -Ov -o variants.vcf
cd ..


# Define datasets and their corresponding read URLs
SAMPLES=(
  "ACBarrie https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz"
  "Alsen https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz"
  "Baxter https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz"
  "Chara https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz"
  "Drysdale https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz"
)

# Create folders for each flower sample and the reference genome
mkdir -p reference_genome
for sample_data in "${SAMPLES[@]}"; do
  sample_name=$(echo "$sample_data" | cut -d' ' -f1)
  mkdir -p "$sample_name"
done

# Download read files and move them to respective folders
for sample_data in "${SAMPLES[@]}"; do
  sample_name=$(echo "$sample_data" | cut -d' ' -f1)
  read_urls=($(echo "$sample_data" | cut -d' ' -f2-))

  for url in "${read_urls[@]}"; do
    filename=$(basename "$url")
    wget -O "$sample_name/$filename" "$url"
  done
done

# Download reference genome
echo "Downloading reference genome..."
REFERENCE_URL="https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/reference.fasta"
wget -O "reference_genome/reference.fasta" "$REFERENCE_URL"

#QUALITYCONTROL-------------------------FASTQC

# Define flower sample folders
SAMPLES=(
  "ACBarrie"
  "Alsen"
  "Baxter"
  "Chara"
  "Drysdale"
)

# Perform FastQC analysis for each flower sample
for folder in "${SAMPLES[@]}"; do
  echo "Running FastQC for $folder..."
  fastqc "$folder"/*.fastq.gz -o "$folder"
done

#FASTP----------------------TRIM

# Perform Fastp analysis for each flower sample
for folder in "${SAMPLES[@]}"; do
  echo "Running Fastp for $folder..."
  fastp -i "$folder/${folder}_R1.fastq.gz" -I "$folder/${folder}_R2.fastq.gz" \
        -o "$folder/${folder}_R1.trimmed.fastq.gz" -O "$folder/${folder}_R2.trimmed.fastq.gz" \
        --html "$folder/${folder}_fastp_report.html"
done

#TRIMMED-----QC

# Perform FastQC analysis for trimmed reads of each flower sample
for folder in "${SAMPLES[@]}"; do
  echo "Running FastQC for trimmed reads of $folder..."
  fastqc "$folder/${folder}_R1.trimmed.fastq.gz" "$folder/${folder}_R2.trimmed.fastq.gz" -o "$folder"
done

#BWA-----INDEXING AND MAPPING

# Loop through each flower sample
for sample_name in "${SAMPLES[@]}"; do
  # Index the reference genome using BWA
  echo "Indexing reference genome for $sample_name..."
  bwa index "reference_genome/reference.fasta"

  # Perform BWA mapping
  echo "Performing BWA mapping for $sample_name..."
  bwa mem -t 4 "reference_genome/reference.fasta" "$sample_name/${sample_name}_R1.trimmed.fastq.gz" "$sample_name/${sample_name}_R2.trimmed.fastq.gz" > "${sample_name}/${sample_name}_alignment.sam"
done

#SAMTOBAM---------SORT

#for flower datasets
# Loop through each flower sample
for sample_name in "${SAMPLES[@]}"; do
  # Convert SAM to BAM and sort
  echo "Converting SAM to BAM and sorting for $sample_name..."
  samtools view -bS "${sample_name}/${sample_name}_alignment.sam" | samtools sort -o "${sample_name}/${sample_name}_alignment_sorted.bam"

  # Index the sorted BAM file
  echo "Indexing sorted BAM file for $sample_name..."
  samtools index "${sample_name}/${sample_name}_alignment_sorted.bam"
done

#VCF FILES-------------VARIANT CALLING

# Define path to the reference genome file
REFERENCE_PATH="reference_genome/reference.fasta"

# Loop through each flower sample
for sample_name in "${SAMPLES[@]}"; do
  # Perform variant calling
  echo "Performing variant calling for $sample_name..."
  bcftools mpileup -Ou -f "${REFERENCE_PATH}" "${sample_name}/${sample_name}_alignment_sorted.bam" | bcftools call -mv -Ov -o "${sample_name}/variants.vcf"
done