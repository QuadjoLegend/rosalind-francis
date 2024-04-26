#!/bin/bash
set -e

# Function to download files with error handling
download_file() {
    url="$1"
    directory="$2"
    filename=$(basename "$url")
    
    if [ ! -f "$directory/$filename" ]; then
        wget -P "$directory" "$url"
    else
        echo "$filename already exists in $directory. Skipping download."
    fi
}

# Create directories for the analysis output
mkdir -p samples fastq_out trimmed_samples fastqc_trimmed_reads variant_calls alignment_results bam_files sorted_bam_files multiqc_report multiqc_report_trimmed


# Download fastq files to the "samples" directory
download_file "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz" "samples"
download_file "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz" "samples"
download_file "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz" "samples"
download_file "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz" "samples"
download_file "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz" "samples"
download_file "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz" "samples"
download_file "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz" "samples"
download_file "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz" "samples"
download_file "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz" "samples"
download_file "https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz" "samples"

# Download Reference Genome
download_file "https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta" "samples"

# Unzip fastq files
gzip -d samples/*.fastq.gz || true

# Run FastQC for all .fastq files in the "samples" directory
fastqc samples/*.fastq -o fastq_out

# Run MultiQC
multiqc fastq_out -o multiqc_report

# Run fastp for each sample
for SAMPLE in "ACBarrie" "Alsen" "Baxter" "Chara" "Drysdale"; do
  fastp \
    -i "samples/${SAMPLE}_R1.fastq" \
    -I "samples/${SAMPLE}_R2.fastq" \
    -o "trimmed_samples/${SAMPLE}_R1.fastq.gz" \
    -O "trimmed_samples/${SAMPLE}_R2.fastq.gz" \
    --html "trimmed_samples/${SAMPLE}_fastp.html" 
done

# Run FastQC for all .fastq files in the "trimmed_samples" directory
fastqc trimmed_samples/*.fastq.gz -o fastqc_trimmed_reads

# Run MultiQC for trimmed reads
multiqc fastqc_trimmed_reads -o multiqc_report_trimmed

# Index the reference FASTA file using BWA if index files are not present
if [ ! -f "samples/reference.fasta.bwt" ]; then
    bwa index samples/reference.fasta
fi

# Run BWA for each sample
for SAMPLE in "ACBarrie" "Alsen" "Baxter" "Chara" "Drysdale"; do
  bwa mem \
    samples/reference.fasta \
    "trimmed_samples/${SAMPLE}_R1.fastq.gz" \
    "trimmed_samples/${SAMPLE}_R2.fastq.gz" \
    > "alignment_results/${SAMPLE}_aligned_reads.sam" || true
done

# Convert SAM to BAM and sort for each sample
for SAMPLE in "ACBarrie" "Alsen" "Baxter" "Chara" "Drysdale"; do
  samtools view -bS "alignment_results/${SAMPLE}_aligned_reads.sam" | \
  samtools sort -o "sorted_bam_files/${SAMPLE}_aligned_reads_sorted.bam" -
  samtools index "sorted_bam_files/${SAMPLE}_aligned_reads_sorted.bam"
done

# Perform variant calling for each sample
for SAMPLE in "ACBarrie" "Alsen" "Baxter" "Chara" "Drysdale"; do
  bcftools mpileup -Ou -f samples/reference.fasta "sorted_bam_files/${SAMPLE}_aligned_reads_sorted.bam" | \
  bcftools call -Ov -mv > "variant_calls/${SAMPLE}_variants.vcf"
done

# Display end message
echo "Pipeline completed successfully."

