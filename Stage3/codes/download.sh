#!/bin/bash
mkdir data
# Set output directory
output_dir="data"

# Text file containing accession numbers (one per line)
accessions_file="samples.txt"

# Loop through each line in the text file
while IFS= read -r accession; do
  # Check if line is empty or a comment (starting with #)
  if [[ -z "$accession" || "$accession" =~ ^# ]]; then
    continue
  fi

  # Check if output files (R1 and R2) already exist (adjust based on naming)
  if [[ -f "$output_dir/$accession_1.fastq" && -f "$output_dir/$accession_2.fastq" ]]; then
    echo "Skipping download: $accession (files already exist)"
    continue
  fi

  # Download FASTQ files with SRA Toolkit
  fastq-dump --split-files -O "$output_dir" "$accession"


done < "$accessions_file"

gzip data/*
echo "Download, conversion, and separation completed!"

