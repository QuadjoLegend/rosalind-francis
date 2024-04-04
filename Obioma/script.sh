#Creating needed directories
mkdir datasets ref vcf mapping qc trim trimmedqc

#Downloading the datasets
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

#Downloading the reference files
wget -P ./ref https://zenodo.org/records/10886725/files/Reference.fasta
wget -P ./ref https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta

#Performing Quality Control
sudo apt-get install fastqc
fastqc ./datasets/*.fastq.gz --outdir qc

#Performing Trimming
SAMPLES=("ACBarrie_R" "Alsen_R" "Baxter_R" "Chara_R" "Drysdale_R" "ERR8774458_")
sudo apt-get install fastp
for smp in "${SAMPLES[@]}"; do
  fastp \
    -i "datasets/${smp}1.fastq.gz" \
    -I "datasets/${smp}2.fastq.gz" \
    -o "trim/${smp}1_trim.fastp.gz" \
    -O "trim/${smp}2_trim.fastp.gz" \
    --html "qc/${smp}_fastp.html" \
    --json "qc/${smp}_fastp.json"
done

#Performing qc on trimmed reads
fastqc trim/*.fastp.gz --outdir trimmedqc

#Indexing reference files
sudo apt-get install bwa
sudo apt-get install samtools
bwa index ref/reference.fasta
samtools faidx ref/reference.fasta
bwa index ref/Reference.fasta
samtools faidx ref/Reference.fasta

#Performing allignments on trimmed reads
bwa mem \
  "ref/Reference.fasta" \
  "trim/ERR8774458_1_trim.fastp.gz" \
  "trim/ERR8774458_2_trim.fastp.gz" \
  > "mapping/ERR8774458_.sam"

for smp in "${SAMPLES[@]}"; do
  if [ "$smp" != "ERR8774458_" ]; then
    bwa mem \
      "ref/reference.fasta" \
      "trim/${smp}1_trim.fastp.gz" \
      "trim/${smp}2_trim.fastp.gz" \
      > "mapping/${smp}.sam"
    fi
    done
    
#Convert sam to bam
for smp in "${SAMPLES[@]}"; do
    samtools view -b -S -o "mapping/${smp}.bam" "mapping/${smp}.sam"
done

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




    
