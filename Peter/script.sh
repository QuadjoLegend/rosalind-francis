mkdir HBD
mkdir HBD/raw_data

#Getting the files 
wget -P HBD/raw_data https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz \
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz \
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz \
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz \
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz \
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz \
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz \
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz \
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz \
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz 


#Quality Control
cd HBD
fastqc raw_data/*.fastq.gz
multiqc raw_data/*_fastqc.gz.zip

#trimming
mkdir trim
SAMPLES=("ACBarrie_R" "Alsen_R" "Baxter_R" "Chara_R" "Drysdale_R")

for smp in "${SAMPLES[@]}"; do
  fastp \
    -i "raw_data/${smp}1.fastq.gz" \
    -I "raw_data/${smp}2.fastq.gz" \
    -o "trim/${smp}1_trim.fastp.gz" \
    -O "trim/${smp}2_trim.fastp.gz" \
    --html "raw_data/${smp}_fastp.html" \
    --json "raw_data/${smp}_fastp.json"
done

#getting reference files
mkdir ref
wget -P ref https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta

##index the ref
cd HBD/ref
bwa index reference.fasta

#mapping
mkdir HBD/alignment
for smp in "${SAMPLES[@]}"; do
  bwa mem \
    ref/reference.fasta \
    "trim/${smp}1_trim.fastp.gz" \
    "trim/${smp}2_trim.fastp.gz" \
    > "alignment/${smp}_aligned.sam" || true
done

# Convert SAM to BAM and sort for each sample
for smp in "${SAMPLES[@]}"; do
      samtools view -b -S -o “alignment/${smp}_aligned.bam” “alignment/${smp}_aligned.sam”
done

#Sort and index bam files
for smp in "${SAMPLES[@]}"; do
    samtools sort \
        "alignment/${smp}_aligned.bam" \
        -o "alignment/${smp}_sorted.bam"
      samtools index alignment//${smp}_sorted.bam
done

##Variant Calling
Mkdir VCF
for smp in "${SAMPLES[@]}"; do
    bcftools mpileup -Ou -f ref/reference.fasta alignment//${smp}_sorted.bam | \
    bcftools call -Ov -mv > VCF/${smp}_variants.vcf
done
echo “Pipeline Done.”
