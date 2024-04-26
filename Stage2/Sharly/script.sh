1.mkdir stage2
2.cd stage2
# To download datasets
mkdir datasets
cd datasets
3.wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz
  wget https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz
 https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz
  wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz
  wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz
  wget https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz
#To download reference genome
4.mkdir ref
5.cd ref
6. wget https://zenodo.org/records/10886725/files/Reference.fasta https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta
7.cd ..
#To perform quality control of datasets
8.sudo apt-get install fastqc
9.fastqc *_R1.fastq.gz
  fastqc *_R2.fastq.gz
 fastqc *.fastq.gz
#To trim datasets
10.sudo apt-get install fastp
11.fastp -i *_1.fastq.gz -I *_2.fastq.gz -o trimmed_R1.fastq -O trimmed_R2.fastq --qualified_quality_phred 20 --html report.html --json report.json
#To align trimmed datasets to reference genome
12.sudo apt-get install bwa
13. bwa index reference.fasta  bwa index Reference.fasta
14.bwa mem ref/Reference.fasta trimmed_R1.fastq trimmed_R2.fastq > ERR8774458.sam
15.bwa mem ref/reference.fasta trimmed_Alsen_R1.fastq trimmed_Alsen_R2.fastq > Alsen.sam
  bwa mem ref/reference.fasta trimmed_Baxter_R1.fastq trimmed_Baxter_R2.fastq > Baxter.sam
  bwa mem ref/reference.fasta trimmed_Chara_R1.fastq trimmed_Chara_R2.fastq > Chara.sam
  bwa mem ref/reference.fasta trimmed_Drysdale_R1.fastq trimmed_Drysdale_R2.fastq > Drysdale.sam
 bwa mem ref/reference.fasta trimmed_R1.fastq trimmed_R2.fastq > ACBarrie.sam
#To convert sam files to bam files
16.sudo apt-get install samtools
17.samtools view -hbo ERR8774458.bam ERR8774458.sam
samtools view -hbo Alsen.bam Alsen.sam
samtools view -hbo Baxter.bam Baxter.sam
samtools view -hbo Chara.bam Chara.sam
samtools view -hbo Drysdale.bam Drysdale.sam
samtools view -hbo ACBarrie.bam ACBarrie.sam
#To perform variant calling
18.sudo apt-get install bcftools
19.bcftools mpileup -Ou -f ref/Reference.fasta ERR8774458.bam | bcftools call -Ov -mv > ERR8774458.vcf
bcftools mpileup -Ou -f ref/reference.fasta Alsen.bam | bcftools call -Ov -mv > Alsen.vcf
bcftools mpileup -Ou -f ref/reference.fasta Baxter.bam | bcftools call -Ov -mv > Baxter.vcf
bcftools mpileup -Ou -f ref/reference.fasta Chara.bam | bcftools call -Ov -mv > Chara.vcf
bcftools mpileup -Ou -f ref/reference.fasta Drysdale.bam | bcftools call -Ov -mv > Drysdale.vcf
bcftools mpileup -Ou -f ref/reference.fasta ACBarrie.bam | bcftools call -Ov -mv > ACBarrie.vcf
