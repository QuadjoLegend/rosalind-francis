#Activating the conda base environment
conda init
#Activating the funtools environment
conda activate funtools
#To view the packages installed
conda list


#Create a Data folder to keep the Dataset    
mkdir Data

#To view the folders,files in the given folder  
ls

#To get inside the  'Data' folder
cd Data 

#Download the Fasta files of  the given  Genome Sequences and the reference sequence from the  link
wget https://zenodo.org/records/10426436/files/ERR8774458_1.fastq.gz
wget https://zenodo.org/records/10426436/files/ERR8774458_2.fastq.gz
wget https://zenodo.org/records/10886725/files/Reference.fasta

#To check the downloaded files
ls
#To  do a Quality check  by using fastqc tool
fastqc ERR8774458_1.fastq.gz  ERR8774458_2.fastq.gz

#To create a folder  named 'Ref' for the reference data file.
mkdir Ref
#Moving the Reference file from current directory to the 'Ref' folder.
mv Reference.fasta /root/Ref

#Multiqc is a summary report of the genome sequences
multiqc ERR8774458_1_fastqc.zip   ERR8774458_2_fastqc.zip

# Fastp is a tool used for read quality control , filtering and adapter trimmming.
fastp -i ERR8774458_1.fastq -I ERR8774458_2.fastq -o ERR8774458_trimmed_1.fastq -O ERR8774458_trimmed_2.fastq

#It is used to create  an index of a reference genome.
#This index helps in the faster performance of the sequence alignment against the reference genome during the mapping process.
bwa index /root/Ref/Reference.fasta
#This performs fast & accurate alignment of sequencing reads to  a reference genome.
bwa mem Reference.fasta ERR8774458_trimmed_1.fastq ERR8774458_trimmed_2.fastq > Alignment.sam

#Samtools is used to convert sam file into bam fileformat.
samtools view -bS Alignment.sam > Alignment.bam

#Create a folder Bam to keep the bam files
mkdir Bam

#This'flagstat'command gives the summary statistics about alignments in a Bam file.
samtools  flagstat   Alignment.bam
#This'sort'command sorts alignments in a BAM file by their coordinates,so improve the efficiency for downstream analysis. 
samtools sort Alignment.bam -o Alignment_sorted.bam
#It is used to create  an index of a reference genome.
samtools index Alignment_sorted.bam
#This'flagstat'command gives the summary statistics about alignments in a sorted Bam file.
samtools  flagstat   Alignment.sorted.bam

#Create a folder'VCF'for the vcf files.
mkdir VCF
#Using freebayes to analyze the the alignments in the BAM file against the reference genome and detect the variants in the output vcf fileformat.
freebayes  -f   /root/Ref/Reference.fasta  -b   /root/Bams/Alignment_sorted.bam  --vcf /root/VCF/Alignment.vcf
#To compress the vcf file 
bgzip   /root/VCF/Alignment.vcf
# Bcftools is used to analyze the variants , SNPs , indels.
#to create an index forthe vcf file. 
bcftools index /root/VCF/Alignment.vcf.gz
#To view the SNPs(single nucleotide polymorphisms) and indels(insertions and deletions) in the vcf file.
#This will not include the lines starting from'#'
bcftools view -v snps   /root/VCF/Alignment.vcf.gz|grep  -v  -c  '^#'
bcftools view -v indels  /root/VCF/Alignment.vcf.gz|grep -v   -c   '^#'

