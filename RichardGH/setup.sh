#!/bin/bash

conda init
conda activate
conda create -n funtools
conda activate funtools
conda install bwa -c bioconda
conda install samtools -c  bioconda
conda install bcftools -c  bioconda
conda install fastp -c  bioconda
conda install fastqc -c  bioconda
conda install bbtools -c AgBiome
