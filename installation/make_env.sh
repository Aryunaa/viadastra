#!/bin/bash
conda create --name myviadastra
conda activate myviadastra
conda install -c anaconda openjdk
conda install -c bioconda bedtools
conda install python=3.6
conda install -c anaconda numpy
conda install -c anaconda pandas scipy statsmodels
#conda install -c bioconda gatk4
conda install -c bioconda picard
conda install -c au-eoed gnu-parallel
conda install -c omnia subprocess32
conda install -c bioconda pysam
#conda install -c anaconda configparser
#conda install -c conda-forge r-sys
#conda install -c anaconda seaborn
#conda install -c bioconda bcftools
#conda install -c bioconda pyvcf