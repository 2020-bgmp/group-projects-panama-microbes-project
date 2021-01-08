#!/bin/bash

#Used to generate quality report of raw sequencing data

###load fastqc
module load fastqc/0.11.5

###create output directory
mkdir fqcOut

###run fastqc on all fastq files within current directory
fastqc -t 4 -o./fqcOut –noextract *.fastq.gz    

###load multiqc
ml MultiQC/1.3-Python-3.6.1

###compile into a multiqc report
###multiqc automatically recognizes .html files
multiqc .
