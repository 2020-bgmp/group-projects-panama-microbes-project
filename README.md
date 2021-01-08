# Panama-microbes-project

Pipeline for processing and analysis of fungal ITS amplicon sequencing data of forward reads only. Used to analyze soil samples from across the isthumus of Panama for the McGuire lab at the University of Oregon.   

<img width="610" alt="workflow" src="https://user-images.githubusercontent.com/54604213/104035462-03bbb480-5187-11eb-9791-fca049770839.png">

## Part (a) - data proccessing &#8594; taxonomy assignment/normalization

All scripts run on TALAPAS, the University of Oregon supercomputer.

**Tools/versions:**  
FASTQC - 0.11.5  
QIIME - 1.9.1  
DADA2 - 1.10.1  
ShortRead - 1.40.0  
BioStrings - 2.48.0  
Cutadapt - 1.15  
Funguild - 1.1  

## Part (b) - hierarchical modelling/diversity analysis

All scripts run on local machines in R.

**Packages/versions:**  
tidyverse - 1.3.0  
lme4 - 1.1.26  
DESeq2 - 1.28.1  
phyloseq - 1.32.0  
vegan - 2.5.7  
goeveg - 0.4.2  
MASS - 7.3.53  

 

 


