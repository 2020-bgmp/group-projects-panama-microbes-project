#!/bin/bash

#Used to demultiplex samples using QIIME

###changer header lines of barcode files to be compatible with QIIME
sed 's/2:N:0/1:N:0/g' TAAGGCGA_S1_L001_R2_001.fastq.gz > TAAGGCGAbc_fix.fastq.gz
sed 's/2:N:0/1:N:0/g' CGTACTAG_S2_L001_R2_001.fastq.gz > CGTACTAGbc_fix.fastq.gz

###correct mapping files to be compatible with QIIME
sed 's/\_/\./g' map_file_nx1.txt > map_file_nx1.corrected.txt
sed 's/\_/\./g' map_file_nx2.txt > map_file_nx2.corrected.txt

###load easybuild and QIIME
ml easybuild  icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132 
ml QIIME

###check that mapping files are correctly formatted
validate_mapping_file.py -m map_file_nx1.corrected.txt
validate_mapping_file.py -m map_file_nx2.corrected.txt

###create directories for demultiplexed output
mkdir nx1f_demult
mkdir nx1r_demult
mkdir nx2f_demult
mkdir nx2r_demult





######run QIIME command using SBATCH, can split up forward and reverse commands to decrease time. Each set takes <30 minutes######

unset I_MPI_PMI_LIBRARY #QIIME has issues with MPI setting in easybuild, so unset it

#demultiplex forward reads
/usr/bin/time -v \
split_libraries_fastq.py -o nx2f_demult/ -i data/CGTACTAG_S2_L001_R1_001.fastq.gz\
 -b data/CGTACTAGbc_fix.fastq.gz -m data/map_file_nx2.corrected.txt\
 --rev_comp_mapping_barcode -q 0 -r 999 -n 999 -p 0.0001 --store_demultiplexed_fastq
 
 #demultiplex reverse reads
 /usr/bin/time -v \
split_libraries_fastq.py -o nx2r_demult/ -i data/CGTACTAG_S2_L001_R3_001.fastq.gz\
 -b data/CGTACTAGbc_fix.fastq.gz -m data/map_file_nx2.corrected.txt\
 --rev_comp_mapping_barcode -q 0 -r 999 -n 999 -p 0.0001 --store_demultiplexed_fastq
 
##############################################







######run these commands within each output demultiplexed directory (eg. nx1f_demult, nx1r_demult, etc.)######

###check number of reads to make sure they match for forward and reverse reads
grep -c -e '^@[A-Z].*M0' seqs.fastq

###create directory for files split by sample
mkdir splitseqs

###split files on sample ID using QIIME
split_sequence_file_on_sample_ids.py -i seqs.fastq --file_type fastq -o splitseqs/ 

###rename files to inclue R1 and R2 in respective directories 
rename .fastq _R1.fastq *
rename .fastq _R2.fastq *

###in main directory, create directory to hold all files
mkdir combinedseqs

###back in each demultiplex directory, mv files to combined directory
mv * ../combinedseqs

##############################################




###in combinedseqs directory, compress all files
gzip *
