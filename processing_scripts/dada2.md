## DADA2 Pipeline for fungal ITS sequencing data (forward reads only)
Adapted from https://benjjneb.github.io/dada2/ITS_workflow.html  

Important notes on the ITS region from https://benjjneb.github.io/dada2/ITS_workflow.html:  

"Unlike the 16S rRNA gene, the ITS region is highly variable in length. The commonly amplified ITS1 and ITS2 regions range from 200 - 600 bp in length. This length variation is biological, not technical, and arises from the high rates of insertions and deletions in the evolution of this less conserved gene region.  

The length variation of the ITS region has significant consequences for the filtering and trimming steps of the standard DADA2 workflow. First, truncation to a fixed length is no longer appropriate, as that approach remove real ITS variants with lengths shorter than the truncation length. Second, primer removal is complicated by the possibility of some, but not all, reads extending into the opposite primer when the amplified ITS region is shorter than the read length."  

## Preparation

Load packages and make note of versions:
```
library(dada2)
packageVersion("dada2")
library(Biostrings)
packageVersion("Biostrings")
library(ShortRead)
packageVersion("ShortRead")
```
Define path to point to the directory containing FASTQ files:  
```
path="./"
```
Parse out sample names from fastq file names. Assumes file names are of the format "SAMPLENAME_R1.fastq.gz"
```
name.fS <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
```
## Identify primers

Store ITS primer sequences:
```
FWD <- "CTTGGTCATTTAGAGGAAGTAA"
REV <- "GCTGCGTTCTTCATCGATGC"
```

Verify primers and their orientation in the data:
```
allOrients = function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
```
Pre-filter to remove sequences containing "N"
```
names.fS.filtN_fw <- file.path(path, "filtN_fw", basename(name.fS)) # path to put N filtered reads path+filtN_fw --forward reads

#filter N sequences from file
filterAndTrim(name.fS, names.fS.filtN_fw, maxN = 0, multithread = TRUE,compress=TRUE)

pathfiltN_fw="./filtN_fw/"
#list.files(pathfiltN_fw)
name.fS <- sort(list.files(pathfiltN_fw, pattern = "_R1.fastq.gz", full.names = TRUE))

names.fS.filtN_fw = file.path(path, "filtN_fw", basename(name.fS))
```
Count the number of times the primers appear in the reads:
```
primeHits <- function(primer, file) {
  hits <- vcountPattern(primer, sread(readFastq(file)), fixed = FALSE)
  return(sum(hits > 0))
}
rbind(FWD.in.FWD = sapply(FWD.orients, primeHits, file = names.fS.filtN_fw[[1]]), 
      REV.in.FWD = sapply(REV.orients, primeHits, file = names.fS.filtN_fw[[1]]))
```
## Remove primers

Load cutadapt:  
```
cutadapt <- "cutadapt" 
system2(cutadapt, args = "--version") #  system2 allows you to Run shell commands from R --will look for version of cutadapt
```
Create output file for trimmed reads,define parameters for cutadapt, and run:
```
cutPath <- file.path(path, "cutadapt_fw") 
if(!dir.exists(cutPath)) dir.create(cutPath) 
fs.cut <- file.path(cutPath, basename(name.fS))

FWD.RC <- dada2:::rc(FWD) #get the reverse compliment of fwd primer
REV.RC <- dada2:::rc(REV) #get reverse compliment of reverse primer

# set the flags for cutadapt to trim the FWD and the reverse-complement of the reverse primer off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# set the flags for cutadapt to trim the FWD and the reverse-complement of the reverse primer off of R2 (reverse reads)

# Run Cutadapt
for(i in seq_along(name.fS)) {
  system2(cutadapt, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fs.cut[i],  # output files
                             names.fS.filtN_fw[i])) # input files
}
```
Sanity check - Count the presence of primers after cutadapt removal: (Should be zero)
```
rbind(FWD.in.FWD = sapply(FWD.orients, primeHits, file = fs.cut[[1]]), 
      REV.in.FWD = sapply(REV.orients, primeHits, file = fs.cut[[1]]))
```
Run filter to move files that are too short to 2short directory:
```
system("mkdir cutadapt2Short_fw")
system("./cutadaptFilter_fw.sh") # cutadaptFilter.sh needs to be in the path directory from above.
```
Pull out sample names from cutadapt-ed files:
```
cutFs <- sort(list.files(cutPath, pattern = "_R1.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: McG.L2NLitter_R2.fastq
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
```
## Filter and trim

Assign file names for filtered output:
```
filtFWD <- file.path(cutPath, "filtered", basename(cutFs))
```
Filter and trim:
```
out <- filterAndTrim(cutFs, filtFWD, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = F, multithread = TRUE)  #filter the reads

saveRDS(out, "./out_2filttrim_fw.rds") # look at how many reads were filtered
```
