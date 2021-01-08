## DADA2 Pipeline for fungal ITS sequencing data
Adapted from https://benjjneb.github.io/dada2/ITS_workflow.html  



```
library(dada2)
packageVersion("dada2")
library(Biostrings)
packageVersion("Biostrings")
library(ShortRead)
packageVersion("ShortRead")

path="./"

name.fS <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))

FWD = "CTTGGTCATTTAGAGGAAGTAA"
REV = "GCTGCGTTCTTCATCGATGC"

allOrients = function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients = allOrients(FWD)
REV.orients = allOrients(REV)
FWD.orients

names.fS.filtN_fw = file.path(path, "filtN_fw", basename(name.fS)) # path to put N filtered reads path+filtN_fw --forward reads

#filter N sequences from file
filterAndTrim(name.fS, names.fS.filtN_fw, maxN = 0, multithread = TRUE,compress=TRUE)

pathfiltN_fw="./filtN_fw/"
#list.files(pathfiltN_fw)
name.fS <- sort(list.files(pathfiltN_fw, pattern = "_R1.fastq.gz", full.names = TRUE))

names.fS.filtN_fw = file.path(path, "filtN_fw", basename(name.fS))

primeHits <- function(primer, file) {
  hits <- vcountPattern(primer, sread(readFastq(file)), fixed = FALSE)
  return(sum(hits > 0))
}
rbind(FWD.in.FWD = sapply(FWD.orients, primeHits, file = names.fS.filtN_fw[[1]]), 
      REV.in.FWD = sapply(REV.orients, primeHits, file = names.fS.filtN_fw[[1]]))

cutadapt <- "cutadapt" 
system2(cutadapt, args = "--version") #  system2 allows you to Run shell commands from R --will look for version of cutadapt

cutPath <- file.path(path, "cutadapt_fw") # make a file for cutadapt output
if(!dir.exists(cutPath)) dir.create(cutPath) # create the folder for cutadapt output
fs.cut <- file.path(cutPath, basename(name.fS)) #put forward reads in cutadapt output folder

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

rbind(FWD.in.FWD = sapply(FWD.orients, primeHits, file = fs.cut[[1]]), 
      REV.in.FWD = sapply(REV.orients, primeHits, file = fs.cut[[1]]))

system("mkdir cutadapt2Short_fw")
system("./cutadaptFilter_fw.sh") # cutadaptFilter.sh needs to be in the path directory from above.
`

cutFs <- sort(list.files(cutPath, pattern = "_R1.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: McG.L2NLitter_R2.fastq
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

filtFWD <- file.path(cutPath, "filtered", basename(cutFs))


out <- filterAndTrim(cutFs, filtFWD, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = F, multithread = TRUE)  #filter the reads

saveRDS(out, "./out_2filttrim_fw.rds") # look at how many reads were filtered
```
