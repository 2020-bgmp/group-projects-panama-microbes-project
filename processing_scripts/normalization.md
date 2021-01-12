# Normalization with Phyloseq in R

Using information in the "sample_sums" output from DADA2, a cutoff is chosen to rarify at an even depth. This is done by creating a phyloseq object and then using the phyloseq function "rarefy_even_depth". 

Rarefying is a classic method of normalizing that chooses a sample size and randomly subsamples that number of sequences from each sample. NOTE: The default phyloseq function samples with replacement.

This was done on TALAPAS, the University of Oregon's supercomputer.

Create phyloseq object:
```
library(dplyr)
library(phyloseq)
library(DESeq2)
library(vegan)

metadata <- read_tsv("Panama_metadata_full_v1.txt")

#contains OTU IDs and functional groups
funguild <- read.table("./nodrysamples/transformed_ASVfunguildtable.guilds.txt", row.names=1, header=T, sep='\t')

#contains OTU Ids and sample names
otus <- read.table("./nodrysamples/PanamaPrecip_ITS_ASVs.txt", row.names=1, header=T)

otus$otuID <- rownames(otus)
fung_sub <- funguild %>% select(Guild)
fung_sub$otuID <- rownames(fung_sub)

#combine taxonomy and guild data for future use in analysis
otu_taxa <- merge(otus, fung_sub, by=0, all=TRUE)
rownames(otu_taxa) <- otu_taxa$Row.names
otu_taxa <- otu_taxa %>% select(!c(Row.names, otuID.y, otuID.x))

#create taxonomy matrix
taxa <- as.matrix(otu_taxa[,1:7])
taxa <- cbind(taxa, otu_taxa$Guild)
colnames(taxa)[8] <- c('Guild')

#subset metadata to include only variables of interest
meta_sub <- subset(metadata, select = c("#sampleID", "PlotCreated", "pH.water", "ResinP.mg_kg", "MAP", "Latitude", "Longitude"))

#format dataframe to match other dataframes for merging
names(meta_sub)[1] <- "Sample_ID"
meta_sub$Sample_ID <- gsub("_", "\\.", meta_sub$Sample_ID)

#create ASV dataframe
otumap <- otus[,-c(1:7)] #check column names and make sure you're pulling out the ASV id column through all taxonomy columns 
otumap2 <- otumap %>% select(one_of(as.character(meta_sub$Sample_ID))) #check the dimensions to make sure you're not losing any samples from metadata sample ID errors
sample_names <- colnames(otumap2)

#more formatting to make sure everything matches
otut <- data.frame(t(otumap2))
rownames(meta_sub) <- meta_sub$Sample_ID
meta2 <- meta_sub[rownames(otut),] #just wanted to get everything in the same order
rownames(meta2) <- meta2$Sample_ID

#make a phyloseq object
OTU <- otu_table(otumap2, taxa_are_rows = TRUE) 
TAX <- tax_table(taxa)
physeq <- phyloseq(OTU,TAX)
sampledata <- sample_data(as.data.frame(meta2,row.names=sample_names(physeq), stringsAsFactors=FALSE) %>% drop_na())
ps <- phyloseq(OTU, sampledata, TAX)
```
Rarefy to even depth and save output: 
```
set.seed(111) #make sure it's reproducible
###changed sample size here to Krista's request of 10,000
physeq2 <- rarefy_even_depth(physeq1, sample.size = 10000) 
saveRDS(physeq2, "./rarefied_physeq.rds") #save the phyloseq object to make downstream analyses easier

rarefied_taxa <- data.frame(tax_table(physeq2))
ot <- data.frame(otu_table(physeq2))
rarefy_tab <- cbind(rarefied_taxa,ot)

write.table(rarefy_tab, file="PanamaPrecip_ITS_ASVs_r10000.txt", sep="\t", row.names=FALSE, quote=FALSE)
```
