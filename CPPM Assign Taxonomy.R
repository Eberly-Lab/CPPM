
# Adapted from https://github.com/jkzorz/SituSeq

######
#parameters to set before running
subsample_depth = 1000 #each sample will be randomly subsampled to this number of reads, prior to taxonomic assignment (after filtering and trimming). For no subsampling see Nanopore_no_rarefaction.R under "backups" 
path_to_taxonomy_database = "silva_nr99_v138.1_train_set.fa.gz" #change to location of taxonomy database in relation to working directory (easiest to copy taxonomy database to working directory)
path_to_working_directory = "." #leave as a "." if you want to set your working directory manually in RStudio "Session"--> "Set Directory" --> "Choose Directory"
minBoot = 50 #Set the minBoot parameter for assignTaxonomy. minBoot refers to the minimum bootstrapping support required to return a taxonomic classification. Choose a number between 0-100, with 100 being the most stringent. 
######

#in R 
#set working directory 
setwd(path_to_working_directory)

#load packages in this order to avoid masking issues
library(ShortRead)
library(dada2)
library(tidyverse)

#save path to object
path = getwd()

#fastq filenames have format: 
#barcode01_combined.fastq 
fnFs = sort(list.files(path, pattern="_combined.fastq", full.names = TRUE))

#extract sample names, assuming filenames have format: #samplename_XXX.fastq
sample.names = sapply(strsplit(basename(fnFs), "\\."), `[`, 1)

#path for filtered and trimmed reads 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq"))
names(filtFs) = sample.names

#import sequences and assign taxonomy - with subsetting to subsampling depth
#this will create a csv file for each sample with the sequence and its assigned taxonomy
fastaRef <- "/Path/to/Silva_138.1/silva_nr99_v138.1_train_set.fa.gz"
for (fastq in filtFs) {
print(fastq)
seqs = getSequences(fastq)
sub = sample(1:length(seqs), subsample_depth, replace=FALSE) 
seq2 = seqs[sub]
tax_rc = assignTaxonomy(seq2, fastaRef, multithread=TRUE, tryRC = TRUE, minBoot = minBoot)
base = basename(fastq)
samples = gsub("_filt.fastq", "", base)
write.csv(tax_rc, paste('tax', samples, 'csv', sep = '.' ))
}

