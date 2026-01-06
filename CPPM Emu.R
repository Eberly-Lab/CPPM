setwd("set/path")


# Import and analyze taxonomic assignments of Nanopore sequence data performed in Emu
# https://github.com/treangenlab/emu
# Outputs from Emu are .tsv files with taxonomy and abundance
# The following workflow needs to be run for each location-year (sequencing run)

library(phyloseq)
library(microbiome)
library(dplyr)

df <- read.table(file = "emu-combined-tax_id-counts.tsv", sep = '\t', header = TRUE)
df <- data.frame(df)

# Removed unassigned taxa
df <- df[-c(1467),] # Adjust to match number of rows. This will be different for each sequencing run
colnames(df)

rownames(df) <- df$tax_id

# Add refseqs
library(Biostrings) 

fasta_file <- "../species_taxid.fasta"
fasta_sequences <- readDNAStringSet(fasta_file, format = "fasta")

fasta_tax_ids <- data.frame(names(fasta_sequences))

fasta_info <- data.frame(
  header = names(fasta_sequences),
  sequence = as.character(fasta_sequences),
  tax_id = NA,
  organism_name = NA,
  gene_name = NA,
  protein_name = NA,
  stringsAsFactors = FALSE
)

fasta_info$tax_id <- sapply(strsplit(fasta_info$header, ":"), `[`, 1)
fasta_info <- fasta_info[,-c(1,4:6)]

fasta_info <- fasta_info %>% distinct(tax_id, .keep_all = TRUE)

row.names(fasta_info) <- fasta_info$tax_id
row.names(df) <- df$tax_id

df2 <- as.data.frame(fasta_info)
df2 <- df2[,1, drop = FALSE]
df1 <- merge(df, df2, by = 0, all.x=TRUE)

df1 <- df1 %>% distinct(sequence, .keep_all = TRUE)
rownames(df1) <- df1$sequence
colnames(df1)
df1 <- df1[,-c(1,91)] # Adjust to match number of columns

# Clean up column names
df <- df1
names(df)
colnames(df)[1] <- 'TaxID'

tax <- df[,c(2:8)]
# Reorder and rename taxa columns
tax <- tax[,c(7,6,5,4,3,2,1)]
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

otu <- df[,c(9:89)]

########
# Replace NAs with 0 in otu table
otu[is.na(otu)] <- 0

str(otu)

sd <- read.csv("sample_data.csv") # e.g., "sample_data_BZ2022.csv"
rownames(sd) <- sd$Sample

# Make phyloseq object
taxmat <- as.matrix(tax)
otumat <- as.matrix(otu)
otu1 <- t(otu)

# Convert otu values to integers
otu_mat <- otu
otu_mat[,c(1:81)] <- lapply(otu_mat[,c(1:81)], as.integer) # Adjust to match number of columns

str(otu_mat)
class(otu_mat)
# Make a new sample data table to account for missing samples in otu table
df2 <- merge(sd, otu1, by = 0)
sd1 <- df2[,2:31]
row.names(sd1) <- sd1[,1]
sd1$Row.names <-NULL
sd1 <- sample_data(sd1)

OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
ps <- phyloseq(OTU, TAX, sd1)
ps

saveRDS(ps, "emu_phyloseq_bz_2022.RDS") 

# After running the above workflow for each location year, you will have 4 phyloseq objects:
# emu_phyloseq_bz_2022.RDS
# emu_phyloseq_bz_2023.RDS
# emu_phyloseq_sarc_2022.RDS
# emu_phyloseq_sarc_2023.RDS





