setwd("Your Directory")

# Adapted from https://github.com/jkzorz/SituSeq
#load packages (in this order to avoid masking issues)
library(ShortRead)
library(dada2)
library(tidyverse)

######


# Nanopore generates multiple fastq files for each barcode which need to be merged
# Run the following workflow prior to running Emu (https://github.com/treangenlab/emu)

# ---------------------------------
# Merge fastq.gz files per barcode 
# ---------------------------------

# --- 1) Set your top-level directory ---
# This is the root under which you want to search for fastq_pass directories from the nanopore sequencing run.
root_dir <- "path/to/fastq_pass"  # <-- change this

# --- 2) Output directory at the top level ---
out_root <- file.path(root_dir, "fastq_merged")
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

# --- 3) Find all directories named "fastq_pass" at any depth ---
all_dirs <- list.dirs(root_dir, recursive = TRUE, full.names = TRUE)
fastq_pass_dirs <- all_dirs[basename(all_dirs) == "fastq_pass"]

message(sprintf("Found %d fastq_pass directories.", length(fastq_pass_dirs)))
if (length(fastq_pass_dirs) == 0) {
  stop("No fastq_pass directories found under: ", root_dir)
}

# --- 4) Helper function ---
merge_fastq_gz <- function(input_files, output_file, chunk_size = 100000L) {
  if (length(input_files) == 0) {
    stop("No input FASTQ files provided to merge.")
  }
  # Ensure output directory exists
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Open output gz text connection
  con_out <- gzfile(output_file, open = "wt")

  on.exit({
    try(close(con_out), silent = TRUE)
  }, add = TRUE)

  total_lines <- 0L

  for (f in input_files) {
    # Read gz text connection in chunks
    con_in <- gzfile(f, open = "rt")
    repeat {
      lines <- readLines(con_in, n = chunk_size, warn = FALSE)
      if (length(lines) == 0) break
      writeLines(lines, con_out, sep = "\n", useBytes = TRUE)
      total_lines <- total_lines + length(lines)
    }
    close(con_in)
  }

  return(total_lines)
}

# --- 5) Build a map: barcode name -> all .fastq.gz files across all fastq_pass ---
barcode_files <- list()

for (fp in fastq_pass_dirs) {
  # Immediate subfolders like barcode01, barcode02, ...
  subdirs <- list.dirs(fp, recursive = FALSE, full.names = TRUE)

  # Match 'barcode' followed by digits (e.g., barcode01, barcode12, barcode001)
  barcode_dirs <- subdirs[grepl("^barcode\\d+$", basename(subdirs))]
  if (length(barcode_dirs) == 0) next

  for (bd in barcode_dirs) {
    bc <- basename(bd)

    # Files directly in barcode dir (typical ONT layout)
    fastqs_direct <- list.files(bd, pattern = "\\.fastq\\.gz$", full.names = TRUE, recursive = FALSE)
    # Also include nested files if present (some pipelines nest)
    fastqs_nested <- list.files(bd, pattern = "\\.fastq\\.gz$", full.names = TRUE, recursive = TRUE)

    fastqs <- sort(unique(c(fastqs_direct, fastqs_nested)))
    if (length(fastqs) == 0) {
      message(sprintf("[INFO] No .fastq.gz files found under %s", bd))
      next
    }

    if (is.null(barcode_files[[bc]])) barcode_files[[bc]] <- character(0)
    barcode_files[[bc]] <- c(barcode_files[[bc]], fastqs)
  }
}

# Deduplicate and order paths deterministically for each barcode
for (nm in names(barcode_files)) {
  barcode_files[[nm]] <- sort(unique(barcode_files[[nm]]))
}

unique_barcodes <- sort(names(barcode_files))
message(sprintf("Discovered %d unique barcodes: %s",
                length(unique_barcodes),
                if (length(unique_barcodes) > 0) paste(unique_barcodes, collapse = ", ") else "NONE"))

if (length(unique_barcodes) == 0) {
  stop("No barcode directories with .fastq.gz files found under fastq_pass.")
}

# --- 6) Merge per barcode into one output at top-level fastq_merged ---
results <- data.frame(
  barcode       = character(),
  files_merged  = integer(),
  lines_written = integer(),
  output_file   = character(),
  stringsAsFactors = FALSE
)

for (bc in unique_barcodes) {
  files <- barcode_files[[bc]]
  if (length(files) == 0) next

  message(sprintf("Preparing %s: %d files", bc, length(files)))
  message(paste("  Example inputs:", paste(head(files, 3), collapse = " | ")))

  out_file <- file.path(out_root, paste0(bc, "_merged.fastq.gz"))

  # Avoid creating empty outputs: only call helper if we have inputs
  lines_written <- merge_fastq_gz(files, out_file, chunk_size = 100000L)

  # Quick size check
  sz <- file.info(out_file)$size
  message(sprintf("  -> Wrote %s (%s bytes, %s lines)", basename(out_file),
                  format(sz, big.mark = ","),
                  format(lines_written, big.mark = ",")))

  results <- rbind(
    results,
    data.frame(
      barcode       = bc,
      files_merged  = length(files),
      lines_written = lines_written,
      output_file   = out_file,
      stringsAsFactors = FALSE
    )
  )
}

# Use FastQC to check read quality and get total read count
setwd("path/to/fastq")

library(fastqcr)
dir.create("fastqc")
fastqc(fq.dir = "fastq_combined", # FASTQ files directory
       qc.dir = "fastqc", # Results direcory
       threads = 10                    # Number of threads
       )
qc <- qc_aggregate("fastqc")
summary(qc)

count_df <- subset(qc, module == "Basic Statistics")
count_df$tot.seq <- as.integer(count_df$tot.seq)
sum(count_df$tot.seq) 



# Quality filtering with dada2
minLength = 1200 #Removes reads shorter than this length. Minimum length is enforced AFTER trimming
maxLength = 1800 #Removes reads longer than this length. Maximum length is enforced BEFORE trimming 
trimLeft = 100 #The number of nucleotides to remove from the start of each read. Must cover the primer length
trimRight = 100 #The number of nucleotides to remove from the end of each read. Must cover the primer length

# 11.27.24 Tried adjusting min and max length to get more reads
minLength = 300
maxLength = 4000

######
#in R

# Nanopore generates multiple fastq files for each barcode
# Combine all fastq files for each barcode
setwd("path/to/fastq_pass")
folders <- list.files(pattern = "barcode" )

for (directory in folders) {
  print(directory) 
  files = list.files(path = paste(directory, "/", sep = ""), pattern = ".fastq*") 
  print(files)
  fout = file.path(paste(directory, "combined.fastq.gz", sep = "_"))
  for (fl in files) {
    fq = readFastq(paste(directory,"/",fl, sep = ""))
    writeFastq(fq, fout, mode="a")
  }
}

#save path to object
path = 'path/to/fastq_pass'

# Combined fastq filenames have format: 
#barcode01_combined.fastq 
fnFs = sort(list.files(path, pattern="_combined.fastq", full.names = TRUE))

#Extract sample names, assuming filenames have format: #samplename_XXX.fastq
sample.names = sapply(strsplit(basename(fnFs), "\\."), `[`, 1)
sample.names = sapply(strsplit(sample.names, "\\_"), `[`, 1)

#filter and trim reads- create new paths for new files 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq"))
names(filtFs) = sample.names

#filter and trim command - will create new fastq files with filtered and trimmed sequences
out = filterAndTrim(fnFs, filtFs, trimLeft = trimLeft, trimRight = trimRight, maxLen = maxLength, minLen = minLength,  truncQ = 0, compress = FALSE)

#see how many reads were lost and write to csv file
head(out,12) 
write.csv(out, "Filtered_sequence_summary.csv")

