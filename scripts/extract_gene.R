#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

file_loc <- args[1]
gene_name <- args[2]
nt_or_aa <- args[3]
output_dir <- args[4]

gene_name2 <- gsub(" |,", "_", gene_name)

# Libraries

if (!require("pacman")) install.packages("pacman")
pacman::p_load(Biostrings,
               seqinr)

## Functions
# Identifies filenames in input folder
file_names <- function(filepath, pattern) {
  files <- list.files(path = filepath,
                      pattern = pattern,
                      recursive = TRUE)
  return(files)
}

# Import data from assembled_genes.fa.gz files created by ariba
get_seq_data <- function(filepath, type = "aa") {
  
  if (type == "aa") {
    files <- file_names(filepath, ".faa")
    
    data_list <- lapply(files,
                        FUN = function(file) {
                          readAAStringSet(paste0(filepath, "/", file))
                        })
  }
  
  if (type == "nt") {
    files <- file_names(filepath, ".ffn")
    
    data_list <- lapply(files,
                        FUN = function(file) {
                          readDNAStringSet(paste0(filepath, "/", file))
                        })
  }

  names(data_list) <- files
  return(data_list)
}

# Filter out specific gene of interest
filter_gene <- function(fasta_file, gene) {
  gene_pos <- which(grepl(gene, names(fasta_file), ignore.case = TRUE) == TRUE)
  gene_seq <- fasta_file[gene_pos]
  return(gene_seq)
}

# Append sequences to fasta file
append_fasta_file <- function(list, gene) {
  i <- 1
  for (x in list) {
    write.fasta(x, names(list[i]), paste0(output_dir, "/", gene, "_", nt_or_aa, "_sequences.fa"), open = "a")
    i <- i + 1
  }
}

## Run functions
all_sequences <- get_seq_data(file_loc, nt_or_aa)
gene_sequences <- lapply(all_sequences, function(x) filter_gene(x, gene_name))
append_fasta_file(gene_sequences, gene_name2)
