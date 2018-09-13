#!/usr/bin/env Rscript

## Script that soft links files to a given directory,
## by matching files in a list of file names.

args <- commandArgs(trailingOnly = TRUE)

file_list <- args[1]
file_dir <- args[2]
output_dir <- args[3]
pattern_setting <- args[4]
recursive_setting <- args[5]

# functions
## get file names in file_dir
get_files <- function(file_dir, pattern = "fastq.gz", recursive = FALSE) {
  if (recursive == FALSE) {
    total_files <- list.files(file_dir, pattern = pattern)
  }
  if (recursive == TRUE) {
    total_files <- list.files(file_dir, pattern = pattern, recursive = TRUE)
  }
  return(total_files)
}

## match file names in two vectors
match_files <- function(file_list, total_files) {
  final_files <- c()
  for (i in file_list) {
    for (j in total_files) {
      if (grepl(i,j) == TRUE) {
        final_files <- c(final_files, j)
      }
    }
  }
  return(final_files)
}

## soft-link matched files to output directory
symlink_files <- function(matched_files, output_dir) {
  for (i in matched_files) {
    file.symlink(paste0(file_dir, "\\", i), output_dir)
  }
}

# import file list, tab separated file with a single column
keep_files <- read.table(file_list,
                    header = TRUE,
                    sep = "\t",
                    stringsAsFactors = FALSE)

keep_files <- keep_files[,]

# get file names in file_dir
total_files <- get_files(file_dir,
                         pattern = pattern_setting,
                         recursive = recursive_setting)

# match file names in file_list to total_files
matched_files <- match_files(keep_files, total_files)

# symlink files from matched_files to output_dir
symlink_files(matched_files, output_dir)
