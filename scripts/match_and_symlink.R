#!/usr/bin/env Rscript

## Script that soft links files to a given directory,
## by matching files to ID's in a vector

args <- commandArgs(trailingOnly = TRUE)

file_list <- as.character(args[1])
file_dir <- args[2]
output_dir <- args[3]
pattern_setting <- args[4]
recursive_setting <- args[5]

# convert "recursive_setting" to a true TRUE/FALSE value
if (grepl("true", recursive_setting, ignore.case = TRUE) == TRUE) {
  recursive_setting <- TRUE
} else {
  recursive_setting <- FALSE
}

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
  n_unmatched <- c()
  total <- length(file_list)
  
  # match files to file list to identify full file names
  for (i in file_list) {
    for (j in total_files) {
      if (grepl(i,j) == TRUE) {
        final_files <- c(final_files, j)
      }
    }
  }
  
  # create matching string to identify how many entries in
  # file list were identified in the files
  final_files_string <- paste(final_files, collapse = " ")
  for (i in file_list) {
    if (grepl(i,final_files_string) == FALSE) {
      n_unmatched <- c(n_unmatched, i)
    }
  }
  
  # convert to number of unique entries
  n_unmatched <- length(unique(n_unmatched))
  
  # identify number of matched entries in file list
  matched <- total - n_unmatched
  
  # print results
  if (n_unmatched > 0) {
    print(paste0("Successfully matched ", matched, " of ", total, " ID strings."))
    print("Wrote output file 'unmatched_files.txt' with unmatched IDs.")
    
    unmatched_files <- total_files[which(final_files %in% total_files)]
    
    write.table(unmatched_files,
                "unmatched_files.txt",
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
  } else {
    print("Successfully identified all ID strings!")
  }
  return(final_files)
}

## soft-link matched files to output directory
symlink_files <- function(matched_files, output_dir) {
  for (i in matched_files) {
    file.symlink(paste0(file_dir, "/", i), output_dir)
  }
}

# import file list, tab separated file with a single column
keep_files <- read.table(file_list,
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
