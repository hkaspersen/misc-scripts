#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

file_loc <- args[1]
output_dir <- args[2]
study <- args[3]
name <- args[4]
instrument <- args[5]
library_name <- args[6]
library_source <- args[7]
library_selection <- args[8]
description <- args[9]

# Libraries
library(dplyr)
library(tidyr)
library(tibble)

# Functions
func_paste <- function(x)
  paste(unique(x[!is.na(x)]), collapse = ",")


fix_samples <- function(df) {
  df <- df %>%
    select(-ref) %>%
    gather(field_name, field_value) %>%
    mutate(field_name = ifelse(
      grepl("fastq",
            field_name,
            ignore.case = T) == T,
      "FASTQ",
      field_name
    ))
  return(df)
}


# Create files
manifest_list <- read.table(
  file_loc,
  sep = "\t",
  header = F,
  stringsAsFactors = F
) %>%
  mutate(
    STUDY = study,
    NAME = name,
    INSTRUMENT = instrument,
    LIBRARY_NAME = library_name,
    LIBRARY_SOURCE = library_source,
    LIBRARY_SELECTION = library_selection,
    DESCRIPTION = description
  ) %>%
  mutate(ref = sub("_L00._R._001.fastq.gz", "", V1),
         SAMPLE = ref) %>%
  group_by(ref) %>%
  summarise_all(list(func_paste)) %>%
  separate(V1,
           into = c("FASTQ1",
                    "FASTQ2"),
           sep = ",") %>%
  select(
    ref,
    STUDY,
    SAMPLE,
    NAME,
    INSTRUMENT,
    LIBRARY_NAME,
    LIBRARY_SOURCE,
    LIBRARY_SELECTION,
    FASTQ1,
    FASTQ2
  ) %>%
  split(f = .$ref)

# Gather into two columns per file
file_list <- lapply(manifest_list,
                    function(x)
                      fix_samples(x))

# Write manifest files to disk
for (file in 1:length(file_list)) {
  write.table(
    file_list[[file]],
    paste0(output_dir, "/", names(file_list)[file], ".txt"),
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}
