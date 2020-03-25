#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

report_loc <- args[1]
pattern <- args[2]
filter_value <- args[3]
organism <- args[4]
output_dir <- args[5]

filter_value <- as.numeric(filter_value)

# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(viridis)
library(impoRt)

# Functions
func_paste <- function(x) paste(unique(x[!is.na(x)]), collapse = ", ")

## Identifies the largest value in the cell and removes the other values
scan_max <- function(x) max(scan(text = x, 
                                 what = "",
                                 sep = ",",
                                 quiet = TRUE,
                                 strip.white = TRUE))

'%not_in%' <- Negate('%in%')

# Run analyses
print("Reading data...")

## Import data
mash_raw <- get_data(report_loc,
                     pattern,
                     delim = "\t",
                     col_names = FALSE) %>%
  rename("identity" = X1,
         "shared_hashes" = X2,
         "median_multiplicity" = X3,
         "p_value" = X4,
         "query_id" = X5,
         "query_comment" = X6) %>%
  mutate(p_value = round(as.numeric(p_value), 5),
         ref = sub("_mash.out", "", ref))

## Wrangle data
mash_results <- mash_raw %>%
  separate(shared_hashes, c("min","max"), sep = "/", remove = FALSE) %>%
  mutate(test = as.numeric(min)/as.numeric(max)) %>%
  select(-c(min, max)) %>%
  filter(test >= filter_value) %>%
  mutate(species = query_comment %>%
           str_replace("\\[.*?\\] ", "") %>%
           sub("^.+_.+\\.[0-9] (.*?), .+", "\\1", .) %>%
           word(1,2))

## Match species name for correct query
species_id <- mash_results %>%
  mutate(species_test = grepl(organism, species, ignore.case = TRUE)) %>%
  filter(species_test == TRUE) %>%
  select(species) %>%
  mutate(query = organism) %>%
  summarise_all(list(func_paste))

organism <- unlist(strsplit(species_id$species, ", ", fixed = TRUE))

print("Creating output files...")

# Plotting
mash_plot <- ggplot(mash_results, aes(ref, as.numeric(identity), fill = species))+
  geom_point(pch = 21, size = 2)+
  geom_hline(yintercept = 0.95, alpha = 0.5)+
  scale_fill_brewer(type = "div", palette = 2)+
  labs(y = "Identity",
       x = "Samples")+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.ticks.x = element_blank())

# Tables
contaminated_df <- mash_results %>%
  filter(species %not_in% organism,
         identity >= 0.95) %>%
  select(ref)
contaminated_ids <- contaminated_df$ref

if (length(contaminated_ids) == 0) {
  print("No contamined isolates were identified within set limits.")
} else {
  contam_ids <- suppressWarnings(
    mash_results %>%
      select(
        -c(
          shared_hashes,
          query_id,
          query_comment,
          median_multiplicity,
          test,
          p_value
        )
      ) %>%
      filter(ref %in% contaminated_ids) %>%
      mutate(id2 = 1:n()) %>%
      spread(species, identity, fill = NA) %>%
      select(-id2) %>%
      group_by(ref) %>%
      summarise_all(list(func_paste)) %>%
      mutate_at(vars(-ref), list(sapply(., scan_max))) %>%
      mutate_all(list(gsub("^$", NA, .)))
  )
  
  write.table(contam_ids,
              paste0(output_dir,
                     "/contaminated_samples_report.txt"),
              sep = "\t",
              row.names = FALSE)
}

## Create reports

mash_report <- suppressWarnings(
  mash_results %>%
    select(
      -c(
        shared_hashes,
        query_id,
        query_comment,
        median_multiplicity,
        test,
        p_value
      )
    ) %>%
    mutate(id2 = 1:n()) %>%
    spread(species, identity, fill = NA) %>%
    select(-id2) %>%
    group_by(ref) %>%
    summarise_all(funs(func_paste)) %>%
    mutate_at(vars(-ref), funs(sapply(., scan_max))) %>%
    mutate_all(funs(gsub("^$", NA, .)))
)

cont_species <- paste0(names(mash_report)[-1], collapse = ", ")
print(paste0("Species identified: ", cont_species))

# Write to file

write.table(mash_report,
            paste0(output_dir,
                   "/full_mash_report.txt"),
            sep = "\t",
            row.names = FALSE)

invisible(suppressGraphics(
  ggsave(paste0(output_dir,
                "/mash_plot.png"),
         mash_plot,
         device = "png",
         dpi = 100,
         height = 14,
         width = 16)
)
)

print("Analysis complete!")
