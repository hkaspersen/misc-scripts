#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

report_loc <- args[1]
organism <- args[2]
output_dir <- args[3]

# Libraries

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,
               dplyr,
               tidyr,
               stringr,
               purrr,
               viridis,
               svglite,
               R.devices)

# Functions
func_paste <- function(x) paste(unique(x[!is.na(x)]), collapse = ", ")

## Identifies the largest value in the cell and removes the other values
scan_max <- function(x) max(scan(text = x, 
                                 what = "",
                                 sep = ",",
                                 quiet = TRUE,
                                 strip.white = TRUE))

## Identifies filenames in input folder
file_names_mash <- function(filepath) {
  files <- list.files(path = filepath, pattern = "_mash.out")
  return(files)
}

## Import ariba data from report.tsv from chosen database used in ariba
get_mash_data <- function(filepath) {
  files <- file_names_mash(filepath)
  
  data_list <- lapply(files,
                      FUN = function(file) {
                        read.delim(
                          paste0(filepath, "/", file),
                          stringsAsFactors = F,
                          header = FALSE, # no header in files
                          quote = "", # disable quoting, avoid EOS error
                          sep = "\t"
                        )
                      })
  
  names(data_list) <- files # set name of each df in list 
  data <- bind_rows(data_list, .id = "ref") %>% # bind to data frame
    rename("identity" = V1,
           "shared_hashes" = V2,
           "median_multiplicity" = V3,
           "p_value" = V4,
           "query_id" = V5,
           "query_comment" = V6) %>%
    mutate(p_value = round(p_value, 5))
  return(data)
}

# Run analyses

print("Reading data...")

## Import data
mash_raw <- get_mash_data(report_loc) %>%
  mutate(ref = sub("_L00[0-9]_R[0-9]_00[0-9].fastq.gz_mash.out", "", ref))

## Wrangle data
mash_results <- mash_raw %>%
  separate(shared_hashes, c("min","max"), sep = "/", remove = FALSE) %>%
  mutate(test = as.numeric(min)/as.numeric(max)) %>%
  select(-c(min, max)) %>%
  filter(test >= 0.1) %>%
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
  summarise_all(funs(func_paste))

organism <- species_id$species

print("Creating output files...")

# Plotting
mash_plot <- ggplot(mash_results, aes(ref, identity, fill = species))+
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
  filter(species != organism,
         identity >= 0.95) %>%
  select(ref)
contaminated_ids <- contaminated_df$ref

## Create reports
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
    summarise_all(funs(func_paste)) %>%
    mutate_at(vars(-ref), funs(sapply(., scan_max))) %>%
    mutate_all(funs(gsub("^$", NA, .)))
)

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
write.table(contam_ids,
            paste0(output_dir,
                   "/contaminated_samples_report.txt"),
            sep = "\t",
            row.names = FALSE)

write.table(mash_report,
            paste0(output_dir,
                   "/full_mash_report.txt"),
            sep = "\t",
            row.names = FALSE)

invisible(suppressGraphics(
  ggsave(paste0(output_dir,
                "/mash_plot.svg"),
         mash_plot,
         device = "svg",
         dpi = 100,
         height = 14,
         width = 16)
  )
)

print("Analysis complete!")