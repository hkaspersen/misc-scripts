#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

report_loc <- args[1]
output_dir <- args[2]

# Libraries

install.packages("pacman")
pacman::p_load(fastqcr,
               dplyr,
               ggplot2,
               tidyr,
               viridis,
               ggsci,
               scales,
               svglite,
               R.utils) 

# Functions

# Function that creates the output directory
check_dir <- function(output_dir) {
  folder <- paste0("fastqc_results_", Sys.Date())
  dir.create(file.path(output_dir, folder), showWarnings = FALSE)
}

# Function that lists all the file names in the input folder
file_names <- function(filepath, folder) {
  files <- list.files(path = paste0(filepath, "/", folder), pattern = "_fastqc.zip")
  return(files)
}

# Function that searches recursively for all filenames with fastqc.zip
file_names_recursive <- function(filepath) {
  files <- list.files(path = filepath, pattern = "_fastqc.zip", recursive = TRUE)
  return(files)
}

# Lists folders in the input folder
folder_names <- function(filepath) {
  folders <- list.files(path = filepath)
  return(folders)
}

# Identifies the names of the files and groups them accoring to their respective folders
get_grouping_variable <- function(path, folder) {
  files <- file_names(path, folder)
  df <- data.frame(files = files, group = folder, stringsAsFactors = FALSE)
  df$files <- gsub("(.*?)_fastqc.zip", "\\1", df$files)
  colnames(df) <- c("ref","group")
  return(df)
}

# Creates a data frame with grouping information of the reports
create_group_df <- function(path) {
  folders <- folder_names(path)
  if (length(folders) > 1) {
    df <- lapply(folders, function(folder) get_grouping_variable(path, folder))
    df <- bind_rows(df)
  } else {
    df <- get_grouping_variable(path, folders)
  }
  return(df)
}

# Function that filter out duplicated counts
filter_counts <- function(df) {
  df <- df %>%
    mutate(dupl = duplicated(Count)) %>%
    filter(dupl == FALSE)
  return(df)
}

# Function that prepares the sequence length data for plotting
prepare_seq_len_data <- function(list) {
  x <- split(list$sequence_length_distribution, list$sequence_length_distribution$group)
  x <- lapply(x, filter_counts)
  x <- bind_rows(x)
  return(x)
}

# Function that imports and wrangles fastqc data
get_fastqc_data <- function(filepath) {
  folders <- folder_names(filepath)
  
  get_files <- file_names_recursive(filepath)
  
  data_list <- lapply(get_files,
                      FUN = function(file) {
                        qc_read(paste0(filepath, "/", file),
                                modules = "all",
                                verbose = FALSE)
                      })
  names(data_list) <- gsub("(.*?)/(.*?)_fastqc.zip", "\\2", get_files)
  data_list <- purrr::transpose(data_list)
  
  data_list$sequence_length_distribution <- NULL
  data_list$kmer_content <- NULL
  
  list_names <- names(data_list)
  list_numbers <- 1:length(list_names)
  for (i in list_numbers) {
    assign(list_names[i], bind_rows(data_list[[i]], .id = "ref"))
  }
  df_list <- list(summary,
                  basic_statistics,
                  per_base_sequence_quality,
                  per_tile_sequence_quality,
                  per_sequence_quality_scores,
                  per_base_sequence_content,
                  per_sequence_gc_content,
                  per_base_n_content,
                  sequence_duplication_levels,
                  overrepresented_sequences,
                  adapter_content,
                  total_deduplicated_percentage)
  names(df_list) <- list_names
  df_list$basic_statistics <- df_list$basic_statistics %>%
    spread(Measure,Value) %>%
    left_join(group_df, by = "ref")
  return(df_list)
}

# Function that saves plots
save_plots <- function(plot, title, height, width) {
  ggsave(paste0(output_dir, "/", title, ".svg"),
         plot,
         dpi = 100,
         device = "svg",
         units = "cm",
         height = height,
         width = width)
}

# Function that creates and saves plots from fastqc data
create_plots <- function(df_list) {
  p1 <- df_list$adapter_content %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length", "group")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    gather(key,
           value, -c(ref,
                     Position,
                     group,
                     seqlen)) %>%
    ggplot(aes(factor(
      Position,
      levels = unique(Position),
      ordered = TRUE
    ), value, color = key)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(outlier.size = 0.5) +
    labs(
      x = "Position in Read",
      y = "Percent (%) Adapter Content",
      color = NULL,
      title = "Adapter content"
    ) +
    scale_colour_jama() +
    scale_y_continuous(limits = c(0, 20)) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom"
    ) +
    facet_wrap(~ group, scales = "free", dir = "v")
  
  p2 <- df_list$per_base_sequence_content %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length", "group")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    gather(key, value, -c(ref, Base, group, seqlen)) %>%
    ggplot(aes(factor(
      Base, levels = unique(Base), ordered = TRUE
    ), value, color = key)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(outlier.size = 0.5) +
    labs(
      x = "Position in Read",
      y = "Percent (%)",
      color = NULL,
      title = "Per base sequence content"
    ) +
    theme_classic() +
    scale_color_jama() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom"
    ) +
    facet_wrap( ~ group, scales = "free", dir = "v")
  
  p3 <- df_list$sequence_duplication_levels %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length", "group")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    gather(key, value, -c(ref, `Duplication Level`, group, seqlen)) %>%
    ggplot(aes(
      factor(
        `Duplication Level`,
        levels = unique(`Duplication Level`),
        ordered = TRUE
      ),
      value,
      fill = key
    )) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = c("#ef8a62",
                                 "#67a9cf")) +
    theme_classic() +
    labs(x = "Duplication Level",
         y = "Percent (%) of Sequences",
         fill = NULL,
         title = "Sequence duplication levels") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.4
          )) +
    facet_wrap( ~ group, scales = "free")
  
  # p4 <- df_list %>%
  #   prepare_seq_len_data() %>%
  #   left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
  #   rename(seqlen = "Sequence length") %>%
  #   ggplot(aes(factor(Length,
  #                     ordered = TRUE,
  #                     levels = unique(Length)),
  #              Count)) +
  #   geom_boxplot(outlier.size = 0.5) +
  #   scale_y_continuous(labels = comma) +
  #   scale_fill_viridis(discrete = TRUE) +
  #   labs(x = "Read Size",
  #        y = "Total sequence length",
  #        title = "Sequence length per read size") +
  #   guides(fill = FALSE) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(
  #     angle = 90,
  #     hjust = 1,
  #     vjust = 0.4
  #   )) +
  #   facet_wrap(~seqlen, scales = "free")
  
  p5 <- df_list$per_sequence_quality_scores %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length", "group")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    ggplot(aes(factor(Quality), Count, fill = factor(Quality))) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(outlier.size = 0.5) +
    scale_y_continuous(labels = comma) +
    scale_x_discrete(breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40)) +
    scale_fill_viridis(discrete = TRUE) +
    labs(x = "Quality",
         y = "Number of reads",
         title = "Per sequence quality scores") +
    guides(fill = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 12)) +
    facet_wrap( ~ group)
  
  p6 <- df_list$per_sequence_gc_content %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length", "group")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    ggplot(aes(factor(`GC Content`), Count)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(outlier.size = 0.5) +
    labs(x = "GC content (%)",
         y = "Number of reads",
         title = "Per sequence GC content") +
    scale_y_continuous(labels = comma) +
    scale_x_discrete(breaks = as.character(seq(
      from = 0, to = 100, by = 10
    ))) +
    theme_classic() +
    facet_wrap( ~ group, scales = "free")
  
  p7 <- df_list$per_base_n_content %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length", "group")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    ggplot(aes(factor(
      Base, levels = unique(Base), ordered = TRUE
    ), `N-Count`)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(fill = "#e6e6e6",
                 outlier.size = 0.5) +
    labs(x = "Position in read",
         title = "Per base N content") +
    guides(fill = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    facet_wrap( ~ group, scales = "free", dir = "v")
  
  p8 <- df_list$total_deduplicated_percentage %>%
    gather(key = "ref", value = `1`) %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length", "group")], by = "ref") %>%
    rename(seqlen = "Sequence length",
           perc = `1`) %>%
    mutate(perc = as.numeric(perc)) %>%
    ggplot(aes(factor(group), perc)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(fill = "#e6e6e6",
                 outlier.size = 0.5) + 
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = "Group",
         y = "Total percentage of deduplicated reads",
         title = "Total deduplicated percentage") +
    theme_classic()
  
  # This produces a huge figure - excluded for now
  # p9 <- df_list$per_tile_sequence_quality %>%
  #   left_join(., df_list$basic_statistics[, c("ref", "Sequence length")], by = "ref") %>%
  #   rename(seqlen = "Sequence length") %>%
  #   ggplot(aes(factor(Base,
  #                     levels = unique(Base),
  #                     ordered = TRUE),
  #              factor(Tile),
  #              fill = Mean)) +
  #   geom_tile() +
  #   scale_fill_viridis(limits = c(-10, 5),
  #                      breaks = c(-10, -8, -6, -4, -2, 0, 2, 4),
  #                      guide=guide_colourbar(ticks=T,nbin=50,barheight=.5,label=T,barwidth=10)) +
  #   labs(y = "Tiles",
  #        x = "Position in read",
  #        title = "Per tile sequence quality (mean)",
  #        fill = NULL) +
  #   theme_classic() +
  #   theme(legend.position="bottom",
  #         legend.justification="center",
  #         legend.direction="horizontal",
  #         legend.text=element_text(color="grey20"),
  #         axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         axis.text.y = element_text(size = 6)) +
  #   facet_wrap(~seqlen, scales = "free")
  
  p10 <- df_list$per_base_sequence_quality %>%
    left_join(., df_list$basic_statistics[, c("ref", "Sequence length", "group")], by = "ref") %>%
    rename(seqlen = "Sequence length") %>%
    mutate(Base = factor(Base,
                         levels = unique(Base),
                         ordered = TRUE)) %>%
    group_by(group) %>%
    mutate(xmax = length(unique(Base)) + 1) %>%
    ungroup() %>%
    ggplot(aes(Base,
               Mean)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(outlier.size = 0.4,
                 fill = "#7f7f7f") +
    geom_hline(aes(yintercept = 28),
               color = "green") +
    labs(x = "Position in read",
         y = "Sequence quality",
         title = "Per base mean sequence quality") +
    scale_y_continuous(limits = c(0, 42)) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    facet_wrap(~ group, scales = "free", dir = "v")
  
  save_plots(p1, "adapter_content", 25, 35)
  save_plots(p2, "per_base_sequence_content", 25, 35)
  save_plots(p3, "sequence_duplication_levels", 25, 30)
  save_plots(p5, "per_sequence_quality_scores", 25, 35)
  save_plots(p6, "per_sequence_gc_content", 25, 30)
  save_plots(p7, "per_base_n_content", 25, 35)
  save_plots(p8, "total_deduplicated_percentage", 20, 20)
  save_plots(p10, "per_base_mean_sequence_quality", 25, 25)
}

# Create grouping data frame
group_df <- create_group_df(report_loc)

# Check output directory
check_dir(output_dir)
output_dir <- paste0(output_dir, "/fastqc_results_", Sys.Date())

# Data and plotting
df_list <- get_fastqc_data(report_loc)
create_plots(df_list)
