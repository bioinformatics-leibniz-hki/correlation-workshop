#!/usr/bin/env Rscript

#
#  Setup script
#

# laod libraries 
# The order matters: put most important
# packages at the bottom to overwrite the namespace
library(curatedMetagenomicData)
library(phyloseq)
library(SpiecEasi)
library(patchwork)
library(ggraph)
library(drake)
library(tidygraph)

library(tidyverse)
library(magrittr)
library(dplyr)

filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate

# source plans
c("src/lib", "src/plans/") %>%
  list.files(pattern = "\\.R$", full.names = TRUE) %>%
  walk(source)

# merge all objects ending with _plan to the master plan
plan <-
  ls() %>%
  keep(~ .x %>% str_ends("_plan")) %>%
  map(~ .x %>% parse(text = .) %>% eval()) %>%
  bind_plans()

# Number of parallel threads
jobs <- as.integer(system("nproc", intern = T))

options(
  mc.cores = jobs,
  clustermq.scheduler = "multicore"
)

# Set directory for external R packages not included in the install script
dir.create("/analysis/R-site-library")
.libPaths(c("/analysis/R-site-library/", .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths()[1], Sys.getenv("R_LIBS"), sep = .Platform$path.sep))

# make fastspar available
Sys.setenv(PATH = paste0("/miniconda3/envs/fastspar/bin/", ":", Sys.getenv("PATH")))

#
# project specific theme setup
#

# Set default theme
theme_my <- function() {
  ggplot2::theme_minimal(base_size = 20) +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(size = 0.8),
      axis.line.y = ggplot2::element_line(size = 0.8),
      axis.ticks = ggplot2::element_line(colour = "black", size = 0.8),
      axis.text = ggplot2::element_text(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = element_blank()
    )
}

ggplot2::theme_set(theme_my())

# Feature type colors
feature_type_colors <- c(
  taxon = "#e76f51",
  pathway = "#2a9d8f"
)

scale_fill_feature_type <- partial(
  ggplot2::scale_fill_manual,
  values = feature_type_colors
)

scale_color_feature_type <- partial(
  ggplot2::scale_color_manual,
  values = feature_type_colors
)

# disease colors
disease_colors <- c(
  healthy = "#90be6d",
  CRC = "#f9844a"
)

scale_fill_disease <- partial(
  ggplot2::scale_fill_manual,
  values = disease_colors
)

scale_color_disease <- partial(
  ggplot2::scale_color_manual,
  values = disease_colors
)
