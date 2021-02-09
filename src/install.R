#!/usr/bin/env Rscript

#
# Installing R packages inside the docker image
#

options(
  Ncpus = as.integer(system("nproc", intern = T))
)

devtools::install_cran(
  pkgs = c(
    "drake",
    "here",
    "qs",
    "clustermq",
    "Hmisc",
    "ppcor",
    "visNetwork",
    "tidygraph",
    "ggpubr",
    "styler",
    "patchwork",
    "ggraph"
  ),
  upgrade = "never"
)

devtools::install_bioc(
  c(
    "3.11/curatedMetagenomicData",
    "3.11/phyloseq",
    "3.11/ComplexHeatmap",
    "3.11/banocc"
  ),
  upgrade = "never"
)

devtools::install_github(
  "zdk123/SpiecEasi",
  ref = "master",
  upgrade = "always"
)
