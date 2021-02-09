#!/usr/bin/env Rscript

#
# Setup and configure drake 
#

source(here::here("src/setup.R"))

cache <- drake::new_cache(".drake")
drake::drake_cache(here::here(".drake"))$unlock()

make(
  # plan
  plan = plan,

  parallelism = "clustermq",
  
  # allow multiple works to access global env at the same time
  lock_envir = FALSE,
  
  recoverable = TRUE,
  recover = TRUE,
  
  # caching
  cache = cache,
  garbage_collection = TRUE,
  memory_strategy = "preclean",
  
  # config
  seed = 1337,
  keep_going = TRUE,
  jobs = jobs,
  
  # output
  verbose = 1,
  format = "qs",
  log_make = "log/drake.log"
)

rmarkdown::render(
  input = "src/report.Rmd",
  knit_root_dir = ".",
  output_dir = ".",
  output_file = "report.html",
)
