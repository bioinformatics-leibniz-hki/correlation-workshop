correlate <- function(data, method = "spearman", ...) {
  switch(method,
    "sparcc" = correlate_sparcc(data = data, ...),
    "mb" = correlate_spiec_easi_pulsar(data = data, method = "mb", ...),
    "glasso" = correlate_spiec_easi_pulsar(data = data, method = "glasso", ...),
    "pearson" = correlate_pearson(data = data, ...),
    "spearman" = correlate_spearman(data = data, ...),
    "banocc" = correlate_banocc(data = data, ...),
    stop(stringr::str_glue("method {method} is not implemented!"))
  )
}

correlate_banocc <- function(
  data,
  compiled_banocc_model = NULL,
  conf_alpha = 0.05,
  n = rep(0, ncol(data)),
  L = 10 * diag(ncol(data)),
  a = 0.5,
  b = 0.01,
  threads = getOption("mc.cores", 1L),
  chains = 4,
  iter = 50,
  warmup = floor(iter / 2),
  thin = 1,
  init = NULL,
  control = NULL,
  verbose = FALSE,
  num_level = 0, ...) {
  if (is.null(compiled_banocc_model)) {
    compiled_banocc_model <- rstan::stan_model(model_code = banocc::banocc_model)
  }

  b_fit <- banocc::run_banocc(
    C = data,
    compiled_banocc_model = compiled_banocc_model,
    n = n,
    L = L,
    a = a,
    b = b,
    cores = threads,
    chains = chains,
    iter = iter,
    warmup = warmup,
    thin = thin,
    init = init,
    control = control,
    verbose = verbose,
    num_level = num_level
    
  )

  b_output <- banocc::get_banocc_output(banoccfit = b_fit, conf_alpha = conf_alpha)

  b_output %>%
    purrr::pluck("Estimates.median") %>%
    tibble::as_tibble(rownames = "from") %>%
    tidyr::pivot_longer(-from, names_to = "to", values_to = "estimate") %>%
    dplyr::mutate(conf_alpha = conf_alpha) %>%
    tidyr::replace_na(list(estimate = 0)) %>%
    dplyr::filter(estimate != 0) %>%
    as_tbl_graph(...)
}

#' Coabundance analysis using SparCC as implemented in fastspar
#'
#' This implementation is way faster than `correlate_spiec_easi_sparcc` but requires linux and the external shell command `fastspar`.
#' @param data integer matrix of abundance count data. One sample per row and one taxon per column
correlate_fastspar <- function(data, iterations = 50, exclude_iterations = 10, bootstraps = 200, threads = getOption("mc.cores"), ...) {
  system <- function(...) base::system(ignore.stdout = TRUE, ignore.stderr = TRUE, ...)

  threads <- min(parallel::detectCores(), threads)

  # sanity checks
  if (class(data) != "matrix") stop("data must be of type matrix")
  if (str_glue("fastspar --version") %>% system() != 0) {
    stop("Command fastspar not found")
  }

  dir <- tempfile(pattern = "fastspar")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  data_path <- str_glue("{dir}/data.tsv")

  data %>%
    t() %>%
    as_tibble(rownames = "#OTU ID") %>%
    write_tsv(data_path)

  paste(
    "fastspar",
    "--yes",
    "--iterations", iterations,
    "--exclude_iterations", exclude_iterations,
    "--otu_table", data_path,
    "--correlation", paste0(dir, "/cor.tsv"),
    "--covariance", paste0(dir, "median_covariance.tsv"),
    "--threads", threads,
    sep = " "
  ) %>%
    system()


  paste0(dir, "/bootstraps_counts") %>% dir.create()
  paste0(dir, "/bootstraps_cor") %>% dir.create()

  paste(
    "fastspar_bootstrap",
    "--otu_table", data_path,
    "--number", bootstraps,
    "--prefix", paste0(dir, "/bootstraps_counts", "/data"),
    sep = " "
  ) %>%
    system()

  paste(
    "parallel",
    "--jobs", threads,

    "fastspar",
    "--yes",
    "--otu_table {}",
    "--correlation", paste0(dir, "/bootstraps_cor/cor_{/}"),
    "--covariance", paste0(dir, "/bootstraps_cor/cov_{/}"),
    "--iterations", iterations,
    "--exclude_iterations", exclude_iterations,
    ":::",
    paste0(dir, "/bootstraps_counts/*"),
    sep = " "
  ) %>%
    system()

  paste(
    "fastspar_pvalues",
    "--otu_table", data_path,
    "--correlation", paste0(dir, "/cor.tsv"),
    "--prefix", paste0(dir, "/bootstraps_cor/cor_data_"),
    "--permutations", bootstraps,
    "--outfile", paste0(dir, "/pvals.tsv"),
    sep = " "
  ) %>%
    system()

  pval_tbl <-
    paste0(dir, "/pvals.tsv") %>%
    read_tsv(col_types = cols(.default = "c")) %>%
    dplyr::rename(from = `#OTU ID`) %>%
    tidyr::pivot_longer(-from, names_to = "to", values_to = "p.value")

  cor_tbl <-
    paste0(dir, "/cor.tsv") %>%
    readr::read_tsv(col_types = cols(.default = "c")) %>%
    dplyr::rename(from = `#OTU ID`) %>%
    tidyr::pivot_longer(-from, names_to = "to", values_to = "estimate")

  paste0("rm -rf ", dir) %>% system()

  pval_tbl %>%
    dplyr::inner_join(cor_tbl, by = c("from", "to")) %>%
    # only keep triangle
    dplyr::mutate(comp = from %>% map2_chr(to, ~ c(.x, .y) %>%
      sort() %>%
      paste0(collapse = "-"))) %>%
    dplyr::group_by(comp) %>%
    dplyr::slice(1) %>%
    dplyr::select(-comp) %>%
    readr::type_convert() %>%
    dplyr::ungroup() %>%
    dplyr::filter(from != to) %>%
    dplyr::mutate(q.value = p.adjust(p.value, method = "fdr")) %>%
    dplyr::select(from, to, p.value, q.value, estimate) %>%
    as_tbl_graph(...)
}

#' Coabundance analysis using SparCC
#' @param implementation Character indicating the implementation of SparCC algorithm to use. One of "fastspar", "spiec_easi"
#' @param ... further arguments passed to the sparcc correlation functions
correlate_sparcc <- function(implementation = "spiec_easi", ...) {
  switch(implementation,
    "fastspar" = correlate_fastspar(...),
    "spiec_easi" = correlate_spiec_easi_sparcc(...),
    stop(stringr::str_glue("{implementation} must be one of fastspar, auto, or spiec_easi"))
  )
}

correlate_spiec_easi_pulsar <- function(data, method = "mb", nlambda = 300, sel.criterion = "stars", pulsar.select = TRUE, pulsar.params = list(
                                          thresh = 0.05,
                                          subsample.ratio = 0.8,
                                          ncores = getOption("mc.cores"),
                                          rep.num = 20,
                                          seed = 1337
                                        ), ...) {
  spiec.easi(
    data = data, method = method, nlambda = nlambda, sel.criterion = sel.criterion,
    pulsar.select = pulsar.select, pulsar.params = pulsar.params
  ) %>%
    as_tbl_graph(...)
}


#' Coabundance analysis using SparCC as implemented in SpiecEasi
#' @param data integer matrix of abundance count data. One sample per row and one taxon per column
correlate_spiec_easi_sparcc <- function(data, iterations = 10, bootstraps = 200, threads = getOption("mc.cores"), ...) {
  sparcc_boot <-
    SpiecEasi::sparccboot(
      data = data,
      R = bootstraps,
      ncpus = threads,
      sparcc.params = list(
        iter = iterations,
        inner_iter = iterations,
        th = 0.1
      )
    )

  sparcc_pval <-
    sparcc_boot %>%
    SpiecEasi::pval.sparccboot(sided = "both")

  list(boot = sparcc_boot, pval = sparcc_pval) %>%
    structure(class = "spiec_easi_sparcc_res") %>%
    as_tbl_graph(...)
}

correlate_spiec_easi <- function(
                                 data,
                                 method = "glasso",
                                 sel.criterion = "stars",
                                 verbose = TRUE,
                                 pulsar.select = TRUE,
                                 pulsar.params = list(
                                   thresh = 0.05,
                                   subsample.ratio = 0.8,
                                   ncores = getOption("mc.cores"),
                                   rep.num = 20,
                                   seed = 1337
                                 ),
                                 icov.select = pulsar.select,
                                 icov.select.params = pulsar.params,
                                 lambda.log = TRUE,
                                 ...) {
  SpiecEasi::spiec.easi(
    data = data,
    method = method,
    sel.criterion = sel.criterion,
    verbose = verbose,
    pulsar.select = pulsar.select,
    pulsar.params = pulsar.params,
    icov.select = icov.select,
    icov.select.params = icov.select.params,
    lambda.log = lambda.log
  ) %>% as_tbl_graph(...)
}

correlate_pearson <- function(data, ...) {
  data %>%
    Hmisc::rcorr(type = "pearson") %>%
    as_tbl_graph(...)
}


correlate_spearman <- function(data, ...) {
  data %>%
    Hmisc::rcorr(type = "spearman") %>%
    as_tbl_graph(...)
}
