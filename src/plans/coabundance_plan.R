#!/usr/bin/env Rscript

#
# Coabundance analysis of abundance profiles
#

coabundance_plan <- drake_plan(
  intra_feature_type_cor_methods = c("spearman", "sparcc", "banocc"),

  max_pval = 0.05,
  min_abs_estimate = 0.1,

  # parameters of function SpiecEasi::spiec.easi
  spiec_easi_params = list(
    method = "mb",
    nlambda = 30,
    lambda.min.ratio = 1e-2,
    sel.criterion = "stars",
    pulsar.select = TRUE,
    pulsar.params = list(
      thresh = 0.05,
      subsample.ratio = 0.8,
      ncores = 32,
      rep.num = 20,
      seed = 1337
    )
  ),

  sparcc_params = list(
    iterations = 10,
    bootstraps = 100
  ),

  banocc_params = list(
    compiled_banocc_model = compiled_banocc_model,
    conf_alpha = 0.05,
    chains = 5,
    iter = 10e3,
    warmup = 5e3,
    thin = 10
  ),

  compiled_banocc_model = rstan::stan_model(model_code = banocc::banocc_model),

  # List of abundance matrices needed for coabundance functions
  # within feature_type (e.g. only taxa or pathways)
  intra_feature_type_mats = target(
    {
      res <-
        subset_abundances %>%
        pull(abundance) %>%
        first() %>%
        select(sample_id, feature, abundance) %>%
        pivot_wider(
          names_from = feature,
          values_from = abundance,
          values_fill = list(abundance = 0)
        ) %>%
        set_rownames(.$sample_id) %>%
        select(-sample_id) %>%
        as.matrix()

      subset_abundances %>%
        discard(is.list) %>%
        mutate(mat = list(res))
    },
    dynamic = map(subset_abundances)
  ),

  # List of abundance matrices needed for banocc coabundance function
  # within feature_type (e.g. only taxa or pathways)
  intra_feature_type_banocc_mats = target(
    {
      res <-
        subset_abundances %>%
        pull(abundance) %>%
        first() %>%
        select(sample_id, feature, abundance) %>%

        # banocc requires abundances summed to 1 for each sample
        group_by(sample_id) %>%
        mutate(abundance = abundance / sum(abundance)) %>%
        ungroup() %>%
        pivot_wider(
          names_from = feature,
          values_from = abundance,
          values_fill = list(abundance = 0)
        ) %>%
        set_rownames(.$sample_id) %>%
        select(-sample_id) %>%
        as.matrix()

      subset_abundances %>%
        discard(is.list) %>%
        mutate(mat = list(res))
    },
    dynamic = map(subset_abundances)
  ),

  # List of abundance matrices needed for coabundance functions
  # between feature_types (e.g. both taxa and pathways)
  inter_feature_type_mats = target(
    {
      subset_abundances <-
        subset_abundances %>%
        filter(disease == diseases) %>%
        unnest(abundance)

      # get samples of subset group having positive abundance
      # in both feature types (both taxa and pathway counts)
      subset_samples <-
        subset_abundances %>%
        filter(abundance > 0) %>%
        distinct(feature_type, sample_id) %>%
        group_by(sample_id) %>%
        count() %>%
        filter(n == 2) %>%
        pull(sample_id)

      res <-
        subset_abundances %>%
        filter(sample_id %in% subset_samples) %>%
        group_by(feature_type) %>%
        nest() %>%
        transmute(
          # long tibble to matrix for each feature_type
          mat = data %>% map(~ {
            res <-
              .x %>%
              select(sample_id, feature, abundance) %>%
              pivot_wider(
                names_from = feature,
                values_from = abundance,
                values_fill = list(abundance = 0)
              ) %>%
              select(-sample_id) %>%
              as.matrix()

            rownames(res) <- .x$sample_id %>% unique()

            res
          })
        ) %>%
        pull(mat)

      tibble(
        disease = diseases,
        mat = list(res)
      )
    },
    dynamic = map(diseases)
  ),

  # calcualte coabundances between features of different type
  # e.g. a coabundance between a pathway and a taxon
  inter_feature_type_coabundances = target(
    {
      cancel_if(inter_feature_type_mats$mat[[1]] %>% length() != 2)

      res <-
        spiec_easi_params %>%
        inset2("data", inter_feature_type_mats$mat[[1]]) %>%
        inset2("method", "mb") %>%
        inset2("nodes", lineages) %>%
        do.call(what = correlate)

      inter_feature_type_mats %>%
        discard(is.list) %>%
        mutate(
          method = spiec_easi_params$method,
          feature_type = "all",
          coabundance = list(res)
        )
    },
    dynamic = map(inter_feature_type_mats),
    # use internal parallelization instead
    hpc = FALSE
  ),

  # calculate coabundances within features of a feature type
  # e.g. taxon - taxon coabundances
  intra_feature_type_coabundances = target(
    {
      params <- switch(intra_feature_type_cor_methods,
        "sparcc" = {
          sparcc_params %>%
            inset2("data", intra_feature_type_mats$mat[[1]]) %>%
            inset2("method", "sparcc")
        },
        "banocc" = {
          banocc_params %>%
            inset2("data", {
              # use corresponding banocc matrix
              # because sample sum must be 1.
              # Doing this with inner joins
              intra_feature_type_banocc_mats %>%
                rename(banocc_mat = mat) %>%
                inner_join(intra_feature_type_mats) %>%
                pull(banocc_mat) %>%
                first()
            })
        },
        {
          list(
            data = intra_feature_type_mats$mat[[1]],
            method = intra_feature_type_cor_methods
          )
        }
      ) %>%
        inset2("max_pval", max_pval) %>%
        inset2("nodes", lineages) %>%
        inset2("min_abs_estimate", min_abs_estimate)

      res <-
        params %>%
        do.call(what = correlate)

      intra_feature_type_mats %>%
        discard(is.list) %>%
        mutate(method = intra_feature_type_cor_methods, coabundance = list(res))
    },
    dynamic = cross(intra_feature_type_mats, intra_feature_type_cor_methods),
    # use internal parallelization instead
    hpc = FALSE
  ),

  # all coabundances combined
  coabundances = {
    bind_rows(
      inter_feature_type_coabundances,
      intra_feature_type_coabundances
    )
  }
)
