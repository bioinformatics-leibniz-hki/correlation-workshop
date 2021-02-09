#!/usr/bin/env Rscript

#
# Get data about feature abundance and samples
#

data_plan <- drake_plan(
  # projects present in column dataset_name in table curatedMetagenomicData::combined_metadata
  projects = c("NielsenHB_2014", "ZellerG_2014"),
  
  # disease names present in column disease in table curatedMetagenomicData::combined_metadata
  diseases = c("CRC", "healthy"),
  
  # Works only if metaphlan data is available
  feature_types = c("taxon", "pathway"),

  # Subset: number of features to draw per subset (e.g. disease)
  n_features = 20,

  # Subset: number of samples to draw per subset (e.g. disease)
  n_samples = 20,

  # Min percentage of samples the taxon must be present to be included
  min_prevalence_perc = 10,

  # Sample meta data table
  samples = {
    curatedMetagenomicData::combined_metadata %>%
      as_tibble() %>%
      filter(
        disease %in% diseases &
          dataset_name %in% projects
      ) %>%
      group_by(disease) %>%
      sample_n(n_samples) %>%
      ungroup() %>%
      transmute(
        sample_id = sprintf("S%02i", row_number()),
        sample_alias = sampleID,
        project = dataset_name,
        body_site, disease, age, age_category, gender, country
      ) %>%
      mutate(across(where(is.character), as.factor))
  },
  
  # Abundance table in long format
  all_abundances = target(
    {
      curated_metagenomic_data_feature_type <-
        feature_types %>% recode(taxon = "metaphlan_bugs_list", pathway = "pathabundance_relab")
      
      raw_abundance <-
        paste(projects, curated_metagenomic_data_feature_type, "stool()", sep = ".") %>%
        map(~ eval(parse(text = .x))) %>%
        pluck(1) %>%
        exprs() %>%
        as_tibble(rownames = "feature") %>%
        gather(sample_alias, abundance, -feature) %>%
        # subset and annotate samples
        # inner join to only get subset of samples
        inner_join(samples, by = "sample_alias") %>%
        mutate(feature_type = feature_types)

      # clean names
      abundance <-
        switch(feature_types,
          "taxon" = {
            raw_abundance %>%
              # only keep metaphlan species counts
              filter(
                feature %>% str_detect("s__[A-Za-z_]+$")
              ) %>%
              # cleaning names
              mutate(
                lineage = feature,
                feature = feature %>%
                  str_extract("s__[A-Za-z_]+$") %>%
                  str_remove("^s__") %>%
                  str_replace_all("_", " ")
              )
          },
          "pathway" = {
            raw_abundance %>%
              # discard species contribution of the pathway
              filter(! feature %>% str_detect("\\|")) %>%
              # discard unannotated pathways
              filter(! feature %in% c("UNMAPPED", "UNINTEGRATED")) %>%
              mutate(feature = feature %>% str_remove("^[A-z-+0-9]+: "))
          }
        )

      # no nested table here, just one very long one
      # e.g. just return the abundance tibble
      # meta data is not required because the sample_id is unique
      abundance
    },
    dynamic = cross(feature_types, projects)
  ),
  
  lineages = {
    all_abundances %>%
      distinct(feature_type, feature, lineage) %>%
      mutate(lineage = lineage %>% str_remove_all("[a-z]__") %>% str_replace_all("_", " ")) %>%
      separate(
        col = lineage,
        sep = "\\|",
        into = c("kingdom", "phylum", "order", "class", "family", "genus", "species")
      )
  },
  
  subset_abundances = {
    all_abundances %>%
      group_by(feature_type, disease) %>%
      nest() %>%
      rename(abundance = data) %>%
      mutate(
        abundance = abundance %>% map(~ {
          subset_samples <- .x$sample_id %>% unique()
          
          prevalent_features <-
            .x %>%
            filter(abundance > 0) %>%
            # count samples per feature
            distinct(sample_id, feature) %>%
            group_by(feature) %>%
            count() %>%
            # filter by prevalence
            filter(n / length(subset_samples) * 100 >= min_prevalence_perc) %>%
            arrange(-n) %>%
            pull(feature)
          
          if(length(prevalent_features) < n_features) {
            str_glue(
              "Not enough prevalent taxa! Found {length(prevalent_features)} ",
              "but {n_features} were requested."
            ) %>%
              stop()
          }
          
          selected_features <-
            prevalent_features %>%
            head(n_features)
          
          .x %>% filter(feature %in% selected_features)
        }
        )
      )
  },
  
  write_subsets = {
    dir.create("data", showWarnings = FALSE)
    
    all_abundances %>%
      select(sample_id, feature, abundance) %>%
      pivot_wider(names_from = feature, values_from = abundance) %>%
      write_csv(file_out("data/features.csv"))
    
    samples %>%
      write_csv(file_out("data/samples.csv"))
    
    lineages %>%
      write_csv(file_out("data/lineages.csv"))
  }
)
