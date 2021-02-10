#!/usr/bin/env Rscript

#
# Create plots for the report
#

plots_plan <- drake_plan(
  node_topology_metrics = c("degree", "closeness", "betweeness"),
  
  abundance_heatmap_plot = {
    abundances_tbl <-
      subset_abundances %>%
      unnest(abundance) %>%

      # normalization
      group_by(sample_id, feature_type) %>%
      mutate(abundance = scale(abundance)) %>%
      
      # long to wide table
      ungroup() %>%
      select(feature, sample_id, abundance) %>%
      pivot_wider(
        names_from = feature,
        values_from = abundance,
        values_fill = list(abundance = 0)
      )
    
    abundance_mat <-
      abundances_tbl %>%
      ungroup() %>%
      select(-sample_id) %>%
      as.matrix()
    rownames(abundance_mat) <- abundances_tbl$sample_id
    
    lineages_tbl <-
      tibble(feature = colnames(abundance_mat)) %>%
      filter(feature != "sample_id") %>%
      left_join(lineages, by = "feature")
    
    samples_tbl <-
      tibble(sample_id = rownames(abundance_mat)) %>%
      left_join(samples, by = "sample_id")
    
    viridis_col <- circlize::colorRamp2(c(-5, 0, 5), c("#440154FF", "#21908CFF", "#FDE725FF"))
    
    ComplexHeatmap::Heatmap(
      matrix = abundance_mat %>% t(),
      name = "Abundance (z-score))",
      
      col = viridis_col,
      column_split = samples$disease,
      row_split = lineages_tbl$feature_type,

      left_annotation = ComplexHeatmap::HeatmapAnnotation(
        which = "row",
        Phylum = lineages_tbl$phylum
      ),
      top_annotation = ComplexHeatmap::HeatmapAnnotation(
        which = "column",
        Disease = samples_tbl$disease,
        Age = samples_tbl$age,
        Country = samples_tbl$country,
        Gender = samples_tbl$gender
      )
    )
  },

  # plot node topology for each project and feature_type
  intra_feature_type_topology_plots = target(
    {
      res <-
        intra_feature_type_coabundances %>%
        filter(method == intra_feature_type_cor_methods) %>%

        # get one graph containing all edges
        mutate(
          graph = list(coabundance, disease, feature_type) %>% pmap(~ {
            ..1 %>%
              activate(nodes) %>%
              mutate(disease = ..2, feature_type = ..3)
          })
        ) %>%
        pull(graph) %>%
        reduce(tidygraph::graph_join) %>%

        # get nodes table
        activate(nodes) %>%
        as_tibble() %>%
        select(c(feature, feature_type, disease, node_topology_metrics)) %>%
        pivot_longer(cols = node_topology_metrics, names_to = "metric", values_to = "value") %>%

        # plot
        ggplot(aes(disease, value, fill = disease)) +
        geom_boxplot() +
        facet_wrap(feature_type ~ metric, scales = "free_y") +
        scale_fill_disease() +
        labs(
          # title = "Node topology of coabundance network",
          # subtitle = str_glue("Method: {intra_feature_type_cor_methods}"),
          x = "Feature type",
          fill = "Disease",
          y = "Topology"
        )

      tibble(
        method = intra_feature_type_cor_methods,
        plot = list(res)
      )
    },
    dynamic = map(intra_feature_type_cor_methods)
  ),

  intra_feature_type_network_plots = target(
    {
       res <- switch (intra_feature_type_coabundances$feature_type[[1]],
         "taxon" = {
            intra_feature_type_coabundances$coabundance[[1]] %>%
               activate(nodes) %>%
               left_join(lineages) %>%
               to_directed() %>%
               ggraph(layout = "igraph", algorithm = "fr", weights = abs(estimate)) +
               geom_node_point(aes(color = phylum)) +
               # geom_node_text(aes(label = phylum), size = 3) +
               geom_edge_link(aes(color = log10(estimate))) +
               scale_edge_color_gradient2(low = "blue", mid = "transparent", high = "red", midpoint = 0, limits = c(-1, 1)) +
               theme_void()
         },
         "pathway" = {
            intra_feature_type_coabundances$coabundance[[1]] %>%
               activate(nodes) %>%
               to_directed() %>%
               ggraph(layout = "igraph", algorithm = "fr", weights = abs(estimate)) +
               geom_node_point() +
               geom_edge_link(aes(color = estimate)) +
               scale_edge_color_gradient2(low = "blue", mid = "transparent", high = "red", midpoint = 0, limits = c(-1, 1)) +
               theme_void()
         }
       )

      intra_feature_type_coabundances %>%
        discard(is.list) %>%
        mutate(
          plot = list(res)
        )
    },
    dynamic = map(intra_feature_type_coabundances)
  )
)
