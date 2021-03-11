#' Create a coabundance object
#' @param edges tibble with columns from and to describing edges
#' @param method character of correlation method used
#' @param  max_pval maximal p-value. NULL to skip filtering
coabundance <- function(edges, nodes = NULL, method = NULL, max_pval = 0.05, min_abs_estimate = NULL, ...) {
  if (!is.null(max_pval) & "p.value" %in% colnames(edges)) {
    edges <-
      edges %>%
      filter(p.value <= max_pval)
  } else {
    warning("Ignore option max_pval: Not applicable")
  }
  
  if (!is.null(min_abs_estimate)) {
    edges <-
      edges %>%
      filter(abs(estimate) >= min_abs_estimate)
  }
  
  if (!is.null(nodes)) {
    features <- edges$from %>%
      union(edges$to) %>%
      unique()
    nodes <- nodes %>% dplyr::filter(feature %in% features)
  }
  
  edges <- edges %>% arrange(from, to)
  graph <- tidygraph::tbl_graph(edges = edges, nodes = nodes, directed = FALSE)
  
  res <-
    graph %>%

		# ensure simple graph
		tidygraph::activate(edges) %>%
		tidygraph::filter(! edge_is_multiple()) %>%
    
		tidygraph::to_undirected() %>%
    filter_graph(max_pval = max_pval, min_abs_estimate = min_abs_estimate) %>%
    topologize_graph() %>%
    annotate_node_attributes_in_edges()
  
  attr(res, "method") <- method
  res
}

coabundance.tbl_graph <- function(graph, nodes = NULL, method = NULL, max_pval = 0.05, min_abs_estimate = NULL, ...) {
  edges <- graph %>%
    annotate_node_attributes_in_edges() %>%
    activate(edges) %>%
    as_tibble() %>%
    select(from = from_feature, to = to_feature)
  coabundance(edges = edges, nodes = nodes, method = method, max_pval = max_pval, min_abs_estimate = min_abs_estimate, ...)
}

as_tbl_graph.spiec_easi_sparcc_res <- function(cor_res, ...) {
  taxa <- cor_res$boot$data %>% colnames()
  
  edges <-
    tidyr::expand_grid(from = taxa, to = taxa) %>%
    dplyr::mutate(comp = from %>% map2_chr(to, ~ c(.x, .y) %>%
                                             sort() %>%
                                             paste0(collapse = ""))) %>%
    dplyr::group_by(comp) %>%
    dplyr::slice(1) %>%
    dplyr::filter(from != to) %>%
    dplyr::ungroup() %>%
    dplyr::select(-comp) %>%
    dplyr::mutate(
      estimate = cor_res$pval$cors,
      p.value = cor_res$pval$pvals,
      q.value = p.adjust(p.value, method = "fdr")
    )
  
  coabundance(cor_res = cor_res, edges = edges, method = "sparcc", ...)
}

as_tbl_graph.rcorr <- function(cor_res, nodes = NULL, method = "rcorr", ...) {
  edges <-
    cor_res %>%
    broom::tidy() %>%
    dplyr::rename(from = column1, to = column2) %>%
    dplyr::mutate(q.value = p.adjust(p.value, method = "fdr"))
  
  coabundance(cor_res = cor_res, edges = edges, nodes = nodes, method = method, ...)
}

as_tbl_graph.tbl_df <- function(cor_res, nodes = NULL, method = NULL, ...) {
  if (!all(c("from", "to", "estimate") %in% colnames(cor_res))) {
    stop("cor_res must have at least columns from and to")
  }
  
  # keep only one row per comparison
  reduced_cor_res <-
    cor_res %>%
    dplyr::mutate(comp = from %>% map2_chr(to, ~ c(.x, .y) %>% sort() %>% paste0(collapse = "-"))) %>%
    dplyr::group_by(comp) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-comp)
  
  coabundance(edges = reduced_cor_res, nodes = nodes, method = method, ...)
}

as_tbl_graph.pulsar.refit <- function(cor_res, nodes = NULL, method = NULL, ...) {
  used_taxa <- cor_res$est$data %>% colnames()
  
  if (is.null(method)) {
    method <- cor_res$est$method
  }
  
  get_estimate <- switch(method, "mb" = SpiecEasi::getOptBeta, "glasso" = SpiecEasi::getOptCov)
  
  edges <-
    cor_res %>%
    get_estimate() %>%
    as.matrix() %>%
    as_tibble(rownames = "from") %>%
    pivot_longer(-one_of("from"), names_to = "to", values_to = "estimate") %>%
    mutate(to = to %>% str_remove("^V")) %>%
    readr::type_convert(col_types = cols(from = col_integer(), to = col_integer(), estimate = col_double()))
  
  graph <-
    cor_res %>%
    getRefit() %>%
    as.matrix() %>%
    adj2igraph() %>%
    as_tbl_graph() %>%
    mutate(feature = used_taxa) %>%
    activate(edges) %>%
    left_join(edges, by = c("from", "to")) %>%
    activate(nodes)
  
  cur_nodes <-
    graph %>%
    activate(nodes) %>%
    as_tibble()
  
  edges <-
    edges %>%
    left_join(cur_nodes %>% dplyr::rename(from_feature = feature), by = c("from" = "name")) %>%
    left_join(cur_nodes %>% dplyr::rename(to_feature = feature), by = c("to" = "name")) %>%
    dplyr::select(from = from_feature, to = to_feature, estimate) %>%
    dplyr::filter(estimate != 0 & from != to)
  
  coabundance(cor_res = cor_res, edges = edges, nodes = nodes, method = method, ...)
}

as_tbl_graph.coabundance <- function(x, nodes = NULL, method = NULL, ...) {
  if (is.null(method)) {
    method <- x$method
  }
  
  if (is.null(nodes)) {
    return(x)
  }
  
  edges <-
    x$graph %>%
    activate(edges) %>%
    as_tibble()
  
  nodes <-
    x$graph %>%
    activate(nodes) %>%
    as_tibble() %>%
    rename(feature = name) %>%
    left_join(nodes, by = "feature")
  
  
  coabundance(cor_res = x$result, edges = edges, nodes = nodes, method = method, ...)
}

as_tbl_graph.default <- function(cor_res, edges, nodes = NULL, method = NULL, ...) {
  if (!is.null(nodes)) {
    taxa <- edges$from %>%
      union(edges$to) %>%
      unique()
    nodes <- nodes %>% dplyr::filter(feature %in% taxa)
  }
  
  edges <- edges %>% arrange(from, to)
  graph <- tidygraph::tbl_graph(edges = edges, nodes = nodes, directed = FALSE)
  
  res <- list(graph = graph, result = cor_res, method = method)
  class(res) <- "coabundance"
  res %>% as_tbl_graph.coabundance(...)
}

#' Join two coabundances together to form an ensemble graph
#' 
#' @param x,y tbl_graphs to join. x wil be used to set the estimate
#' @param method joining method to combine the edges. Either 'union' or 'intersection'
ensemble_coabundance <- function(x, y, method = "intersection", name_x = "x", name_y = "y") {
  min_multiedges <- switch (method,
                            "intersection" = 2,
                            "union" = 1,
                            stop("Method not implemented.")
  )
  
  annotate_method <- function(x, name = "x") {
    state <- active(x)
    
    x %>%
      activate(edges) %>%
      mutate(method = name) %>%
      activate(!!sym(state)) # ensure idempotency
  }
  
  selected_method <- x %>% attr("method") %>% {.x <- .; ifelse(is.null(.x), "x", .x)}
  x <- x %>% annotate_method(name = name_x)
  y <- y %>% annotate_method(name = name_y)
  
  minimize_columns <-
    . %>%
    activate(edges) %>%
    select(any_of(c("from", "to", "p.value", "estimate", "method"))) %>%
    activate(nodes) %>%
    select(feature)
  
  # join graphs
  list(x, y) %>%
    map(minimize_columns) %>%
    reduce(graph_join) %>%
    annotate_node_attributes_in_edges() %>%
    
    # convert to tibble
    # required for pivoting
    activate(edges) %>%
    as_tibble() %>%
    select(from = from_feature, to = to_feature, estimate, p.value, method) %>%
    
    # only keep multiedges
    # more flexible than tidygraph::edge_is_multiple
    group_by(from, to) %>%
    mutate(n = n()) %>%
    arrange(from, to) %>%
    filter(n >= min_multiedges) %>%
    select(-n) %>%
    ungroup() %>%
    
    # polish
    pivot_wider(names_from = method, values_from = c(estimate, p.value)) %>%
    mutate(
      estimate := !!sym(paste0("estimate_", selected_method)),
      p.value := !!sym(paste0("p.value_", selected_method))
    ) %>%
    tidygraph::tbl_graph(edges = .) %>%
    activate(nodes) %>%
    rename(feature = name) %>%
    topologize_graph() %>%
    annotate_node_attributes_in_edges()
}

topologize_graph <- function(graph) {
  orig_state <- graph %>% active()
  
  graph <-
    graph %>%
    activate(nodes) %>%
    mutate(
      degree = centrality_degree(),
      component = group_components(),
      closeness = centrality_closeness(),
      betweeness = centrality_betweenness()
    )
  
  graph <- graph %>% activate(!!orig_state)
  
  graph
}

#' Filter a coabundance object
#'
#' @param x coabundance object
#' @param max_pval maximum p value to keep an edge. Will be ignores if the p.value is not available
#' @param min_abs_estimate minimal absolute value of the estimate to keep an edge. This is to filter edges with a low effect size.
#' @param remove_isolated_nodes TRUE if nodes not being part of any edge after filtering should be removed, FALSE otherwise.
#' @param recalculate_topology TRUE if node topology metrics e.g. centrality scores should be recalculated after filtering, FALSE otherwise.
filter_graph <- function(graph, max_pval = 0.05, min_abs_estimate = NULL, remove_isolated_nodes = TRUE, recalculate_topology = TRUE) {
  orig_state <- graph %>% active()
  
  edge_colnames <- graph %>%
    activate(edges) %>%
    as_tibble() %>%
    colnames()
  
  if (!is.null(max_pval)) {
    if (!"p.value" %in% edge_colnames) {
      warning("Ignore option max_pval: Not applicable")
    } else {
      graph <-
        graph %>%
        activate(edges) %>%
        tidygraph::filter(p.value <= max_pval)
    }
  }
  
  if (!is.null(min_abs_estimate)) {
    graph <-
      graph %>%
      activate(edges) %>%
      filter(abs(estimate) >= min_abs_estimate)
  }
  
  if (remove_isolated_nodes) {
    graph <-
      graph %>%
      activate(nodes) %>%
      filter(!node_is_isolated())
  }
  
  graph <-
    graph %>%
    activate(nodes) %>%
    filter(!tidygraph::node_is_isolated())
  
  graph <- graph %>% activate(!!orig_state)
  
  graph
}

# as_tbl_graph.tbl_graph <- function(x, cor_res = NULL, method = NULL) {
#   nodes <-
#     x %>%
#     tidygraph::activate(nodes) %>%
#     as_tibble()
#
#   edges <-
#     x %>%
#     tidygraph::activate(edges) %>%
#     as_tibble() %>%
#     dplyr::select(-from, -to) %>%
#     dplyr::rename(from = from_feature, to = to_feature)
#
#   coabundance(cor_res = cor_res, edges = edges, nodes = nodes, method = method)
# }

correlate_plot <- function(x) {
  # see http://www.fabiocrameri.ch/colourmaps.php
  vik_colors <- c(
    "#001261",
    "#033E7D",
    "#1E6F9D",
    "#71A8C4",
    "#C9DDE7",
    "#EACEBD",
    "#D39774",
    "#BE6533",
    "#8B2706",
    "#590008"
  )
  
  x %>%
    tidygraph::to_directed() %>% # needed for plotting
    ggraph::ggraph() +
    ggraph::geom_edge_link(aes(color = estimate)) +
    ggraph::geom_node_point() +
    ggraph::scale_edge_color_gradientn(colors = vik_colors, limits = c(-1, 1)) +
    ggplot2::theme_void()
}

annotate_node_attributes_in_edges <- function(graph) {
  state <- graph %>% active()
  
  graph %>%
    tidygraph::activate(edges) %>%
    # ensure idempotency
    tidygraph::select(-matches("^(from|to)_")) %>%
    tidygraph::left_join(
      graph %>%
        tidygraph::activate(nodes) %>%
        as_tibble() %>%
        mutate(name = row_number()) %>%
        rename_all(~ .x %>% paste0("from_", .)),
      by = c("from" = "from_name")
    ) %>%
    tidygraph::left_join(
      graph %>%
        tidygraph::activate(nodes) %>%
        as_tibble() %>%
        mutate(name = row_number()) %>%
        rename_all(~ .x %>% paste0("to_", .)),
      by = c("to" = "to_name")
    ) %>%
    activate(!!state)
}
