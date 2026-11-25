#' Run the full analysis for a group
#'
#' Convenience wrapper that runs activity budgets, species-typical behavioural
#' diversity, the spread of participation index, zone electivity, an
#' association network and (if the ethogram has directed behaviours) an
#' interaction network, for one sampling method.
#'
#' @inheritParams behaviour_diversity
#' @inheritParams spi
#' @param animals Optional animal attributes table (column `animal`).
#' @param ... Passed to [association_matrix()] (e.g. `method`, `index`,
#'   `max_dist`, `map_size`).
#' @return A list of class `ethogram_analysis`.
#' @examples
#' res <- analyse_group(example_scans, example_ethogram, example_zones,
#'                      animals = example_animals)
#' res$diversity
#' res$spi
#' @export
analyse_group <- function(obs, ethogram, zones = NULL, animals = NULL,
                          sampling = c("auto", "scan", "focal"), n_boot = 0, ...) {
  obs <- as_obs(obs)
  eth <- as_ethogram(ethogram)
  sampling <- .choose_sampling(obs, match.arg(sampling))
  out <- list(sampling = sampling)
  out$budget <- behaviour_budget(obs, eth, sampling = sampling)
  out$diversity <- behaviour_diversity(obs, eth, sampling = sampling, n_boot = n_boot)
  out$diversity_group <- behaviour_diversity(obs, eth, by = NULL, sampling = sampling, n_boot = n_boot)
  if (!is.null(zones)) {
    out$spi <- spi(obs, zones, sampling = sampling)
    out$spi_group <- spi(obs, zones, by = NULL, sampling = sampling)
    out$zone_use <- zone_use(obs, zones, by = NULL, sampling = sampling)
  }
  ids <- if (!is.null(animals)) animals$animal else NULL
  am <- association_matrix(obs, sampling = sampling, animals = ids, ...)
  out$association <- network_metrics(am, attributes = animals)
  if (any(eth$directed)) {
    im <- tryCatch(interaction_matrix(obs, eth, sampling = sampling, animals = ids),
                   error = function(e) NULL)
    if (!is.null(im) && sum(im) > 0) out$interaction <- network_metrics(im, directed = TRUE, attributes = animals)
  }
  class(out) <- "ethogram_analysis"
  out
}

#' @export
print.ethogram_analysis <- function(x, ...) {
  cat("<ethogram_analysis> sampling:", x$sampling, "\n\n")
  cat("Species-typical behavioural diversity (per animal):\n")
  print(x$diversity[c("animal", "richness", "prop_repertoire", "shannon", "evenness", "prop_abnormal")],
        row.names = FALSE, digits = 3)
  if (!is.null(x$spi)) {
    cat("\nSpread of participation index:\n")
    print(x$spi, row.names = FALSE, digits = 3)
  }
  cat("\nAssociation network:\n")
  print(x$association$network, row.names = FALSE, digits = 3)
  if (!is.null(x$interaction)) {
    cat("\nInteraction network (directed):\n")
    print(x$interaction$network, row.names = FALSE, digits = 3)
  }
  invisible(x)
}
