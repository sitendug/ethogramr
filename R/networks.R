#' Association indices from scan samples
#'
#' Builds a symmetric association matrix from scans. In each scan, two animals
#' are "associated" if they are in the same zone (`method = "zone"`) or within
#' `max_dist` of each other on the enclosure map (`method = "distance"`).
#' Focal/video data are first converted to point samples every `focal_step`
#' seconds with [focal_to_scans()].
#'
#' For a dyad A-B: \eqn{x} = scans with A and B associated; \eqn{y_{AB}} =
#' scans with both seen but apart; \eqn{y_A}, \eqn{y_B} = scans with only A or
#' only B seen.
#' \describe{
#'   \item{SRI}{Simple ratio index \eqn{x / (x + y_{AB} + y_A + y_B)}.}
#'   \item{HWI}{Half-weight index \eqn{x / (x + y_{AB} + (y_A + y_B)/2)},
#'     preferred when one animal is more often missed than the other.}
#' }
#' An animal is "seen" in a scan if its record is in sight (and, for
#' `"distance"`, has x/y coordinates; for `"zone"`, has a zone).
#'
#' @inheritParams behaviour_budget
#' @param method `"zone"` or `"distance"`.
#' @param index `"SRI"` or `"HWI"`.
#' @param max_dist Distance threshold for `method = "distance"`, in the units
#'   of `map_size`.
#' @param map_size Width and height of the enclosure map, e.g. `c(60, 40)`
#'   metres. App coordinates are proportions, so the default `c(1, 1)` gives
#'   distances as a proportion of map width/height.
#' @param animals Optional vector of all animal IDs (to include animals never
#'   seen).
#' @param focal_step Sampling interval (s) used to convert focal data.
#' @return A symmetric matrix of association indices (0-1; diagonal 0) with
#'   attributes `counts` (list of x, yAB, yA, yB matrices), `n_scans` and
#'   `index`.
#' @references Cairns, S. J. & Schwager, S. J. (1987). A comparison of
#'   association indices. Animal Behaviour 35: 1454-1469.
#'   \doi{10.1016/S0003-3472(87)80018-0}
#'
#'   Hoppitt, W. J. E. & Farine, D. R. (2018). Association indices for
#'   quantifying social relationships: how to deal with missing observations
#'   of individuals or groups. Animal Behaviour 136: 227-238.
#'   \doi{10.1016/j.anbehav.2017.08.029}
#' @examples
#' association_matrix(example_scans)
#' @export
association_matrix <- function(obs, method = c("zone", "distance"),
                               index = c("SRI", "HWI"), max_dist = NULL,
                               map_size = c(1, 1), animals = NULL,
                               sampling = c("auto", "scan", "focal"),
                               focal_step = 10, include_recipient = TRUE) {
  obs <- as_obs(obs)
  method <- match.arg(method)
  index <- match.arg(index)
  sampling <- .choose_sampling(obs, match.arg(sampling))
  if (sampling == "focal") obs <- focal_to_scans(obs, step = focal_step)
  d <- obs[.is_scan(obs), , drop = FALSE]
  if (!include_recipient) d <- d[d$role != "recipient", , drop = FALSE]
  ids <- sort(unique(c(animals, d$animal[!is.na(d$animal)])))
  d <- d[!.is_oos(d$behaviour, d$behaviour_code) & !is.na(d$animal), , drop = FALSE]
  if (method == "zone") {
    d <- d[!is.na(d$zone), , drop = FALSE]
  } else {
    if (is.null(max_dist)) stop("Set max_dist for method = \"distance\".", call. = FALSE)
    d <- d[!is.na(d$x) & !is.na(d$y), , drop = FALSE]
  }
  d$unit <- .scan_unit(d)
  d <- d[!duplicated(paste(d$unit, d$animal)), , drop = FALSE]
  units <- unique(d$unit)
  n <- length(ids)
  S <- matrix(0, length(units), n, dimnames = list(units, ids))
  S[cbind(match(d$unit, units), match(d$animal, ids))] <- 1
  X <- matrix(0, n, n, dimnames = list(ids, ids))
  if (method == "zone") {
    for (zn in unique(d$zone)) {
      Z <- matrix(0, length(units), n)
      dz <- d[d$zone == zn, , drop = FALSE]
      Z[cbind(match(dz$unit, units), match(dz$animal, ids))] <- 1
      X <- X + crossprod(Z)
    }
  } else {
    for (u in units) {
      du <- d[d$unit == u, , drop = FALSE]
      if (nrow(du) < 2) next
      dm <- as.matrix(stats::dist(cbind(du$x * map_size[1], du$y * map_size[2])))
      a <- (dm <= max_dist) * 1
      i <- match(du$animal, ids)
      X[i, i] <- X[i, i] + a
    }
  }
  B <- crossprod(S)
  dimnames(B) <- list(ids, ids)
  yab <- B - X
  seen <- diag(B)
  ya <- matrix(seen, n, n) - B       # A seen, B not
  yb <- t(ya)                        # B seen, A not
  den <- if (index == "SRI") X + yab + ya + yb else X + yab + 0.5 * (ya + yb)
  ai <- ifelse(den > 0, X / den, NA_real_)
  dimnames(ai) <- list(ids, ids)
  diag(ai) <- 0
  diag(X) <- 0
  attr(ai, "counts") <- list(x = X, yAB = yab, yA = ya, yB = yb)
  attr(ai, "n_scans") <- length(units)
  attr(ai, "index") <- index
  class(ai) <- c("ethogram_matrix", "matrix", "array")
  ai
}

#' Directed interaction matrix
#'
#' Counts (or sums the duration of) directed social behaviours from each actor
#' (rows) to each recipient (columns), using the `recipients` column. Only
#' actor rows are used, so recipient rows in scans are not double counted.
#'
#' @inheritParams behaviour_budget
#' @param ethogram Optional ethogram; if given and `behaviours` is `NULL`, all
#'   behaviours with `directed = TRUE` are used.
#' @param behaviours Behaviour names or codes to include (e.g. `"Allogroom"`).
#'   If `NULL` and no ethogram is given, every record with recipients is used.
#' @param measure `"count"` (scans; bouts and event counts for focal data) or
#'   `"duration"` (seconds, focal states only).
#' @param rate If `TRUE`, divide each actor's row by its observation effort:
#'   in-sight hours (focal) or in-sight scans (scan), giving interactions per
#'   hour or per scan.
#' @param animals Optional vector of all animal IDs.
#' @return A square matrix, actors in rows and recipients in columns.
#' @examples
#' interaction_matrix(example_focal, example_ethogram, behaviours = "Allogroom")
#' @export
interaction_matrix <- function(obs, ethogram = NULL, behaviours = NULL,
                               measure = c("count", "duration"),
                               sampling = c("auto", "scan", "focal"),
                               rate = FALSE, animals = NULL) {
  obs <- as_obs(obs)
  measure <- match.arg(measure)
  sampling <- .choose_sampling(obs, match.arg(sampling))
  d <- if (sampling == "scan") obs[.is_scan(obs), , drop = FALSE] else
    obs[!.is_scan(obs) & obs$record_type %in% c("state", "event"), , drop = FALSE]
  effort_d <- d[!.is_oos(d$behaviour, d$behaviour_code) & d$role != "recipient", , drop = FALSE]
  ids <- sort(unique(c(animals, d$animal, unlist(.split_ids(d$recipients)))))
  ids <- ids[!is.na(ids)]
  d <- d[d$role != "recipient", , drop = FALSE]
  if (is.null(behaviours) && !is.null(ethogram)) {
    eth <- as_ethogram(ethogram)
    behaviours <- c(eth$behaviour[eth$directed], eth$code[eth$directed])
    if (!length(behaviours)) stop("No behaviours with directed = TRUE in the ethogram.", call. = FALSE)
  }
  if (!is.null(behaviours)) {
    b <- tolower(behaviours[!is.na(behaviours)])
    d <- d[tolower(d$behaviour) %in% b | tolower(d$behaviour_code) %in% b, , drop = FALSE]
  }
  rec <- .split_ids(d$recipients)
  keep <- lengths(rec) > 0
  d <- d[keep, , drop = FALSE]; rec <- rec[keep]
  if (sampling == "scan") {
    if (measure == "duration") stop("measure = \"duration\" needs focal data.", call. = FALSE)
    d$w <- 1
  } else if (measure == "duration") {
    st <- d$record_type == "state"
    d <- d[st, , drop = FALSE]; rec <- rec[st]
    d$w <- d$duration_s
  } else {
    d$w <- ifelse(d$record_type == "event", d$count, 1)
    st <- d$record_type == "state"
    first <- !st | !duplicated(.bout_key(d))   # one count per bout
    d <- d[first, , drop = FALSE]; rec <- rec[first]
  }
  M <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  if (nrow(d)) {
    actor <- rep(d$animal, lengths(rec))
    w <- rep(d$w, lengths(rec))
    r <- unlist(rec)
    ok <- actor != r & !is.na(w)
    agg <- tapply(w[ok], list(factor(actor[ok], ids), factor(r[ok], ids)), sum)
    agg[is.na(agg)] <- 0
    M[] <- agg
  }
  if (rate) {
    eff <- if (sampling == "scan") {
      tapply(rep(1, nrow(effort_d)), effort_d$animal, sum)
    } else {
      e <- effort_d[effort_d$record_type == "state", , drop = FALSE]
      tapply(e$duration_s, e$animal, sum, na.rm = TRUE) / 3600
    }
    eff <- as.numeric(eff[ids])
    M <- sweep(M, 1, ifelse(is.na(eff) | eff == 0, NA, eff), "/")
    attr(M, "effort") <- stats::setNames(eff, ids)
  }
  class(M) <- c("ethogram_matrix", "matrix", "array")
  attr(M, "measure") <- if (rate) paste(measure, if (sampling == "scan") "per scan" else "per hour") else measure
  M
}

#' Make a matrix symmetric
#'
#' @param m Square matrix.
#' @param how `"sum"`, `"mean"`, `"max"` or `"min"` of `m[i, j]` and `m[j, i]`.
#' @return Symmetric matrix.
#' @export
symmetrise <- function(m, how = c("sum", "mean", "max", "min")) {
  how <- match.arg(how)
  t_m <- t(m)
  out <- switch(how, sum = m + t_m, mean = (m + t_m) / 2,
                max = pmax(m, t_m), min = pmin(m, t_m))
  dimnames(out) <- dimnames(m)
  out
}

#' Social network metrics
#'
#' Builds a weighted network from an association or interaction matrix and
#' returns node-level and network-level metrics.
#'
#' Node metrics: `degree` (number of partners) and `strength` (sum of edge
#' weights), plus `in_`/`out_` versions for directed networks;
#' `eigenvector` centrality, `betweenness` (edge weights converted to distances
#' as 1/weight) and weighted `clustering` (Barrat et al. 2004), and the
#' `community` from Louvain modularity optimisation. Eigenvector, clustering
#' and communities use the undirected version of a directed network (weights
#' summed).
#'
#' Network metrics: number of nodes and edges, `density`, `mean_strength`,
#' `differentiation` (coefficient of variation of the dyadic weights, a
#' measure of how unequal relationships are), `mean_clustering`, `modularity`
#' and `n_communities`.
#'
#' @param m Square matrix from [association_matrix()] or
#'   [interaction_matrix()].
#' @param directed `NULL` (detect from symmetry), `TRUE` or `FALSE`.
#' @param attributes Optional data frame of animal attributes (sex, age...)
#'   joined to the node table.
#' @param id_col Name of the ID column in `attributes`.
#' @param seed Random seed for the stochastic steps (Louvain communities and
#'   the eigenvector start vector),
#'   so results are reproducible; the user's random number state is restored
#'   afterwards. `NULL` uses the current random state.
#' @return An object of class `ethogram_network`: a list with `nodes`,
#'   `network`, `graph` (an igraph object) and `matrix`.
#' @examples
#' net <- network_metrics(association_matrix(example_scans), attributes = example_animals)
#' net$nodes
#' net$network
#' @export
network_metrics <- function(m, directed = NULL, attributes = NULL, id_col = "animal",
                            seed = 1) {
  m <- matrix(as.numeric(m), nrow(m), ncol(m), dimnames = dimnames(m))
  if (nrow(m) != ncol(m)) stop("The matrix must be square.", call. = FALSE)
  if (is.null(rownames(m))) rownames(m) <- colnames(m) <- paste0("A", seq_len(nrow(m)))
  if (anyNA(m)) {
    message(sum(is.na(m)), " NA value(s) (dyads never observed) set to 0.")
    m[is.na(m)] <- 0
  }
  diag(m) <- 0
  if (is.null(directed)) directed <- !isSymmetric(unname(m))
  g <- igraph::graph_from_adjacency_matrix(m, mode = if (directed) "directed" else "upper",
                                           weighted = TRUE, diag = FALSE)
  gu <- if (directed) .undirected(g) else g
  ew <- igraph::E(g)$weight
  ewu <- igraph::E(gu)$weight
  ids <- igraph::V(g)$name
  nodes <- data.frame(animal = ids, stringsAsFactors = FALSE)
  if (directed) {
    nodes$out_degree <- igraph::degree(g, mode = "out")
    nodes$in_degree <- igraph::degree(g, mode = "in")
    nodes$out_strength <- igraph::strength(g, mode = "out", weights = ew)
    nodes$in_strength <- igraph::strength(g, mode = "in", weights = ew)
  }
  nodes$degree <- igraph::degree(gu)
  nodes$strength <- igraph::strength(gu, weights = ewu)
  has_edges <- igraph::ecount(g) > 0
  nodes$eigenvector <- if (has_edges) .with_seed(seed, igraph::eigen_centrality(gu, weights = ewu)$vector) else NA_real_
  nodes$betweenness <- if (has_edges) igraph::betweenness(g, directed = directed, weights = 1 / ew) else 0
  cl <- if (has_edges) igraph::transitivity(gu, type = "barrat", weights = ewu, isolates = "NaN") else rep(NaN, length(ids))
  nodes$clustering <- ifelse(is.nan(cl), NA_real_, cl)
  comm <- if (has_edges) .with_seed(seed, igraph::cluster_louvain(gu, weights = ewu)) else NULL
  nodes$community <- if (!is.null(comm)) as.integer(igraph::membership(comm)) else NA_integer_
  rownames(nodes) <- NULL
  if (!is.null(attributes)) {
    attributes <- as.data.frame(attributes)
    if (!id_col %in% names(attributes)) stop("Column '", id_col, "' not in attributes.", call. = FALSE)
    names(attributes)[names(attributes) == id_col] <- "animal"
    nodes <- merge(attributes, nodes, by = "animal", all.y = TRUE, sort = FALSE)
    nodes <- nodes[match(ids, nodes$animal), ]
    rownames(nodes) <- NULL
  }
  vals <- if (directed) m[row(m) != col(m)] else m[upper.tri(m)]
  network <- data.frame(
    nodes = igraph::vcount(g), edges = igraph::ecount(g), directed = directed,
    density = igraph::edge_density(g),
    mean_strength = mean(nodes$strength),
    differentiation = if (mean(vals) > 0) stats::sd(vals) / mean(vals) else NA_real_,
    mean_clustering = mean(nodes$clustering, na.rm = TRUE),
    modularity = if (!is.null(comm)) igraph::modularity(comm) else NA_real_,
    n_communities = if (!is.null(comm)) length(comm) else NA_integer_
  )
  structure(list(nodes = nodes, network = network, graph = g, matrix = m),
            class = "ethogram_network")
}

#' @export
print.ethogram_matrix <- function(x, digits = 3, ...) {
  extra <- c(index = attr(x, "index"), measure = attr(x, "measure"),
             scans = attr(x, "n_scans"))
  if (length(extra)) cat(paste(names(extra), extra, sep = ": ", collapse = " | "), "\n")
  m <- matrix(as.numeric(x), nrow(x), ncol(x), dimnames = dimnames(x))
  print(round(m, digits), ...)
  invisible(x)
}

#' @export
print.ethogram_network <- function(x, ...) {
  cat("<ethogram_network>", if (x$network$directed) "directed" else "undirected",
      "|", x$network$nodes, "animals,", x$network$edges, "edges\n\n")
  print(x$network, row.names = FALSE, digits = 3)
  cat("\n")
  print(x$nodes, row.names = FALSE, digits = 3)
  invisible(x)
}

#' Plot a social network
#'
#' @param x An `ethogram_network` from [network_metrics()] or a square matrix.
#' @param colour_by Node column to colour by (e.g. `"sex"`, `"community"`).
#' @param size_by Node column for node size (default `"strength"`).
#' @param edge_scale Multiplier for edge widths.
#' @param min_weight Hide edges below this weight.
#' @param layout An igraph layout function.
#' @param seed Random seed so the layout is reproducible (the user's random
#'   number state is restored afterwards).
#' @param ... Passed to [igraph::plot.igraph()].
#' @return The network object, invisibly.
#' @examples
#' net <- network_metrics(association_matrix(example_scans), attributes = example_animals)
#' plot_network(net, colour_by = "sex")
#' @export
plot_network <- function(x, colour_by = "community", size_by = "strength",
                         edge_scale = 8, min_weight = 0,
                         layout = igraph::layout_with_fr, seed = 1, ...) {
  if (!inherits(x, "ethogram_network")) x <- network_metrics(x)
  g <- x$graph
  g <- igraph::delete_edges(g, igraph::E(g)[igraph::E(g)$weight < min_weight])
  nd <- x$nodes
  pal <- c("#2a78b5", "#e07b39", "#4ca36b", "#c44e52", "#8172b2", "#937860",
           "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd")
  cols <- "#2a78b5"
  if (!is.null(colour_by) && colour_by %in% names(nd)) {
    f <- factor(nd[[colour_by]])
    cols <- pal[(as.integer(f) - 1) %% length(pal) + 1]
  }
  sz <- 18
  if (!is.null(size_by) && size_by %in% names(nd)) {
    s <- nd[[size_by]]
    sz <- if (max(s, na.rm = TRUE) > 0) 10 + 20 * s / max(s, na.rm = TRUE) else 15
  }
  w <- igraph::E(g)$weight
  wd <- if (length(w) && max(w) > 0) 0.5 + edge_scale * w / max(w) else 1
  xy <- .with_seed(seed, layout(g))
  igraph::plot.igraph(g, layout = xy, vertex.color = cols, vertex.size = sz,
                      vertex.label.color = "black", vertex.frame.color = "white",
                      edge.width = wd, edge.color = grDevices::adjustcolor("grey30", 0.6),
                      edge.arrow.size = 0.4, edge.curved = if (x$network$directed) 0.2 else 0, ...)
  if (!is.null(colour_by) && colour_by %in% names(nd)) {
    f <- factor(nd[[colour_by]])
    graphics::legend("bottomleft", legend = levels(f), pt.bg = pal[seq_along(levels(f))],
                     pch = 21, pt.cex = 1.5, bty = "n", title = colour_by)
  }
  invisible(x)
}

#' Convert a network matrix to an edge list
#'
#' Returns a long data frame with one row per dyad, ready for ggplot2 (e.g. a
#' heatmap with `geom_tile()`) or for export.
#'
#' @param m Square matrix from [association_matrix()] or [interaction_matrix()].
#' @param directed `NULL` (detect from symmetry), `TRUE` or `FALSE`. For
#'   undirected matrices each dyad appears once.
#' @param drop_zero Drop dyads with weight 0 (and NA).
#' @param full If `TRUE`, return every ordered pair (both triangles and the
#'   diagonal), which suits a heatmap.
#' @return Data frame with `from`, `to`, `weight`.
#' @examples
#' head(matrix_to_edges(association_matrix(example_scans)))
#' @export
matrix_to_edges <- function(m, directed = NULL, drop_zero = TRUE, full = FALSE) {
  m <- matrix(as.numeric(m), nrow(m), ncol(m), dimnames = dimnames(m))
  if (is.null(directed)) directed <- !isSymmetric(unname(m))
  ids <- rownames(m)
  keep <- if (full) matrix(TRUE, nrow(m), ncol(m)) else if (directed) row(m) != col(m) else upper.tri(m)
  idx <- which(keep, arr.ind = TRUE)
  out <- data.frame(from = ids[idx[, 1]], to = ids[idx[, 2]], weight = m[idx],
                    stringsAsFactors = FALSE)
  if (drop_zero && !full) out <- out[!is.na(out$weight) & out$weight > 0, , drop = FALSE]
  out <- out[order(out$from, out$to), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Network layout as data frames for ggplot2
#'
#' Computes node positions and returns two plain data frames: `nodes` (node
#' metrics plus `x`, `y`) and `edges` (`from`, `to`, `weight` plus `x`, `y`,
#' `xend`, `yend`), so a network can be drawn with `geom_segment()` and
#' `geom_point()` without any network package.
#'
#' @param x An `ethogram_network` ([network_metrics()]) or a square matrix.
#' @param layout An igraph layout function.
#' @param seed Random seed so the layout is reproducible (the user's random
#'   number state is restored afterwards).
#' @return A list with data frames `nodes` and `edges`.
#' @examples
#' net <- network_metrics(association_matrix(example_scans), attributes = example_animals)
#' lay <- network_layout(net)
#' head(lay$edges)
#' @export
network_layout <- function(x, layout = igraph::layout_with_fr, seed = 1) {
  if (!inherits(x, "ethogram_network")) x <- network_metrics(x)
  xy <- .with_seed(seed, layout(x$graph))
  ids <- igraph::V(x$graph)$name
  nodes <- x$nodes
  pos <- data.frame(animal = ids, x = xy[, 1], y = xy[, 2], stringsAsFactors = FALSE)
  nodes <- merge(nodes, pos, by = "animal", sort = FALSE)
  edges <- matrix_to_edges(x$matrix, directed = x$network$directed)
  i <- match(edges$from, pos$animal); j <- match(edges$to, pos$animal)
  edges$x <- pos$x[i]; edges$y <- pos$y[i]
  edges$xend <- pos$x[j]; edges$yend <- pos$y[j]
  list(nodes = nodes, edges = edges)
}
