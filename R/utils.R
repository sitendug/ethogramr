# Internal helpers ---------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

.oos_labels <- c("oos", "out of sight", "out-of-sight", "not visible")

# TRUE for out-of-sight (or missing) behaviour records
.is_oos <- function(behaviour, code = NULL) {
  b <- tolower(trimws(as.character(behaviour)))
  r <- is.na(behaviour) | b %in% .oos_labels
  if (!is.null(code)) r <- r | (!is.na(code) & toupper(trimws(code)) == "OOS")
  r
}

# Split "A; B, C" style ID lists into character vectors
.split_ids <- function(x) {
  lapply(as.character(x), function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return(character(0))
    v <- trimws(strsplit(s, "[;,|]")[[1]])
    v[nzchar(v)]
  })
}

.as_logical <- function(x) {
  if (is.logical(x)) return(x)
  v <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(v))
  out[v %in% c("true", "t", "yes", "y", "1")] <- TRUE
  out[v %in% c("false", "f", "no", "n", "0")] <- FALSE
  out
}

.is_scan <- function(obs) {
  (!is.na(obs$record_type) & obs$record_type == "scan") |
    (!is.na(obs$sampling) & obs$sampling %in% c("scan", "scan_from_focal"))
}

# Resolve the sampling method to analyse
.choose_sampling <- function(obs, sampling = c("auto", "scan", "focal")) {
  sampling <- match.arg(sampling)
  sc <- .is_scan(obs)
  has_scan <- any(sc)
  has_cont <- any(!sc & obs$record_type %in% c("state", "event"))
  if (sampling == "auto") {
    if (has_scan && has_cont) {
      stop("The data contain both scan and focal/video records. ",
           "Set sampling = \"scan\" or sampling = \"focal\".", call. = FALSE)
    }
    sampling <- if (has_scan) "scan" else "focal"
  }
  if (sampling == "scan" && !has_scan) stop("No scan records found.", call. = FALSE)
  if (sampling == "focal" && !has_cont) stop("No focal/video records found.", call. = FALSE)
  sampling
}

# Unique ID for a scan (sampling point)
.scan_unit <- function(d) paste(d$session, d$scan_no, sep = "#")

# Unique ID for a focal bout (rows split by a location change share bout_id)
.bout_key <- function(d) {
  b <- ifelse(is.na(d$bout_id), paste0("row", seq_len(nrow(d))), d$bout_id)
  paste(d$session, d$animal, d$behaviour, b, sep = "#")
}

# Apply fun to each group of df defined by the `by` columns and bind results
.by_groups <- function(df, by, fun, ...) {
  if (length(by) == 0) return(fun(df, ...))
  miss <- setdiff(by, names(df))
  if (length(miss)) stop("Grouping column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
  df <- df[stats::complete.cases(df[by]), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  key <- interaction(df[by], drop = TRUE, lex.order = TRUE, sep = "\r")
  parts <- split(df, key)
  out <- lapply(parts, function(d) {
    r <- fun(d, ...)
    if (is.null(r) || nrow(r) == 0) return(NULL)
    cbind(d[rep(1, nrow(r)), by, drop = FALSE], r, row.names = NULL)
  })
  out <- do.call(rbind, out)
  if (!is.null(out)) rownames(out) <- NULL
  out
}

# Match observation behaviours to ethogram rows (by name, then by code)
.match_ethogram <- function(obs, eth, warn = TRUE) {
  m <- match(tolower(trimws(obs$behaviour)), tolower(trimws(eth$behaviour)))
  if (!all(is.na(eth$code))) {
    m2 <- match(toupper(trimws(obs$behaviour_code)), toupper(trimws(eth$code)))
    m[is.na(m)] <- m2[is.na(m)]
  }
  unk <- unique(obs$behaviour[is.na(m) & !.is_oos(obs$behaviour, obs$behaviour_code)])
  if (warn && length(unk)) {
    warning("Behaviour(s) not in the ethogram, treated as not species-typical: ",
            paste(unk, collapse = ", "), call. = FALSE)
  }
  m
}

.undirected <- function(g) {
  f <- if (exists("as_undirected", asNamespace("igraph"))) {
    get("as_undirected", asNamespace("igraph"))
  } else {
    get("as.undirected", asNamespace("igraph"))
  }
  f(g, mode = "collapse", edge.attr.comb = list(weight = "sum", "ignore"))
}

# Evaluate expr with a fixed seed, then restore the user's RNG state
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = genv, inherits = FALSE) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, envir = genv)
    else if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv)
  }, add = TRUE)
  set.seed(seed)
  expr
}
