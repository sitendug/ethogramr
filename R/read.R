.obs_cols <- c(
  project = "character", session = "character", observer = "character",
  session_datetime = "character", sampling = "character", record_type = "character",
  scan_no = "numeric", animal = "character", co_actors = "character",
  recipients = "character", role = "character", behaviour = "character",
  behaviour_code = "character", behaviour_class = "character", zone = "character",
  enrichment_item = "character", x = "numeric", y = "numeric", start_s = "numeric",
  end_s = "numeric", duration_s = "numeric", count = "numeric",
  start_clock = "character", end_clock = "character", note = "character",
  video_file = "character", bout_id = "character"
)

#' Read an Ethogram app CSV export
#'
#' Reads one or more CSV files exported from the Ethogram app (or any file with
#' the same columns) and standardises them with [as_obs()].
#'
#' @param file Path(s) to CSV file(s). Several files are combined.
#' @param map Optional named character vector to rename columns in your own data
#'   to the package's names, e.g. `c(animal = "ID", behaviour = "Behav")`.
#' @param ... Passed to [utils::read.csv()].
#' @return A data frame of observations (see [as_obs()]).
#' @examples
#' f <- system.file("extdata", "ethogram_export_example.csv", package = "ethogramR")
#' obs <- read_ethogram_csv(f)
#' table(obs$sampling)
#' @export
read_ethogram_csv <- function(file, map = NULL, ...) {
  dfs <- lapply(file, function(f) {
    utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE,
                    na.strings = c("", "NA"), ...)
  })
  all_cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(d) {
    for (nm in setdiff(all_cols, names(d))) d[[nm]] <- NA
    d[all_cols]
  })
  as_obs(do.call(rbind, dfs), map = map)
}

#' Standardise observation records
#'
#' Converts a data frame of scan and/or focal records to the column layout used
#' throughout the package (the Ethogram app CSV layout). Only `session`,
#' `animal` and `behaviour` are required; missing columns are added as `NA`.
#'
#' Key columns:
#' \describe{
#'   \item{session}{Session ID.}
#'   \item{sampling, record_type}{`"scan"` rows are instantaneous samples
#'     (one row per animal per scan, numbered by `scan_no`). Focal/video data
#'     have `record_type` `"state"` (with `start_s`, `end_s` or `duration_s`) or
#'     `"event"` (with `count`). `"environment"` rows are ignored.}
#'   \item{animal, behaviour, behaviour_code}{Who did what. `"OOS"` / `"Out of
#'     sight"` marks out-of-sight records.}
#'   \item{role, recipients, co_actors}{For directed (social) behaviours.
#'     `role` is `"actor"` or `"recipient"`; ID lists are separated by `;`, `,`
#'     or `|`.}
#'   \item{zone, x, y}{Location: zone name and map position as proportions.}
#'   \item{bout_id}{Rows of one bout split by a location change.}
#'   \item{enclosure}{Optional; otherwise inferred from zone names.}
#' }
#' If `record_type` is missing it is set to `"scan"` where `scan_no` is given
#' and `"state"` otherwise. Missing `duration_s` is computed from
#' `end_s - start_s`; missing event `count` is set to 1.
#'
#' @param x A data frame.
#' @param map Optional named character vector, `new_name = "old_name"`.
#' @return A data frame.
#' @export
as_obs <- function(x, map = NULL) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!is.null(map)) {
    for (nm in names(map)) {
      if (!map[[nm]] %in% names(x)) stop("Column '", map[[nm]], "' not found.", call. = FALSE)
      names(x)[names(x) == map[[nm]]] <- nm
    }
  }
  req <- c("session", "animal", "behaviour")
  miss <- setdiff(req, names(x))
  if (length(miss)) stop("Missing required column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  for (nm in names(.obs_cols)) {
    if (!nm %in% names(x)) x[[nm]] <- NA
    if (.obs_cols[[nm]] == "numeric") {
      x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
    } else {
      v <- as.character(x[[nm]])
      v[!is.na(v) & !nzchar(trimws(v))] <- NA
      x[[nm]] <- trimws(v)
    }
  }
  x$sampling <- tolower(x$sampling)
  x$record_type <- tolower(x$record_type)
  x$role <- tolower(x$role)
  rt_na <- is.na(x$record_type)
  x$record_type[rt_na] <- ifelse(!is.na(x$scan_no[rt_na]) | x$sampling[rt_na] %in% "scan",
                                 "scan", "state")
  s_na <- is.na(x$sampling)
  x$sampling[s_na] <- ifelse(x$record_type[s_na] == "scan", "scan", "focal")
  x$role[is.na(x$role)] <- "actor"
  d <- is.na(x$duration_s) & !is.na(x$start_s) & !is.na(x$end_s)
  x$duration_s[d] <- x$end_s[d] - x$start_s[d]
  e <- is.na(x$end_s) & !is.na(x$start_s) & !is.na(x$duration_s)
  x$end_s[e] <- x$start_s[e] + x$duration_s[e]
  ev <- x$record_type == "event" & is.na(x$count)
  x$count[ev] <- 1
  x
}

#' Standardise an ethogram
#'
#' @param x Data frame with at least a `behaviour` (or `name`) column. Optional
#'   columns: `code`; `type` (`"state"` or `"event"`); `class` (e.g.
#'   "social", "abnormal"); `directed` (logical, social behaviours with
#'   recipients); `species_typical` (logical); `reference` (proportion of time
#'   in a wild or reference population, used by [behaviour_diversity()]).
#' @details If `species_typical` is missing, behaviours whose name or class
#'   contains "abnormal" or "stereotyp" are marked `FALSE`, all others `TRUE`.
#'   Out-of-sight entries are removed.
#' @return A data frame.
#' @export
as_ethogram <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!"behaviour" %in% names(x)) {
    alt <- intersect(c("name", "behavior", "Behaviour", "Behavior"), names(x))
    if (!length(alt)) stop("The ethogram needs a 'behaviour' column.", call. = FALSE)
    names(x)[names(x) == alt[1]] <- "behaviour"
  }
  x$behaviour <- trimws(as.character(x$behaviour))
  if (!"code" %in% names(x)) x$code <- NA_character_
  if (!"type" %in% names(x)) x$type <- "state"
  x$type <- tolower(as.character(x$type))
  x$type[is.na(x$type)] <- "state"
  if (!"class" %in% names(x)) x$class <- NA_character_
  if (!"directed" %in% names(x)) x$directed <- FALSE
  x$directed <- .as_logical(x$directed)
  x$directed[is.na(x$directed)] <- FALSE
  if (!"species_typical" %in% names(x)) {
    x$species_typical <- !grepl("abnorm|stereotyp", paste(x$class, x$behaviour),
                                ignore.case = TRUE)
  }
  x$species_typical <- .as_logical(x$species_typical)
  x$species_typical[is.na(x$species_typical)] <- TRUE
  if ("reference" %in% names(x)) x$reference <- suppressWarnings(as.numeric(x$reference))
  x <- x[!.is_oos(x$behaviour, x$code), , drop = FALSE]
  if (anyDuplicated(tolower(x$behaviour))) stop("Duplicate behaviour names in the ethogram.", call. = FALSE)
  rownames(x) <- NULL
  x
}

#' Standardise enclosure zones
#'
#' @param x Data frame with a `zone` column and, optionally, `enclosure` and
#'   `area` (any unit, e.g. m^2). An optional `weight` column overrides `area`
#'   as the expected share of use (e.g. resource-weighted zones).
#' @return A data frame with `enclosure`, `zone`, `area`, `weight`.
#' @export
as_zones <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!"zone" %in% names(x)) stop("The zones table needs a 'zone' column.", call. = FALSE)
  x$zone <- trimws(as.character(x$zone))
  if (!"enclosure" %in% names(x)) x$enclosure <- "enclosure"
  x$enclosure <- as.character(x$enclosure)
  if (!"area" %in% names(x)) x$area <- NA_real_
  x$area <- suppressWarnings(as.numeric(x$area))
  if (!"weight" %in% names(x)) x$weight <- x$area
  x$weight <- suppressWarnings(as.numeric(x$weight))
  if (any(is.na(x$weight)) && !all(is.na(x$weight))) {
    stop("Give an area (or weight) for every zone, or for none (equal-sized zones).", call. = FALSE)
  }
  if (any(x$weight <= 0, na.rm = TRUE)) stop("Zone areas/weights must be positive.", call. = FALSE)
  if (anyDuplicated(paste(x$enclosure, x$zone))) stop("Duplicate zone names within an enclosure.", call. = FALSE)
  x[c("enclosure", "zone", "area", "weight", setdiff(names(x), c("enclosure", "zone", "area", "weight")))]
}

#' Convert continuous focal records into instantaneous samples
#'
#' Samples each animal's state (behaviour, zone, x, y, recipients) every `step`
#' seconds through each focal/video session, producing scan-format rows. This
#' lets continuous data be used where point samples are needed, e.g. for
#' association indices. In the Ethogram app all animals are followed at once
#' in a focal session, so every animal has a full timeline.
#'
#' @param obs Observations ([as_obs()]).
#' @param step Sampling interval in seconds.
#' @param offset Time of the first sample (seconds from session start).
#' @return Scan-format observations with `sampling = "scan_from_focal"`.
#' @export
focal_to_scans <- function(obs, step = 10, offset = 0) {
  obs <- as_obs(obs)
  d <- obs[!.is_scan(obs) & obs$record_type == "state" &
             !is.na(obs$start_s) & !is.na(obs$end_s), , drop = FALSE]
  if (!nrow(d)) stop("No focal state records with start and end times.", call. = FALSE)
  k1 <- ceiling((d$start_s - offset) / step)
  k2 <- ceiling((d$end_s - offset) / step) - 1
  k1 <- pmax(k1, 0)
  n <- pmax(k2 - k1 + 1, 0)
  idx <- rep(seq_len(nrow(d)), n)
  k <- unlist(mapply(function(a, b) if (b >= a) a:b else integer(0), k1, k2,
                     SIMPLIFY = FALSE), use.names = FALSE)
  out <- d[idx, , drop = FALSE]
  out$scan_no <- k + 1
  out$start_s <- offset + k * step
  out$end_s <- out$start_s
  out$duration_s <- NA_real_
  out$sampling <- "scan_from_focal"
  out$record_type <- "scan"
  # keep one row per animal per sample (first if overlapping)
  out <- out[!duplicated(paste(out$session, out$scan_no, out$animal)), , drop = FALSE]
  rownames(out) <- NULL
  out
}
