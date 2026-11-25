# Prepare weighted records for budgets/diversity
.prep_weights <- function(obs, sampling, measure, include_recipient) {
  if (sampling == "scan") {
    d <- obs[.is_scan(obs), , drop = FALSE]
    if (!include_recipient) d <- d[d$role != "recipient", , drop = FALSE]
    d$w <- 1
    d$unit <- .scan_unit(d)
  } else {
    d <- obs[!.is_scan(obs) & obs$record_type %in% c("state", "event"), , drop = FALSE]
    if (measure == "duration") {
      d <- d[d$record_type == "state", , drop = FALSE]
      bad <- is.na(d$duration_s)
      if (any(bad)) {
        warning(sum(bad), " focal state row(s) without a duration were dropped.", call. = FALSE)
        d <- d[!bad, , drop = FALSE]
      }
      d$w <- d$duration_s
      d$unit <- .bout_key(d)
    } else {
      st <- d[d$record_type == "state", , drop = FALSE]
      st$unit <- .bout_key(st)
      st <- st[!duplicated(st$unit), , drop = FALSE]
      st$w <- 1
      ev <- d[d$record_type == "event", , drop = FALSE]
      ev$unit <- paste0("ev", seq_len(nrow(ev)))
      ev$w <- ev$count
      d <- rbind(st, ev)
    }
  }
  d$oos <- .is_oos(d$behaviour, d$behaviour_code)
  d
}

#' Activity budgets
#'
#' Proportion of scans (scan sampling) or of observed time (focal sampling) in
#' each behaviour, per animal or other grouping. Out-of-sight records are
#' excluded from the denominator unless `include_oos = TRUE`.
#'
#' @param obs Observations ([as_obs()] or [read_ethogram_csv()]).
#' @param ethogram Optional ethogram ([as_ethogram()]); adds `species_typical`
#'   and lists behaviours never seen with zero values.
#' @param by Grouping columns, e.g. `"animal"`, `c("animal", "session")`, or
#'   `NULL` for the whole group pooled.
#' @param sampling `"auto"`, `"scan"` or `"focal"` (focal includes video coding).
#' @param measure For focal data: `"duration"` (time in states) or
#'   `"frequency"` (number of bouts plus event counts).
#' @param include_recipient For scans, count rows where the animal is the
#'   recipient of a directed behaviour (e.g. being groomed).
#' @param include_oos Keep out-of-sight as a category.
#' @return A data frame: grouping columns, `behaviour`, `value` (scans,
#'   seconds or counts), `prop`, `n_total`.
#' @examples
#' behaviour_budget(example_scans, example_ethogram)
#' @export
behaviour_budget <- function(obs, ethogram = NULL, by = "animal",
                             sampling = c("auto", "scan", "focal"),
                             measure = c("duration", "frequency"),
                             include_recipient = TRUE, include_oos = FALSE) {
  obs <- as_obs(obs)
  sampling <- .choose_sampling(obs, match.arg(sampling))
  measure <- match.arg(measure)
  d <- .prep_weights(obs, sampling, measure, include_recipient)
  if (!include_oos) d <- d[!d$oos, , drop = FALSE]
  d$behaviour[d$oos] <- "OOS"
  eth <- if (!is.null(ethogram)) as_ethogram(ethogram) else NULL
  all_beh <- if (!is.null(eth)) {
    b <- eth$behaviour
    if (sampling == "scan" || measure == "duration") b <- eth$behaviour[eth$type == "state"]
    b
  } else character(0)
  if (!is.null(eth)) {
    m <- .match_ethogram(d, eth, warn = FALSE)
    d$behaviour[!is.na(m)] <- eth$behaviour[m[!is.na(m)]]
  }
  out <- .by_groups(d, by, function(g) {
    v <- tapply(g$w, g$behaviour, sum)
    beh <- union(all_beh, names(v))
    val <- as.numeric(v[beh]); val[is.na(val)] <- 0
    data.frame(behaviour = beh, value = val, prop = val / sum(val),
               n_total = sum(val), stringsAsFactors = FALSE)
  })
  if (!is.null(eth) && !is.null(out)) {
    out$species_typical <- eth$species_typical[match(tolower(out$behaviour), tolower(eth$behaviour))]
  }
  attr(out, "sampling") <- sampling
  attr(out, "measure") <- if (sampling == "scan") "scans" else measure
  out
}

# Diversity statistics from a named weight vector
.div_stats <- function(w, n_rep) {
  w <- w[w > 0]
  N <- sum(w)
  p <- w / N
  S <- length(p)
  H <- if (S > 0) -sum(p * log(p)) else NA_real_
  c(N = N, richness = S,
    repertoire_size = n_rep,
    prop_repertoire = if (n_rep > 0) S / n_rep else NA_real_,
    shannon = H,
    evenness = if (S > 1) H / log(S) else NA_real_,
    evenness_repertoire = if (n_rep > 1) H / log(n_rep) else NA_real_,
    simpson = if (S > 0) 1 - sum(p^2) else NA_real_,
    inv_simpson = if (S > 0) 1 / sum(p^2) else NA_real_)
}

#' Species-typical behavioural diversity
#'
#' Calculates behavioural diversity from scan or focal data, by default using
#' only species-typical behaviours (abnormal/stereotypic behaviours and
#' out-of-sight records are excluded from the diversity calculation but
#' reported separately as `prop_abnormal`).
#'
#' Indices, with \eqn{p_i} the proportion of scans (or time) in behaviour
#' \eqn{i} and \eqn{S} the number of behaviours seen:
#' \describe{
#'   \item{richness}{\eqn{S}, behaviours expressed.}
#'   \item{prop_repertoire}{\eqn{S} divided by the number of species-typical
#'     behaviours in the ethogram (repertoire expressed). For scan data and
#'     focal durations only `type = "state"` behaviours count towards the
#'     repertoire.}
#'   \item{shannon}{Shannon-Wiener \eqn{H' = -\sum p_i \ln p_i}.}
#'   \item{evenness}{Pielou's \eqn{J' = H'/\ln S}.}
#'   \item{evenness_repertoire}{\eqn{H'/\ln R}, with \eqn{R} the ethogram
#'     repertoire size, so groups expressing fewer behaviours score lower.}
#'   \item{simpson, inv_simpson}{\eqn{1-\sum p_i^2} and \eqn{1/\sum p_i^2}.}
#'   \item{prop_abnormal}{Share of all in-sight scans/time spent in
#'     non-species-typical behaviours.}
#'   \item{similarity_ref, shannon_ref}{If a reference budget is supplied
#'     (ethogram column `reference` or argument `reference`): percentage
#'     similarity \eqn{\sum \min(p_i, q_i)} (Schoener's index, 1 = identical
#'     budgets; 1 - Bray-Curtis dissimilarity) and \eqn{H'} of the reference.}
#' }
#'
#' Confidence intervals (`n_boot > 0`) are percentile bootstrap intervals,
#' resampling scans (scan data) or bouts (focal data) with replacement within
#' each group. Successive samples are autocorrelated, so treat these as
#' approximate; resampling sessions is more conservative when there are many.
#'
#' @inheritParams behaviour_budget
#' @param ethogram Ethogram ([as_ethogram()]), required.
#' @param species_typical_only Use only species-typical behaviours.
#' @param reference Optional data frame with `behaviour` and `prop` (e.g. a
#'   wild activity budget). Overrides the ethogram's `reference` column.
#' @param n_boot Number of bootstrap resamples (0 = none).
#' @param boot_unit `"sample"` (scans or bouts) or `"session"`.
#' @param conf Confidence level.
#' @param seed Optional random seed for the bootstrap. The user's random
#'   number state is restored afterwards.
#' @return A data frame with one row per group.
#' @references Shannon, C. E. (1948). A mathematical theory of communication.
#'   Bell System Technical Journal 27: 379-423.
#'
#'   Pielou, E. C. (1966). The measurement of diversity in different types of
#'   biological collections. Journal of Theoretical Biology 13: 131-144.
#'
#'   Miller, L. J., Vicino, G. A., Sheftel, J. & Lauderdale, L. K. (2020).
#'   Behavioral diversity as a potential indicator of positive animal welfare.
#'   Animals 10: 1211. \doi{10.3390/ani10071211}
#' @examples
#' behaviour_diversity(example_scans, example_ethogram)
#' behaviour_diversity(example_focal, example_ethogram, by = NULL)
#' @export
behaviour_diversity <- function(obs, ethogram, by = "animal",
                                sampling = c("auto", "scan", "focal"),
                                measure = c("duration", "frequency"),
                                species_typical_only = TRUE,
                                include_recipient = TRUE,
                                reference = NULL,
                                n_boot = 0, boot_unit = c("sample", "session"),
                                conf = 0.95, seed = NULL) {
  obs <- as_obs(obs)
  eth <- as_ethogram(ethogram)
  sampling <- .choose_sampling(obs, match.arg(sampling))
  measure <- match.arg(measure)
  boot_unit <- match.arg(boot_unit)
  d <- .prep_weights(obs, sampling, measure, include_recipient)
  d <- d[!d$oos, , drop = FALSE]
  m <- .match_ethogram(d, eth)
  d$behaviour[!is.na(m)] <- eth$behaviour[m[!is.na(m)]]
  d$typical <- !is.na(m) & eth$species_typical[ifelse(is.na(m), 1, m)]
  if (boot_unit == "session") d$unit <- d$session

  rep_eth <- eth
  # point samples and durations can only capture states, not events
  if (sampling == "scan" || measure == "duration") rep_eth <- rep_eth[rep_eth$type == "state", ]
  if (species_typical_only) rep_eth <- rep_eth[rep_eth$species_typical, ]
  n_rep <- nrow(rep_eth)

  ref <- NULL
  if (!is.null(reference)) {
    ref <- stats::setNames(as.numeric(reference$prop), trimws(reference$behaviour))
  } else if ("reference" %in% names(eth) && any(!is.na(eth$reference))) {
    ref <- stats::setNames(eth$reference, eth$behaviour)
  }
  if (!is.null(ref)) {
    ref <- ref[!is.na(ref)]
    if (species_typical_only) {
      ref <- ref[tolower(names(ref)) %in% tolower(eth$behaviour[eth$species_typical])]
    }
    ref <- ref / sum(ref)
  }
  alpha <- (1 - conf) / 2

  fun <- function(g) {
    prop_abn <- sum(g$w[!g$typical]) / sum(g$w)
    gg <- if (species_typical_only) g[g$typical, , drop = FALSE] else g
    w <- tapply(gg$w, gg$behaviour, sum)
    if (length(w) == 0) w <- stats::setNames(numeric(0), character(0))
    st <- .div_stats(w, n_rep)
    res <- as.data.frame(as.list(st))
    res$prop_abnormal <- prop_abn
    if (!is.null(ref)) {
      p <- w / sum(w)
      beh <- union(tolower(names(ref)), tolower(names(p)))
      pv <- stats::setNames(rep(0, length(beh)), beh)
      qv <- pv
      pv[tolower(names(p))] <- p
      qv[tolower(names(ref))] <- ref
      res$similarity_ref <- sum(pmin(pv, qv))
      qq <- ref[ref > 0]
      res$shannon_ref <- -sum(qq * log(qq))
    }
    if (n_boot > 0 && nrow(gg) > 0) {
      tab <- tapply(gg$w, list(gg$unit, gg$behaviour), sum)
      tab[is.na(tab)] <- 0
      bs <- replicate(n_boot, {
        i <- sample.int(nrow(tab), nrow(tab), replace = TRUE)
        .div_stats(colSums(tab[i, , drop = FALSE]), n_rep)[c("richness", "shannon", "evenness", "simpson")]
      })
      q <- apply(bs, 1, stats::quantile, probs = c(alpha, 1 - alpha), na.rm = TRUE)
      res$richness_lo <- q[1, "richness"]; res$richness_hi <- q[2, "richness"]
      res$shannon_lo <- q[1, "shannon"]; res$shannon_hi <- q[2, "shannon"]
      res$evenness_lo <- q[1, "evenness"]; res$evenness_hi <- q[2, "evenness"]
      res$simpson_lo <- q[1, "simpson"]; res$simpson_hi <- q[2, "simpson"]
    }
    res
  }
  out <- .with_seed(seed, .by_groups(d, by, fun))
  attr(out, "sampling") <- sampling
  attr(out, "units") <- if (sampling == "scan") "scans" else if (measure == "duration") "seconds" else "counts"
  out
}

#' Event rates
#'
#' Counts of events (e.g. vocalisations, aggression) per hour of in-sight
#' observation time, from focal/video data.
#'
#' @inheritParams behaviour_budget
#' @return Data frame with grouping columns, `behaviour`, `count`,
#'   `hours_observed`, `rate_per_hour`.
#' @export
event_rates <- function(obs, by = "animal") {
  obs <- as_obs(obs)
  d <- obs[!.is_scan(obs), , drop = FALSE]
  st <- d[d$record_type == "state" & !.is_oos(d$behaviour, d$behaviour_code) &
            !is.na(d$duration_s), , drop = FALSE]
  ev <- d[d$record_type == "event", , drop = FALSE]
  if (!nrow(ev)) stop("No event records found.", call. = FALSE)
  if (length(by)) {
    key_st <- do.call(paste, c(st[by], sep = "\r"))
    hours <- tapply(st$duration_s, key_st, sum) / 3600
  } else {
    hours <- sum(st$duration_s) / 3600
  }
  .by_groups(ev, by, function(g) {
    h <- if (length(by)) hours[do.call(paste, c(g[1, by, drop = FALSE], sep = "\r"))] else hours
    h <- as.numeric(h); if (is.na(h)) h <- NA_real_
    v <- tapply(g$count, g$behaviour, sum)
    data.frame(behaviour = names(v), count = as.numeric(v), hours_observed = h,
               rate_per_hour = as.numeric(v) / h, stringsAsFactors = FALSE)
  })
}
