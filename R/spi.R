# Prepare located records (zone known) with enclosure assigned
.prep_zones <- function(obs, zones, sampling, include_recipient, enclosure) {
  z <- as_zones(zones)
  if (!is.null(enclosure)) {
    z <- z[z$enclosure %in% enclosure, , drop = FALSE]
    if (!nrow(z)) stop("No zones for enclosure(s): ", paste(enclosure, collapse = ", "), call. = FALSE)
  }
  d <- .prep_weights(obs, sampling, "duration", include_recipient)
  d <- d[!d$oos, , drop = FALSE]
  n_nozone <- sum(is.na(d$zone))
  if (n_nozone) {
    message(n_nozone, " in-sight record(s) without a zone were excluded.")
    d <- d[!is.na(d$zone), , drop = FALSE]
  }
  has_encl <- "enclosure" %in% names(d) && any(!is.na(d$enclosure))
  if (!has_encl) {
    n_encl <- tapply(z$enclosure, z$zone, function(e) length(unique(e)))
    amb <- intersect(names(n_encl)[n_encl > 1], unique(d$zone))
    if (length(amb)) {
      stop("Zone name(s) used in more than one enclosure: ", paste(amb, collapse = ", "),
           ". Add an 'enclosure' column to the observations or set `enclosure`.", call. = FALSE)
    }
    d$enclosure <- z$enclosure[match(d$zone, z$zone)]
  }
  key_d <- paste(d$enclosure, d$zone, sep = "\r")
  key_z <- paste(z$enclosure, z$zone, sep = "\r")
  unk <- !key_d %in% key_z
  if (any(unk)) {
    if (is.null(enclosure)) {
      warning("Zone(s) not in the zones table were excluded: ",
              paste(unique(d$zone[unk]), collapse = ", "), call. = FALSE)
    }
    d <- d[!unk, , drop = FALSE]
  }
  list(d = d, z = z)
}

.zone_table <- function(g, z) {
  zz <- z[z$enclosure == g$enclosure[1], , drop = FALSE]
  fo <- tapply(g$w, g$zone, sum)
  fo <- as.numeric(fo[zz$zone]); fo[is.na(fo)] <- 0
  wt <- if (all(is.na(zz$weight))) rep(1, nrow(zz)) else zz$weight
  N <- sum(fo)
  data.frame(zone = zz$zone, area = zz$area, observed = fo,
             expected = N * wt / sum(wt),
             prop_observed = if (N > 0) fo / N else NA_real_,
             prop_expected = wt / sum(wt), stringsAsFactors = FALSE)
}

#' Spread of participation index (enclosure use)
#'
#' Calculates the spread of participation index (SPI) for each animal (or the
#' pooled group) in each enclosure. With equal-sized zones this is Dickens'
#' (1955) index; with zone areas (or weights) it is the unequal-zone version
#' of Plowman (2003):
#' \deqn{SPI = \frac{\sum |f_o - f_e|}{2 (N - f_{e,min})}}
#' where \eqn{f_o} is the observed number of scans (or seconds, for focal data)
#' in a zone, \eqn{f_e} the number expected if use were proportional to zone
#' size, \eqn{N} the total and \eqn{f_{e,min}} the expected value for the
#' smallest zone. SPI runs from 0 (use spread evenly in proportion to size) to
#' 1 (only one zone used). Zones never used count as 0. Out-of-sight records
#' and records without a zone are excluded.
#'
#' @inheritParams behaviour_budget
#' @param zones Zones table ([as_zones()]): `zone`, optional `enclosure`,
#'   `area` (or `weight`). Without areas, zones are treated as equal.
#' @param enclosure Optional enclosure name(s) to analyse.
#' @param include_recipient For scans, include recipient rows (located at the
#'   actor's zone by the app).
#' @return A data frame: grouping columns, `enclosure`, `n_zones`,
#'   `zones_used`, `N`, `spi`. The per-zone table is in `attr(x, "zones")`;
#'   see also [zone_use()].
#' @references Dickens, M. (1955). A statistical formula to quantify the
#'   "spread-of-participation" in group discussion. Speech Monographs 22: 28-31.
#'
#'   Plowman, A. B. (2003). A note on a modification of the spread of
#'   participation index allowing for unequal zones. Applied Animal Behaviour
#'   Science 83: 331-336. \doi{10.1016/S0168-1591(03)00142-4}
#' @examples
#' spi(example_scans, example_zones)
#' spi(example_scans, example_zones, by = NULL)  # whole group
#' @export
spi <- function(obs, zones, by = "animal", sampling = c("auto", "scan", "focal"),
                enclosure = NULL, include_recipient = TRUE) {
  obs <- as_obs(obs)
  sampling <- .choose_sampling(obs, match.arg(sampling))
  p <- .prep_zones(obs, zones, sampling, include_recipient, enclosure)
  zt_all <- list()
  out <- .by_groups(p$d, c(by, "enclosure"), function(g) {
    zt <- .zone_table(g, p$z)
    N <- sum(zt$observed)
    fe_min <- min(zt$expected)
    val <- if (N > 0 && nrow(zt) > 1) sum(abs(zt$observed - zt$expected)) / (2 * (N - fe_min)) else NA_real_
    lab <- paste(c(unlist(g[1, by, drop = FALSE]), g$enclosure[1]), collapse = " | ")
    zt_all[[lab]] <<- zt
    data.frame(n_zones = nrow(zt), zones_used = sum(zt$observed > 0), N = N, spi = val)
  })
  out <- out[c(by, "enclosure", "n_zones", "zones_used", "N", "spi")]
  attr(out, "zones") <- zt_all
  attr(out, "units") <- if (sampling == "scan") "scans" else "seconds"
  out
}

#' Zone use and electivity
#'
#' Observed and expected use of each zone, with Ivlev's electivity
#' \eqn{E = (r - p)/(r + p)} and Jacobs' index
#' \eqn{D = (r - p)/(r + p - 2rp)}, where \eqn{r} is the observed and \eqn{p}
#' the expected (area-based) proportion. Both run from -1 (avoided) to +1
#' (preferred); 0 means use in proportion to size.
#'
#' @inheritParams spi
#' @return A data frame with one row per group and zone.
#' @references Ivlev, V. S. (1961). Experimental ecology of the feeding of
#'   fishes. Yale University Press.
#'
#'   Jacobs, J. (1974). Quantitative measurement of food selection. Oecologia
#'   14: 413-417. \doi{10.1007/BF00384581}
#' @examples
#' zone_use(example_scans, example_zones, by = NULL)
#' @export
zone_use <- function(obs, zones, by = "animal", sampling = c("auto", "scan", "focal"),
                     enclosure = NULL, include_recipient = TRUE) {
  obs <- as_obs(obs)
  sampling <- .choose_sampling(obs, match.arg(sampling))
  p <- .prep_zones(obs, zones, sampling, include_recipient, enclosure)
  out <- .by_groups(p$d, c(by, "enclosure"), function(g) {
    zt <- .zone_table(g, p$z)
    r <- zt$prop_observed; q <- zt$prop_expected
    zt$ivlev <- ifelse(r + q > 0, (r - q) / (r + q), NA_real_)
    den <- r + q - 2 * r * q
    zt$jacobs <- ifelse(den > 0, (r - q) / den, NA_real_)
    zt
  })
  attr(out, "units") <- if (sampling == "scan") "scans" else "seconds"
  out
}
