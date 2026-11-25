# Simulate example data in the Ethogram app CSV layout.
# Run from the package root: Rscript data-raw/make_example_data.R
set.seed(2026)

animals <- data.frame(
  animal = c("M1", "M2", "F1", "F2", "F3", "J1"),
  sex = c("M", "M", "F", "F", "F", "M"),
  age_class = c("adult", "adult", "adult", "adult", "subadult", "juvenile"),
  stringsAsFactors = FALSE
)
ids <- animals$animal

ethogram <- data.frame(
  behaviour = c("Rest", "Locomote", "Forage", "Feed", "Vigilant", "Allogroom",
                "Play", "Enrichment use", "Pacing", "Vocalise", "Displace"),
  code = c("RE", "LO", "FO", "FE", "VI", "AG", "PL", "EN", "PA", "VO", "DI"),
  type = c(rep("state", 9), "event", "event"),
  class = c("maintenance", "locomotion", "foraging", "foraging", "vigilance",
            "social", "social", "enrichment", "abnormal", "communication", "agonistic"),
  directed = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE),
  species_typical = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
  reference = c(0.35, 0.15, 0.25, 0.08, 0.07, 0.07, 0.03, NA, 0, NA, NA),
  stringsAsFactors = FALSE
)

# Zones as rectangles on the map (x0, x1, y0, y1 as proportions)
zones <- data.frame(
  enclosure = "Main",
  zone = c("Pool", "Rocks", "Grass", "Shelter", "Feeding", "Back"),
  area = c(150, 300, 800, 100, 200, 450),
  x0 = c(0, 0.25, 0.25, 0, 0.75, 0.75), x1 = c(0.25, 0.75, 0.75, 0.25, 1, 1),
  y0 = c(0, 0, 0.4, 0.4, 0, 0.4), y1 = c(0.4, 0.4, 1, 1, 0.4, 1),
  stringsAsFactors = FALSE
)
zone_pref <- c(Pool = 2, Rocks = 2.5, Grass = 3, Shelter = 2.5, Feeding = 1.5, Back = 0.6)

states <- ethogram$behaviour[ethogram$type == "state"]
base_p <- c(Rest = .33, Locomote = .14, Forage = .2, Feed = .08, Vigilant = .08,
            Allogroom = .08, Play = .03, `Enrichment use` = .04, Pacing = .02)
beh_p <- lapply(ids, function(a) {
  p <- base_p
  if (a == "F3") p["Pacing"] <- .15
  if (a == "J1") p["Play"] <- .1
  if (a == "M2") p["Allogroom"] <- .03
  p / sum(p)
})
names(beh_p) <- ids
subgroup <- c(M1 = "A", F1 = "A", F2 = "A", M2 = "B", F3 = "B", J1 = "B")
partner_w <- function(a) {
  w <- ifelse(subgroup[ids] == subgroup[a], 4, 1)
  names(w) <- ids
  w[a] <- 0
  w / sum(w)
}
xy_in_zone <- function(z) {
  r <- zones[zones$zone == z, ]
  c(round(stats::runif(1, r$x0, r$x1), 3), round(stats::runif(1, r$y0, r$y1), 3))
}
code_of <- function(b) ethogram$code[match(b, ethogram$behaviour)]
class_of <- function(b) ethogram$type[match(b, ethogram$behaviour)]

cols <- c("project", "session", "observer", "session_datetime", "sampling", "record_type",
          "scan_no", "animal", "co_actors", "recipients", "role", "behaviour",
          "behaviour_code", "behaviour_class", "zone", "enrichment_item", "x", "y",
          "start_s", "end_s", "duration_s", "count", "start_clock", "end_clock", "note",
          "video_file", "bout_id")
blank_row <- function() as.list(stats::setNames(rep(NA, length(cols)), cols))

# ---- Scans ---------------------------------------------------------------
scan_rows <- list()
for (s in 1:4) {
  sess <- sprintf("Scan %02d", s)
  for (k in 1:24) {
    zA <- sample(zones$zone, 1, prob = zone_pref)
    zB <- sample(zones$zone, 1, prob = zone_pref)
    rec <- list()
    for (a in sample(ids)) {
      if (!is.null(rec[[a]])) next
      r <- blank_row()
      r$project <- "Example group"; r$session <- sess; r$observer <- "SG"
      r$session_datetime <- sprintf("2026-11-%02d 09:00", 10 + s)
      r$sampling <- "scan"; r$record_type <- "scan"; r$scan_no <- k; r$animal <- a
      r$role <- "actor"
      if (stats::runif(1) < 0.06) {
        r$behaviour <- "OOS"; r$behaviour_code <- "OOS"
        rec[[a]] <- r; next
      }
      b <- sample(states, 1, prob = beh_p[[a]])
      own <- if (subgroup[a] == "A") zA else zB
      z <- if (stats::runif(1) < 0.75) own else sample(zones$zone, 1, prob = zone_pref)
      if (b == "Pacing") z <- "Back"
      if (b == "Feed") z <- "Feeding"
      if (b == "Enrichment use") r$enrichment_item <- "Puzzle feeder"
      r$behaviour <- b; r$behaviour_code <- code_of(b); r$behaviour_class <- class_of(b)
      r$zone <- z; xy <- xy_in_zone(z); r$x <- xy[1]; r$y <- xy[2]
      if (b %in% c("Allogroom", "Play")) {
        p <- sample(ids, 1, prob = partner_w(a))
        r$recipients <- p
        if (is.null(rec[[p]])) {
          rp <- r; rp$animal <- p; rp$role <- "recipient"; rp$recipients <- a
          rp$x <- min(1, r$x + 0.01); rp$enrichment_item <- NA
          rec[[p]] <- rp
        }
      }
      rec[[a]] <- r
    }
    scan_rows <- c(scan_rows, rec)
  }
}
scans <- do.call(rbind.data.frame, c(lapply(scan_rows, as.data.frame, stringsAsFactors = FALSE),
                                     stringsAsFactors = FALSE))

# ---- Focal (all animals followed continuously) -----------------------------
focal_rows <- list()
bout <- 0
for (s in 1:2) {
  sess <- sprintf("Focal %02d", s)
  len <- 1800
  for (a in ids) {
    t <- 0
    # starts out of sight
    t1 <- round(stats::runif(1, 5, 40))
    r <- blank_row(); r$session <- sess; r$animal <- a; r$behaviour <- "OOS"
    r$behaviour_code <- "OOS"; r$start_s <- 0; r$end_s <- t1; bout <- bout + 1
    r$bout_id <- paste0("b", bout); focal_rows[[length(focal_rows) + 1]] <- r
    t <- t1
    z <- sample(zones$zone, 1, prob = zone_pref)
    while (t < len) {
      b <- if (stats::runif(1) < 0.04) "OOS" else sample(states, 1, prob = beh_p[[a]])
      dur <- max(5, round(stats::rexp(1, 1 / 80)))
      end <- min(len, t + dur)
      bout <- bout + 1
      if (b == "Pacing") z <- "Back"
      if (b == "Feed") z <- "Feeding"
      partner <- if (b %in% c("Allogroom", "Play")) sample(ids, 1, prob = partner_w(a)) else NA
      # possibly split the bout by a location change
      cuts <- if (b != "OOS" && end - t > 60 && stats::runif(1) < 0.4) round((t + end) / 2) else numeric(0)
      bounds <- c(t, cuts, end)
      for (i in seq_len(length(bounds) - 1)) {
        if (i > 1) z <- sample(zones$zone, 1, prob = zone_pref)
        r <- blank_row(); r$session <- sess; r$animal <- a; r$behaviour <- b
        r$behaviour_code <- if (b == "OOS") "OOS" else code_of(b)
        r$start_s <- bounds[i]; r$end_s <- bounds[i + 1]; r$bout_id <- paste0("b", bout)
        r$role <- "actor"; r$recipients <- partner
        if (b != "OOS") {
          r$zone <- z; xy <- xy_in_zone(z); r$x <- xy[1]; r$y <- xy[2]
          r$behaviour_class <- "state"
        }
        focal_rows[[length(focal_rows) + 1]] <- r
      }
      if (b != "OOS" && stats::runif(1) < 0.12) {
        e <- blank_row(); e$session <- sess; e$animal <- a; e$record_type <- "event"
        eb <- if (stats::runif(1) < 0.5) "Vocalise" else "Displace"
        e$behaviour <- eb; e$behaviour_code <- code_of(eb); e$behaviour_class <- "event"
        e$start_s <- round(stats::runif(1, t, end)); e$count <- 1; e$role <- "actor"
        if (eb == "Displace") e$recipients <- sample(ids[ids != a], 1)
        e$zone <- z
        focal_rows[[length(focal_rows) + 1]] <- e
      }
      t <- end
    }
  }
}
focal <- do.call(rbind.data.frame, c(lapply(focal_rows, as.data.frame, stringsAsFactors = FALSE),
                                     stringsAsFactors = FALSE))
focal$project <- "Example group"; focal$observer <- "SG"; focal$sampling <- "focal"
focal$session_datetime <- ifelse(focal$session == "Focal 01", "2026-11-16 15:00", "2026-11-17 15:00")
focal$record_type[is.na(focal$record_type)] <- "state"
focal$duration_s <- ifelse(focal$record_type == "state", focal$end_s - focal$start_s, NA)
focal <- focal[order(focal$session, focal$start_s, focal$animal), ]
scans <- scans[cols]; focal <- focal[cols]

export <- rbind(scans, focal)
write.csv(export, "inst/extdata/ethogram_export_example.csv", row.names = FALSE, na = "")

devtools_free_load <- function() {
  for (f in list.files("R", full.names = TRUE)) source(f)
}
devtools_free_load()
example_scans <- as_obs(scans)
example_focal <- as_obs(focal)
example_ethogram <- ethogram
example_zones <- zones  # x0-x1, y0-y1: zone rectangle on the map (proportions)
example_animals <- animals
save(example_scans, file = "data/example_scans.rda", compress = "xz")
save(example_focal, file = "data/example_focal.rda", compress = "xz")
save(example_ethogram, file = "data/example_ethogram.rda", compress = "xz")
save(example_zones, file = "data/example_zones.rda", compress = "xz")
save(example_animals, file = "data/example_animals.rda", compress = "xz")
cat("scans:", nrow(scans), " focal:", nrow(focal), "\n")

# Sample CSV files for users
dir.create("inst/extdata/sample_data", showWarnings = FALSE)
w <- function(x, f) write.csv(x, file.path("inst/extdata/sample_data", f), row.names = FALSE, na = "")
w(scans, "observations_scan.csv")
w(focal, "observations_focal.csv")
w(export, "observations_all_ethogram_export.csv")
w(ethogram, "ethogram.csv")
w(zones, "zones.csv")
w(animals, "animals.csv")
