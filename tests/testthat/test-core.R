scan_df <- function(animal, zone, behaviour = "Rest", scan_no = seq_along(animal), session = "s1") {
  as_obs(data.frame(session = session, sampling = "scan", record_type = "scan",
                    scan_no = scan_no, animal = animal, zone = zone,
                    behaviour = behaviour, stringsAsFactors = FALSE))
}

test_that("SPI matches hand calculation with unequal zones (Plowman 2003)", {
  z <- data.frame(zone = c("A", "B", "C"), area = c(1, 1, 2))
  # all in one zone -> 1
  expect_equal(spi(scan_df(rep("X", 4), rep("A", 4)), z)$spi, 1)
  # proportional to area -> 0
  expect_equal(spi(scan_df(rep("X", 4), c("A", "B", "C", "C")), z)$spi, 0)
  # fo = 2,2,0 ; fe = 1,1,2 ; sum|fo-fe| = 4 ; 2(N - fe_min) = 6
  expect_equal(spi(scan_df(rep("X", 4), c("A", "A", "B", "B")), z)$spi, 4 / 6)
})

test_that("SPI with equal zones is Dickens' index", {
  z <- data.frame(zone = c("A", "B", "C"))
  r <- spi(scan_df(rep("X", 4), c("A", "A", "B", "C")), z)
  fe <- 4 / 3
  expect_equal(r$spi, (abs(2 - fe) + 2 * abs(1 - fe)) / (2 * (4 - fe)))
  expect_equal(r$zones_used, 3)
})

test_that("zone electivity is 0 for proportional use", {
  z <- data.frame(zone = c("A", "B"), area = c(1, 3))
  u <- zone_use(scan_df(rep("X", 4), c("A", "B", "B", "B")), z)
  expect_equal(u$ivlev, c(0, 0))
  expect_equal(u$jacobs, c(0, 0))
})

test_that("diversity indices are correct", {
  eth <- data.frame(behaviour = c("A", "B", "C", "Pace"), class = c(NA, NA, NA, "abnormal"))
  o <- scan_df(rep("X", 5), rep("z", 5), behaviour = c("A", "A", "B", "B", "Pace"))
  d <- behaviour_diversity(o, eth)
  expect_equal(d$richness, 2)
  expect_equal(d$shannon, log(2))
  expect_equal(d$evenness, 1)
  expect_equal(d$simpson, 0.5)
  expect_equal(d$prop_repertoire, 2 / 3)
  expect_equal(d$evenness_repertoire, log(2) / log(3))
  expect_equal(d$prop_abnormal, 1 / 5)
  d2 <- behaviour_diversity(o, eth, species_typical_only = FALSE)
  expect_equal(d2$richness, 3)
})

test_that("OOS is excluded and reference similarity works", {
  eth <- data.frame(behaviour = c("A", "B"), reference = c(0.5, 0.5))
  o <- scan_df(rep("X", 3), c("z", "z", NA), behaviour = c("A", "B", "OOS"))
  d <- behaviour_diversity(o, eth)
  expect_equal(d$N, 2)
  expect_equal(d$similarity_ref, 1)
})

test_that("bootstrap CIs are returned and bracket the estimate", {
  d <- behaviour_diversity(example_scans, example_ethogram, by = NULL, n_boot = 100, seed = 1)
  expect_true(d$shannon_lo <= d$shannon && d$shannon <= d$shannon_hi)
})

test_that("SRI and HWI match hand calculation", {
  o <- scan_df(
    animal = c("P", "Q", "R", "P", "Q", "R", "P", "Q", "R", "P", "Q", "R"),
    zone = c("z1", "z1", "z2", "z1", "z1", NA, "z1", "z2", "z2", NA, "z1", "z1"),
    behaviour = c("Rest", "Rest", "Rest", "Rest", "Rest", "OOS", "Rest", "Rest", "Rest", "OOS", "Rest", "Rest"),
    scan_no = rep(1:4, each = 3))
  sri <- association_matrix(o, index = "SRI")
  hwi <- association_matrix(o, index = "HWI")
  expect_equal(unclass(sri)["P", "Q"], 0.5)
  expect_equal(unclass(sri)["Q", "R"], 0.5)
  expect_equal(unclass(sri)["P", "R"], 0)
  expect_equal(unclass(hwi)["P", "Q"], 2 / 3.5)
  expect_true(isSymmetric(matrix(as.numeric(sri), 3)))
})

test_that("focal_to_scans samples states at fixed intervals", {
  f <- as_obs(data.frame(session = "f", animal = "X", record_type = "state",
                         behaviour = c("Rest", "Feed"), start_s = c(0, 25), end_s = c(25, 40),
                         zone = c("a", "b")))
  s <- focal_to_scans(f, step = 10)
  expect_equal(s$behaviour, c("Rest", "Rest", "Rest", "Feed"))
  expect_equal(s$scan_no, 1:4)
})

test_that("interaction matrix counts bouts once and sums duration", {
  f <- as_obs(data.frame(session = "f", animal = c("X", "X", "Y"), record_type = "state",
                         behaviour = c("Groom", "Groom", "Groom"), recipients = c("Y", "Y", "X; Z"),
                         bout_id = c("b1", "b1", "b2"), start_s = c(0, 30, 0), end_s = c(30, 50, 10)))
  m <- interaction_matrix(f, behaviours = "Groom")
  expect_equal(unclass(m)["X", "Y"], 1)
  expect_equal(unclass(m)["Y", "Z"], 1)
  md <- interaction_matrix(f, behaviours = "Groom", measure = "duration")
  expect_equal(unclass(md)["X", "Y"], 50)
  expect_equal(unclass(md)["Y", "X"], 10)
})

test_that("network metrics run on example data", {
  net <- network_metrics(association_matrix(example_scans), attributes = example_animals)
  expect_s3_class(net, "ethogram_network")
  expect_equal(nrow(net$nodes), 6)
  expect_true("sex" %in% names(net$nodes))
  im <- interaction_matrix(example_focal, example_ethogram)
  dn <- network_metrics(im)
  expect_true(dn$network$directed)
  expect_equal(sum(dn$nodes$out_strength), sum(im))
})

test_that("CSV export reads and mixed sampling needs a choice", {
  f <- system.file("extdata", "ethogram_export_example.csv", package = "ethogramR")
  o <- read_ethogram_csv(f)
  expect_error(spi(o, example_zones), "both scan and focal")
  expect_equal(nrow(spi(o, example_zones, sampling = "scan")), 6)
})

test_that("functions with a seed leave the user's RNG state unchanged", {
  set.seed(42); before <- .Random.seed
  behaviour_diversity(example_scans, example_ethogram, by = NULL, n_boot = 20, seed = 1)
  network_layout(association_matrix(example_scans), seed = 3)
  expect_identical(.Random.seed, before)
})

test_that("matrix_to_edges and network_layout return plotting tables", {
  am <- association_matrix(example_scans)
  e <- matrix_to_edges(am)
  expect_equal(nrow(e), 15)                     # 6 animals, undirected
  expect_equal(nrow(matrix_to_edges(am, full = TRUE)), 36)
  lay <- network_layout(am)
  expect_true(all(c("x", "y", "xend", "yend") %in% names(lay$edges)))
  expect_true(all(c("x", "y", "strength") %in% names(lay$nodes)))
})

test_that("exact rarefied richness matches Hurlbert's formula", {
  # counts 3 A, 1 B ; n = 2: E[S] = (1 - C(1,2)/C(4,2)) + (1 - C(3,2)/C(4,2)) = 1 + (1 - 3/6) = 1.5
  expect_equal(ethogramR:::.rare_richness(c(A = 3, B = 1), 2), 1.5)
  expect_equal(ethogramR:::.rare_richness(c(A = 3, B = 1), 4), 2)
  expect_true(is.na(ethogramR:::.rare_richness(c(A = 3, B = 1), 5)))
})

test_that("rarefying to the full sample size returns the observed values", {
  d <- behaviour_diversity(example_scans, example_ethogram, rarefy = "min", n_rare = 5, seed = 1)
  full <- d$N == d$rarefy_n
  expect_true(any(full))
  expect_equal(d$shannon_rare[full], d$shannon[full])
  expect_equal(d$richness_rare[full], d$richness[full])
  expect_equal(d$shannon_rare_sd[full], 0)
  expect_true(all(d$rarefy_n == min(d$N)))
})

test_that("rarefy_within takes the minimum within each animal", {
  o <- example_scans
  o$period <- ifelse(o$session %in% c("Scan 01", "Scan 02"), "Pre", "Post")
  d <- behaviour_diversity(o, example_ethogram, by = c("animal", "period"),
                           rarefy = "min", rarefy_within = "animal", n_rare = 5, seed = 1)
  mins <- tapply(d$N, d$animal, min)
  expect_equal(d$rarefy_n, as.numeric(mins[d$animal]))
})

test_that("rarefaction is reproducible, order-independent and restores the RNG", {
  set.seed(9); before <- .Random.seed
  a <- behaviour_diversity(example_scans, example_ethogram, rarefy = 50, n_rare = 20, seed = 3)
  expect_identical(.Random.seed, before)
  b <- behaviour_diversity(example_scans[example_scans$animal == "F1", ], example_ethogram,
                           rarefy = 50, n_rare = 20, seed = 3)
  expect_equal(a$shannon_rare[a$animal == "F1"], b$shannon_rare)
  expect_true(all(a$richness_rare <= a$richness + 1e-12))
})

test_that("durations cannot be rarefied and too-small groups give NA", {
  expect_error(behaviour_diversity(example_focal, example_ethogram, rarefy = 10), "counted records")
  expect_warning(d <- behaviour_diversity(example_scans, example_ethogram, rarefy = 1e6, n_rare = 2),
                 "fewer records than the rarefaction size")
  expect_true(all(is.na(d$shannon_rare)))
  f <- behaviour_diversity(example_focal, example_ethogram, measure = "frequency",
                           rarefy = "min", n_rare = 5, seed = 1)
  expect_true(all(!is.na(f$shannon_rare)))
})

test_that("rarefaction_curve rises to the observed richness", {
  rc <- rarefaction_curve(example_scans, example_ethogram, steps = 8, n_rare = 10, seed = 1)
  last <- rc[rc$observed, ]
  d <- behaviour_diversity(example_scans, example_ethogram)
  expect_equal(last$richness[order(last$animal)], d$richness[order(d$animal)])
  expect_true(all(tapply(rc$richness, rc$animal, function(r) all(diff(r) >= -1e-9))))
})
