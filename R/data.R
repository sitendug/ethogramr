#' Example data: simulated group of six animals
#'
#' Simulated observations in the Ethogram app CSV layout, for a mixed-sex group
#' of six animals in one enclosure with six zones. The simulation builds in
#' two preferred-partner subgroups, one animal that paces, and preferences for
#' the shelter and pool.
#'
#' \describe{
#'   \item{example_scans}{Scan samples: 4 sessions x 24 scans.}
#'   \item{example_focal}{Continuous focal records: 2 sessions of 30 min, all
#'     animals followed at once, with states and events.}
#'   \item{example_ethogram}{Ethogram with type, class, directed,
#'     species_typical and a wild `reference` budget.}
#'   \item{example_zones}{Six zones with areas (m^2) and their rectangles on
#'     the enclosure map (`x0`, `x1`, `y0`, `y1`, proportions from the top-left).}
#'   \item{example_animals}{Animal IDs, sex and age class.}
#' }
#' The same data are available as CSV files in
#' `system.file("extdata", "sample_data", package = "ethogramR")`, and the
#' combined app export is in
#' `system.file("extdata", "ethogram_export_example.csv", package = "ethogramR")`.
#' @name example_data
#' @aliases example_scans example_focal example_ethogram example_zones example_animals
#' @keywords datasets
NULL
