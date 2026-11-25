# ethogramR 0.1.2

* Functions that take a `seed` (`behaviour_diversity()`, `network_layout()`,
  `plot_network()`) now restore the user's random number state afterwards.
* DOIs added to references; DESCRIPTION follows CRAN conventions.
* GitHub Actions check on Windows, macOS and Linux.

# ethogramR 0.1.1

* New `matrix_to_edges()` and `network_layout()` return networks as data
  frames for ggplot2.
* Sample CSV files and `analysis_with_ggplot.R` in
  `system.file("extdata", "sample_data", package = "ethogramR")`.
* `example_zones` gains zone rectangles (`x0`, `x1`, `y0`, `y1`) for drawing
  the enclosure map.
* Scan budgets list only state behaviours.

# ethogramR 0.1.0

* First version: import of Ethogram app exports, activity budgets,
  species-typical behavioural diversity, spread of participation index, zone
  electivity, association and interaction networks.
