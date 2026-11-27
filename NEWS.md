# ethogramR 0.1.3

* `behaviour_diversity()` gains rarefaction: `rarefy` (a sample size or
  `"min"`), `rarefy_within` (e.g. compare each animal's periods at the smaller
  of its two sample sizes) and `n_rare`. New columns `rarefy_n`,
  `richness_rare` (exact expected richness, Hurlbert 1971), `shannon_rare` and
  `shannon_rare_sd` (mean and SD of H' over random subsamples).
* New `rarefaction_curve()` returns expected richness and H' across sample
  sizes, to check whether sampling captured the repertoire.
* Rarefaction is refused for durations, which are not counted records.

* pkgdown website with a hex logo, a grouped function reference and a new
  article, "Plotting results with ggplot2". GitHub Actions rebuilds and
  publishes it to GitHub Pages on every push.
* JOSS paper draft in `paper/`, a workflow that builds its PDF, and
  CONTRIBUTING.md.

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
