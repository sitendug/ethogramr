# ethogramR <img src="man/figures/logo.png" align="right" height="139" alt="ethogramR logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/YOUR-GITHUB-USERNAME/ethogramR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YOUR-GITHUB-USERNAME/ethogramR/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**ethogramR** analyses behaviour observations of captive animal groups. Give
it scan and/or focal (continuous or video-coded) data, the ethogram and the
enclosure zones, and it calculates:

| Question | Functions |
|---|---|
| How do animals spend their time? | `behaviour_budget()`, `event_rates()` |
| How diverse is their species-typical behaviour? | `behaviour_diversity()`: richness, share of repertoire expressed, Shannon H′, Pielou J′, Simpson, proportion of abnormal behaviour, similarity to a wild reference budget, bootstrap CIs, rarefaction to equal sample sizes; `rarefaction_curve()` |
| How evenly do they use the enclosure? | `spi()`: spread of participation index (Dickens 1955; Plowman 2003 for unequal zones); `zone_use()`: Ivlev and Jacobs electivity |
| Who associates and interacts with whom? | `association_matrix()` (SRI / HWI from shared zone or map distance), `interaction_matrix()` (directed counts, durations, rates), `network_metrics()`, `plot_network()` |
| All of the above | `analyse_group()` |

It reads the CSV export of the **Ethogram** observation app directly, and any
other data once columns are renamed with `map =`. All results are data frames
that go straight into **ggplot2**; `matrix_to_edges()` and `network_layout()`
turn networks into tables for `geom_tile()` and `geom_segment()`.

Website and documentation: <https://YOUR-GITHUB-USERNAME.github.io/ethogramR/>

## Installation

```r
# install.packages("remotes")
remotes::install_github("YOUR-GITHUB-USERNAME/ethogramR", build_vignettes = TRUE)
```

## Inputs

| Table | Required columns | Optional columns |
|---|---|---|
| Observations | `session`, `animal`, `behaviour` | `sampling`, `record_type` (scan/state/event), `scan_no`, `start_s`, `end_s`, `duration_s`, `count`, `zone`, `x`, `y`, `role`, `recipients`, `bout_id`, `enclosure` |
| Ethogram | `behaviour` | `code`, `type` (state/event), `class`, `directed`, `species_typical`, `reference` |
| Zones | `zone` | `area` (or `weight`), `enclosure`, `x0`, `x1`, `y0`, `y1` |
| Animals | `animal` | any attributes, e.g. `sex`, `age_class` |

Out-of-sight records (`OOS`) are excluded. Behaviours with
`species_typical = FALSE` (e.g. stereotypies) are left out of diversity and
reported as `prop_abnormal`.

## Example

```r
library(ethogramR)
library(ggplot2)

f <- system.file("extdata", "sample_data", package = "ethogramR")
scans    <- read_ethogram_csv(file.path(f, "observations_scan.csv"))
ethogram <- read.csv(file.path(f, "ethogram.csv"))
zones    <- read.csv(file.path(f, "zones.csv"))
animals  <- read.csv(file.path(f, "animals.csv"))

# Species-typical behavioural diversity with bootstrap CIs
div <- behaviour_diversity(scans, ethogram, n_boot = 1000, seed = 1)
ggplot(div, aes(reorder(animal, shannon), shannon)) +
  geom_pointrange(aes(ymin = shannon_lo, ymax = shannon_hi)) +
  coord_flip() + labs(x = NULL, y = "Shannon H'")

# Enclosure use
spi(scans, zones)
zone_use(scans, zones, by = NULL)

# Association network
net <- network_metrics(association_matrix(scans), attributes = animals)
net$nodes
plot_network(net, colour_by = "sex")
```

The full workflow is in the
[Getting started](https://YOUR-GITHUB-USERNAME.github.io/ethogramR/articles/getting-started.html)
article, and
[Plotting results with ggplot2](https://YOUR-GITHUB-USERNAME.github.io/ethogramR/articles/plotting-with-ggplot2.html)
shows ten example plots.

## Contributing

Bug reports, questions and pull requests are welcome; see
[CONTRIBUTING.md](CONTRIBUTING.md).

## Citation

```r
citation("ethogramR")
```

## Key references

- Dickens, M. (1955). A statistical formula to quantify the "spread-of-participation" in group discussion. *Speech Monographs* 22: 28–31.
- Plowman, A. B. (2003). A note on a modification of the spread of participation index allowing for unequal zones. *Applied Animal Behaviour Science* 83: 331–336. <https://doi.org/10.1016/S0168-1591(03)00142-4>
- Miller, L. J. et al. (2020). Behavioral diversity as a potential indicator of positive animal welfare. *Animals* 10: 1211. <https://doi.org/10.3390/ani10071211>
- Cairns, S. J. & Schwager, S. J. (1987). A comparison of association indices. *Animal Behaviour* 35: 1454–1469. <https://doi.org/10.1016/S0003-3472(87)80018-0>
- Hoppitt, W. J. E. & Farine, D. R. (2018). Association indices for quantifying social relationships. *Animal Behaviour* 136: 227–238. <https://doi.org/10.1016/j.anbehav.2017.08.029>

## Licence

MIT © Sitendu Goswami
