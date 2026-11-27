---
title: 'ethogramR: behavioural diversity, enclosure use and social networks from scan and focal observations of captive animal groups'
tags:
  - R
  - animal behaviour
  - animal welfare
  - zoo biology
  - ethology
  - social network analysis
  - enclosure use
authors:
  - name: Sitendu Goswami
    orcid: 0000-0000-0000-0000
    corresponding: true
    affiliation: 1
affiliations:
  - name: Vantara, Jamnagar, Gujarat, India
    index: 1
date: 26 November 2026
bibliography: paper.bib
---

<!--
BEFORE SUBMITTING, check every comment like this one, then delete it.
1. Replace the ORCID above with your own (orcid.org).
2. Confirm the affiliation wording, and add co-authors if anyone else contributed.
3. Rewrite "Research impact statement" with real uses once the package has been
   applied to Vantara data (reports, theses, papers, other institutions).
4. Complete the AI usage disclosure with what YOU reviewed and changed.
5. JOSS (2026 rules) desk-rejects repositories with less than six months of public,
   spread-out development history. Submit no earlier than six months after the
   repository is made public.
-->

# Summary

Zoos, sanctuaries and rescue centres monitor the behaviour of the animals in
their care to judge welfare, test changes to husbandry and design better
enclosures. Most of this monitoring uses two classic methods: scan sampling,
where every animal's behaviour and location is recorded at fixed intervals, and
focal (continuous) sampling, where behaviour is logged as it happens
[@altmann1974]. Three summaries are especially common in welfare work: how
diverse an animal's species-typical behaviour is [@miller2020], how evenly it
uses its enclosure [@brereton2020], and who it associates and interacts with
[@rose2015; @farine2015].

`ethogramR` is an R package that turns raw scan and focal records into these
three summaries. The user supplies four small tables: the observations, the
ethogram (the list of behaviours, marking which are directed at other animals
and which are species-typical), the enclosure zones with their areas, and an
optional list of animals. The package returns activity budgets; behavioural
diversity with bootstrap confidence intervals; the spread of participation
index (SPI) and zone electivity; and association and interaction networks with
standard node- and network-level metrics (\autoref{fig:example}). It reads the
CSV export of the Ethogram observation app directly, and any other data after a
one-line column renaming.

# Statement of need

Each of these measures is simple on paper, but computing them correctly from
field records involves many small decisions that are easy to get wrong and are
rarely reported. Out-of-sight records must be excluded from denominators.
Behaviours that are not species-typical, such as stereotypies, must be handled
consistently, because their inclusion changes diversity scores
[@cronin2019]. The SPI must use expected frequencies weighted by zone area when
zones differ in size [@plowman2003], and zones that were never used must still
count. Association indices must account for animals that were missed in some
scans [@cairns1987; @hoppitt2018], and directed behaviours recorded for both
actor and recipient must not be counted twice.

Done by hand in spreadsheets, these calculations are slow, hard to check and
difficult to repeat consistently across months of monitoring or across
institutions. `ethogramR` is aimed at welfare scientists,
zoo researchers, keepers and students who collect behaviour data on captive
groups and want these standard measures computed reproducibly, with the
assumptions written down in code and documentation. Because every result is
returned as a tidy data frame, it fits directly into routine reporting with
`ggplot2` [@wickham2016].

# State of the field

Data-collection tools such as ZooMonitor [@wark2019], BORIS [@friard2016] and
Zoo-Observer [@hahn2025] record observations and offer built-in summaries such
as time budgets and heat maps. They are primarily data-collection tools,
however, and their summaries cannot easily be scripted, extended or combined
into a reproducible analysis pipeline. In R, general network
packages such as `igraph` [@csardi2006] and the animal-specific packages
`asnipe` [@farine2013] and `ANTs` [@sosa2020] provide excellent network
construction, metrics and permutation tests. However, they start from
group-by-individual or interaction matrices, not from scan and focal records,
and they do not calculate behavioural diversity or enclosure use. Diversity
indices are available in ecology packages, but those treat behaviours as
species and know nothing of ethograms, out-of-sight records or species-typical
repertoires. To our knowledge, no package on CRAN computes the spread of participation
index; the existing package named `spi` computes an unrelated precipitation
index.

`ethogramR` therefore fills the step between the observation sheet and these
established tools, rather than duplicating them. Its association and
interaction matrices are ordinary R matrices, so users who need permutation
tests or multilevel network models can pass them straight to `asnipe` or
`ANTs`, which are better suited to that task than a new implementation would
be [@farine2022].

# Software design

The package is organised around one data model: the column layout of the
Ethogram app's CSV export, in which each row is either a scan sample, a timed
state, an event, or an environmental note. `as_obs()` converts any input to this
layout, so every analysis function works the same way on app exports and on
data from other sources. Several design choices follow from the behaviour of
real observation data.

**Scan and focal data are never silently mixed.** A single export can contain
both, but proportions of scans and proportions of time are different
quantities, so functions ask the user to choose a sampling method.
`focal_to_scans()` converts continuous records into point samples when a
point-sample analysis such as an association index is wanted.

**The ethogram carries the biology.** Marking behaviours as directed or
species-typical in the ethogram, rather than in each function call, keeps
decisions consistent across analyses. Diversity is calculated on species-typical
behaviours by default, with time in other behaviours reported separately as
`prop_abnormal`, following the concerns raised by @cronin2019. Alongside
Shannon's $H'$ [@shannon1948], Pielou's evenness [@pielou1966] and Simpson's
index [@brereton2022diversity], the package reports the share of the
repertoire expressed and $H'$ scaled by the full repertoire, so that a group
showing fewer behaviours scores lower even when its budget is even. If the
ethogram includes a reference budget, for example from wild animals, the
similarity between observed and reference budgets is also returned.

**Enclosure use follows the published formulas exactly.** `spi()` uses Dickens'
formula for equal zones [@dickens1955] and Plowman's unequal-zone modification
when areas are given [@plowman2003]. `zone_use()` adds Ivlev's and Jacobs'
electivity [@jacobs1974] so users can see which zones drive the SPI, as
recommended when choosing enclosure-use indices [@brereton2022indices].

**Networks reuse established tools.** Association matrices use the simple ratio
or half-weight index [@cairns1987; @hoppitt2018] computed from shared zones or
map distances. Interaction matrices use only actor records, so recipient rows
that the app creates for directed behaviours are not counted twice. Network
metrics are computed with `igraph`, including weighted clustering
[@barrat2004] and Louvain communities [@blondel2008], rather than being
re-implemented. Stochastic steps take a `seed` argument and restore the user's
random number state afterwards.

**Outputs are tables, not objects.** Every result is a long-format data frame,
and two helpers, `matrix_to_edges()` and `network_layout()`, convert networks
into tables for plotting. This favours transparency and reuse in reports over a
custom plotting system.

Correctness is checked by unit tests that compare results with values
calculated by hand for small examples, including the SPI, association indices
and diversity indices. The package passes `R CMD check --as-cran`, a
continuous-integration workflow repeats the check on Linux, macOS and Windows,
and the documentation, vignettes and simulated example data are built into a
`pkgdown` website.

![Example output from the package's simulated data for a group of six animals.
(A) Species-typical behavioural diversity with 95% bootstrap intervals,
coloured by the share of scans in abnormal behaviour. (B) Zone electivity
(Jacobs' D) and the spread of participation index for each animal. (C)
Association network from scan co-occurrence (simple ratio index), coloured by
Louvain community.\label{fig:example}](figure1.png)

# Research impact statement

<!-- Replace this section with evidence of real use before submitting. -->

`ethogramR` was developed for behaviour and welfare monitoring of animal groups
at Vantara, Jamnagar, India. It is paired with the Ethogram app used to collect
observations, so data move from collection to analysis without manual
reformatting. The bundled simulated data show that the
package recovers known structure: the association network separates the two
subgroups built into the simulation, and the animal simulated with frequent
pacing has the highest share of abnormal behaviour. No animal data were
collected for this paper; all examples use simulated data.

# AI usage disclosure

<!-- Complete the last sentence with what you personally reviewed and changed. -->

Generative AI was used extensively. Claude Opus 5.5 (Anthropic), used through
claude.ai in November 2026, wrote most of the R code, unit tests,
documentation, simulated example data and `pkgdown` website, and produced the
first draft of this paper, working from the author's specification of the
measures, input tables and data format. Correctness was checked by unit tests
against hand-calculated values for the SPI, association indices and diversity
indices; by `R CMD check --as-cran` and continuous integration; and by
checking each formula and reference against the original publications.
[AUTHOR: describe here what you personally reviewed, tested on real data and
changed in the code and in this paper.]

# Acknowledgements

<!-- Add funding, colleagues and data providers here, or delete this section. -->

# References
