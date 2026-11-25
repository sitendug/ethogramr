# ethogramR: load the sample CSVs, analyse, and plot with ggplot2
# Set the folder that holds the CSV files:
data_dir <- system.file("extdata", "sample_data", package = "ethogramR")
# data_dir <- "~/Downloads/ethogramR_sample_data"   # or your own folder

library(ethogramR)
library(ggplot2)

# ---- 1. Load ----------------------------------------------------------------
scans    <- read_ethogram_csv(file.path(data_dir, "observations_scan.csv"))
focal    <- read_ethogram_csv(file.path(data_dir, "observations_focal.csv"))
ethogram <- read.csv(file.path(data_dir, "ethogram.csv"))
zones    <- read.csv(file.path(data_dir, "zones.csv"))
animals  <- read.csv(file.path(data_dir, "animals.csv"))

# ---- 2. Activity budget: stacked bars ----------------------------------------
budget <- behaviour_budget(scans, ethogram)
ggplot(budget, aes(x = animal, y = prop, fill = behaviour)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Animal", y = "Share of in-sight scans", fill = "Behaviour",
       title = "Activity budget (scan sampling)") +
  theme_minimal()

# ---- 3. Behavioural diversity with bootstrap 95% CIs -------------------------
div <- behaviour_diversity(scans, ethogram, n_boot = 1000, seed = 1)
ggplot(div, aes(x = reorder(animal, shannon), y = shannon)) +
  geom_pointrange(aes(ymin = shannon_lo, ymax = shannon_hi)) +
  coord_flip() +
  labs(x = NULL, y = "Shannon diversity H' (species-typical behaviours)",
       title = "Behavioural diversity") +
  theme_minimal()

# Diversity vs abnormal behaviour
ggplot(merge(div, animals), aes(x = prop_abnormal, y = shannon, colour = sex, label = animal)) +
  geom_point(size = 3) + geom_text(vjust = -1, show.legend = FALSE) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Scans in abnormal behaviour", y = "Shannon H'") +
  theme_minimal()

# ---- 4. Enclosure use: SPI and zone electivity -------------------------------
spi_tab <- spi(scans, zones)
ggplot(spi_tab, aes(x = reorder(animal, spi), y = spi)) +
  geom_col(fill = "#2a78b5") +
  coord_flip(ylim = c(0, 1)) +
  labs(x = NULL, y = "Spread of participation index (0 = even, 1 = one zone)",
       title = "Enclosure use") +
  theme_minimal()

zu <- zone_use(scans, zones)
ggplot(zu, aes(x = zone, y = animal, fill = jacobs)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "#c44e52", mid = "white", high = "#2a78b5",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(x = "Zone", y = "Animal", fill = "Jacobs' D",
       title = "Zone preference (+) and avoidance (-)") +
  theme_minimal()

# ---- 5. Enclosure map with observed locations --------------------------------
# App x/y are proportions from the top-left corner, so reverse the y axis.
ggplot() +
  geom_rect(data = zones, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1),
            fill = "grey95", colour = "grey60") +
  geom_text(data = zones, aes(x = (x0 + x1) / 2, y = (y0 + y1) / 2, label = zone),
            colour = "grey50", size = 5) +
  geom_point(data = subset(scans, !is.na(x)), aes(x = x, y = y, colour = animal),
             alpha = 0.6) +
  scale_y_reverse() + coord_equal() +
  labs(x = NULL, y = NULL, colour = "Animal", title = "Scan locations") +
  theme_void()

# ---- 6. Social networks -------------------------------------------------------
am  <- association_matrix(scans, method = "zone", index = "SRI")
net <- network_metrics(am, attributes = animals)

# Association heatmap
edges_full <- matrix_to_edges(am, full = TRUE)
ggplot(edges_full, aes(x = from, y = to, fill = weight)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(weight, 2)), size = 3) +
  scale_fill_gradient(low = "white", high = "#2a78b5") +
  labs(x = NULL, y = NULL, fill = "SRI", title = "Association index") +
  theme_minimal()

# Network drawn with ggplot
lay <- network_layout(net)
ggplot() +
  geom_segment(data = lay$edges, aes(x = x, y = y, xend = xend, yend = yend,
                                     linewidth = weight), colour = "grey60") +
  geom_point(data = lay$nodes, aes(x = x, y = y, size = strength, colour = sex)) +
  geom_text(data = lay$nodes, aes(x = x, y = y, label = animal)) +
  scale_linewidth(range = c(0.2, 3)) +
  scale_size(range = c(8, 16)) +
  labs(title = "Association network (SRI)") +
  theme_void()

# Node metrics
ggplot(net$nodes, aes(x = reorder(animal, strength), y = strength, fill = factor(community))) +
  geom_col() + coord_flip() +
  labs(x = NULL, y = "Strength", fill = "Community") +
  theme_minimal()

# Directed interactions (grooming time, focal data) with arrows
groom <- interaction_matrix(focal, behaviours = "Allogroom", measure = "duration")
lay_g <- network_layout(network_metrics(groom, directed = TRUE, attributes = animals))
ggplot() +
  geom_curve(data = lay_g$edges, aes(x = x, y = y, xend = xend, yend = yend,
                                     linewidth = weight),
             curvature = 0.2, colour = "grey50",
             arrow = arrow(length = unit(3, "mm"), type = "closed")) +
  geom_point(data = lay_g$nodes, aes(x = x, y = y, colour = sex), size = 12) +
  geom_text(data = lay_g$nodes, aes(x = x, y = y, label = animal)) +
  scale_linewidth(range = c(0.2, 2), name = "Seconds") +
  labs(title = "Grooming network (actor -> recipient)") +
  theme_void()
