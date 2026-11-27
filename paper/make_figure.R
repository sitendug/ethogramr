# Builds paper/figure1.png from the package's simulated example data.
# Run from the package root: Rscript paper/make_figure.R
library(ethogramR)
library(ggplot2)
library(patchwork)

theme_set(theme_minimal(base_size = 11) +
            theme(plot.title = element_text(face = "bold", size = 11),
                  panel.grid.minor = element_blank()))

div <- behaviour_diversity(example_scans, example_ethogram, n_boot = 1000, seed = 1)
div <- merge(div, example_animals)
pA <- ggplot(div, aes(x = reorder(animal, shannon), y = shannon, colour = prop_abnormal)) +
  geom_pointrange(aes(ymin = shannon_lo, ymax = shannon_hi), linewidth = 0.8) +
  coord_flip() +
  scale_colour_gradient(low = "#2a78b5", high = "#c44e52", labels = scales::percent,
                        name = "Abnormal") +
  labs(x = NULL, y = "Shannon H' (95% bootstrap CI)",
       title = "A  Species-typical behavioural diversity")

zu <- zone_use(example_scans, example_zones)
spi_tab <- spi(example_scans, example_zones)
zu$animal <- factor(zu$animal, levels = spi_tab$animal[order(spi_tab$spi)])
pB <- ggplot(zu, aes(x = zone, y = animal, fill = jacobs)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "#c44e52", mid = "white", high = "#2a78b5",
                       midpoint = 0, limits = c(-1, 1), name = "Jacobs' D") +
  geom_text(data = spi_tab, aes(x = 6.75, y = animal, label = sprintf("%.2f", spi)),
            inherit.aes = FALSE, size = 3, hjust = 0) +
  annotate("text", x = 6.75, y = 6.65, label = "SPI", size = 3, hjust = 0, fontface = "bold") +
  coord_cartesian(xlim = c(0.5, 7.2), clip = "off") +
  labs(x = NULL, y = NULL, title = "B  Zone electivity and SPI") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

net <- network_metrics(association_matrix(example_scans, index = "SRI"),
                       attributes = example_animals)
lay <- network_layout(net, seed = 2)
pC <- ggplot() +
  geom_segment(data = lay$edges, aes(x = x, y = y, xend = xend, yend = yend,
                                     linewidth = weight), colour = "grey55", lineend = "round") +
  geom_point(data = lay$nodes, aes(x = x, y = y, size = strength,
                                   fill = factor(community)), shape = 21, colour = "white") +
  geom_text(data = lay$nodes, aes(x = x, y = y, label = animal), size = 3.2) +
  scale_linewidth(range = c(0.2, 3), name = "SRI") +
  scale_size(range = c(8, 14), guide = "none") +
  scale_fill_manual(values = c("#e8a33d", "#4ca36b", "#8172b2", "#2a78b5"), name = "Community") +
  guides(fill = guide_legend(override.aes = list(size = 6))) +
  coord_equal(clip = "off") +
  labs(title = "C  Association network") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11))

fig <- (pA | pB) / pC + plot_layout(heights = c(1, 1.15))
ggsave("paper/figure1.png", fig, width = 9, height = 7.5, dpi = 300, bg = "white")
