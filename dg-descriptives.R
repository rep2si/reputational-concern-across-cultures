# Decriptives plot: histograms and Ligon-Schechter style
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggrepel)
library(broom)

source("./rep2si-common/source_all_submodule.R")
sharepoint_path <- Sys.getenv("REP2SI_DIR")

set_theme(theme_cowplot())
cond_cols <- c(anonymous = "#0072B2", revealed = "#D55E00")

# Convenience: two versions of dataset, one with all site labels marked "total"
d <- readRDS(file.path(sharepoint_path, "rep2si-data", "global-data", "allocs.rds"))
d_tot <- bind_rows(d %>% mutate(site = "total"))

# helper to make coherent histograms
base_hist <- function(data) {
  ggplot(data, aes(x = percGiven, fill = Condition)) +
    geom_histogram(
      aes(y = after_stat(density) / 10),
      binwidth = 0.1,
      position = position_dodge()
    ) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = scales::percent) +
    scale_fill_manual(values = cond_cols) +
    scale_colour_manual(values = cond_cols) +
  labs(
    x = "Share of endowment given",
    y = "Proportion (within condition)"
  )
}

# Plot: Global and site-wise histograms ----------------------------------------
p_total <- d_tot %>%
  base_hist() +
  facet_wrap(~ site, labeller = as_labeller(get_site_labels())) +
  theme(legend.position = "none")

p_rest <- d %>%
  filter(site != "total") %>%
  base_hist() +
  facet_wrap(~ site, nrow = 2, ncol = 4, , labeller = as_labeller(get_site_labels()))

# combine, shared legend at top
p_hists <- p_total / p_rest +
  plot_layout(heights = c(1, 2), guides = "collect", axis_titles = "collect") &
  theme(legend.position = "top")

# Plot: Ligon-Schechter style --------------------------------------------------
summaries <- d %>%
  group_by(Condition, site) %>%
  summarise(tidy(t.test(percGiven)), .groups = "drop") %>%
  select(site, Condition, estimate, conf.low, conf.high) %>%
  pivot_wider(
    id_cols = site,
    names_from = Condition,
    values_from = c(estimate, conf.low, conf.high)
  )

p_ls <- summaries %>% ggplot(
  aes(
    x = estimate_anonymous, xmin = conf.low_anonymous, xmax = conf.high_anonymous,
    y = estimate_revealed, ymin = conf.low_revealed, ymax = conf.high_revealed,
    colour = site
  )
) +
  geom_errorbar(orientation = "x", width = 0, linewidth = 1) +
  geom_errorbar(orientation = "y", width = 0, linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 3) +
  geom_text_repel(aes(label = str_to_upper(site)), box.padding = 0.7, segment.colour = NA) +
  labs(
    x = "Anonymous (Mean Percentage Given ± 95% CI)",
    y = "Revealed (Mean Percentage Given ± 95% CI)"
  ) +
  theme(
    legend.position = "none",
  ) +
  scale_colour_manual(values = get_site_colours()) +
  scale_x_continuous(limits = c(0.04, 0.45), labels = scales::percent) +
  scale_y_continuous(limits = c(0.04, 0.45), labels = scales::percent)

# Combine histograms and Ligon-Schehter ----------------------------------------
combo <- wrap_elements(p_hists) + wrap_elements(p_ls) +
  plot_layout(widths = c(4,3)) +
  plot_annotation(
    title = "Distribution of allocations",
  )

combo
