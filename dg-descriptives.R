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

# Table 2 ----------------------------------------------------------------------
sites <- split(d, d$site)
sites$total <- d_tot

sites <- map(sites, \(x) {
  ll <- list()
  ll$allocations <- x
  return(ll)
} )


allocations_df <- sites %>%
  imap_dfr(~ {
    alloc <- .x$allocations

    tibble(
      site = .y,
      # Overall allocation stats
      mean_percGiven = mean(alloc$percGiven, na.rm = TRUE),
      sd_percGiven = sd(alloc$percGiven, na.rm = TRUE),
      perc_optout = mean(alloc$optedOut, na.rm = TRUE) * 100,
      # By condition stats
      mean_percGiven_anon = mean(alloc$percGiven[alloc$Condition == "anonymous"], na.rm = TRUE),
      mean_percGiven_rev = mean(alloc$percGiven[alloc$Condition == "revealed"], na.rm = TRUE),
      delta_mean_percGiven = mean_percGiven_rev - mean_percGiven_anon,
      mean_increase_rev = mean_percGiven_rev / mean_percGiven_anon,
      mean_increase_percent = 100 * (mean_increase_rev - 1),
      sd_percGiven_anon = sd(alloc$percGiven[alloc$Condition == "anonymous"], na.rm = TRUE),
      sd_percGiven_rev = sd(alloc$percGiven[alloc$Condition == "revealed"], na.rm = TRUE),
      # Pairing and social connection stats
      perc_rev = mean(alloc$Condition == "revealed", na.rm = TRUE) * 100,
      perc_alter_top = mean(alloc$alter_cent_istop, na.rm = TRUE) * 100,
      perc_alter_low = mean(alloc$alter_cent_islow, na.rm = TRUE) * 100,
      perc_alter_rank_below_0_1 = mean(alloc$alter_cent_rank <= 0.1, na.rm = TRUE) * 100,
      perc_alter_rank_mid = mean(alloc$alter_cent_rank > 0.45 & alloc$alter_cent_rank <= 0.55, na.rm = TRUE) * 100,
      perc_alter_rank_above_0_9 = mean(alloc$alter_cent_rank > 0.9, na.rm = TRUE) * 100,
      perc_blood = mean(alloc$blood_rel, na.rm = TRUE) * 100,
      perc_aff = mean(alloc$aff_rel, na.rm = TRUE) * 100,
      perc_na_net = mean(is.na(alloc$dist)) * 100,
      mean_dist = mean(alloc$dist[is.finite(alloc$dist)], na.rm = TRUE),
      max_dist = max(alloc$dist, na.rm = TRUE),
      mean_neigh = mean(alloc$co_neigh, na.rm = TRUE),
      max_neigh = max(alloc$co_neigh, na.rm = TRUE),
      mean_resist = mean(alloc$resist[is.finite(alloc$resist)], na.rm = TRUE),
      max_resist = max(alloc$resist[is.finite(alloc$resist)], na.rm = TRUE),
    )
  })

# Round selected columns
allocations_df <- allocations_df %>%
  mutate(across(where(is.numeric) & !matches("max_dist|max_neigh"), ~ round(.x, 2)))

knitr::kable(t(allocations_df)[-1,],
             col.names = t(allocations_df)[1,],
             caption = "Allocations Descriptives")


# Filter to the variables we actually want to show
allocations_df_simple <- allocations_df %>%
  select(
    site,
    mean_percGiven_anon,
    sd_percGiven_anon,
    mean_percGiven_rev,
    sd_percGiven_rev,
    perc_na_net,
    mean_dist,
    max_dist,
    perc_alter_rank_below_0_1,
    perc_alter_rank_mid,
    perc_alter_rank_above_0_9
  ) %>%
  mutate(
    site = str_to_upper(site),
    site = ifelse(site == "TOTAL", "Total", site)
  ) %>%
  rename(
    `Mean % Given Anon.` = mean_percGiven_anon,
    `SD  % Given Anon.` = sd_percGiven_anon,
    `Mean % Given Reveal.` = mean_percGiven_rev,
    `SD  % Given Reveal.` = sd_percGiven_rev,
    `% Missing Net` = perc_na_net,
    `Mean Dist.` = mean_dist,
    `Max Dist` = max_dist,
    `% Alter Cent. Bottom 10%` = perc_alter_rank_below_0_1,
    `% Alter Cent. Mid 10%` = perc_alter_rank_mid,
    `% Alter Cent. Top 10%` = perc_alter_rank_above_0_9
  )


knitr::kable(
  t(allocations_df_simple)[-1,],
  col.names = t(allocations_df)[1,],
  format = "latex"
)
