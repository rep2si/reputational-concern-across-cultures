library(tidyverse)
library(marginaleffects)
library(brms)
set_theme(cowplot::theme_cowplot())

sharepoint_path <- Sys.getenv("REP2SI_DIR")
plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

m <- readRDS(file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/comp/4_levels_componly.rds"))

# Relative frequency of comprehension levels, one row per donor
m$data %>%
  distinct(site, ID, ego_comp) %>%
  count(site, ego_comp) %>%
  group_by(site) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = ego_comp, y = prop)) +
  geom_col() +
  facet_grid(cols = vars(site)) +
  labs(
    x = "Comprehension level",
    y = "Proportion of donors"
  ) +
  theme(panel.spacing = unit(2, "lines"))

ggsave(file.path(plots_dir, "comp-freq.pdf"), width = 8, height = 4)

# Effect of revealing the allocation (revealed - anonymous), by comprehension level
cond_eff <- avg_comparisons(m, variables = "Condition", by = c("site", "ego_comp"))

cond_eff %>%
  ggplot(aes(y = ego_comp, x = estimate, xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange() +
  facet_grid(rows = vars(site)) +
  labs(
    x = "Effect of revealing allocation (revealed - anonymous)",
    y = "Comprehension level"
  ) +
  theme(panel.spacing = unit(2, "lines"))

ggsave(file.path(plots_dir, "comp-cond-effect.pdf"), width = 8, height = 4)

# Does the Condition effect differ across comprehension levels?
# Difference in effect relative to "fail", within site
avg_comparisons(
  m,
  variables = "Condition",
  by = c("site", "ego_comp"),
  hypothesis = difference ~ reference | site
)
