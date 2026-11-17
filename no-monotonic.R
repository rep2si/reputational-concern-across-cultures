library(tidyverse)
library(brms)
library(broom.mixed)
library(marginaleffects)
source("./rep2si-common/source_all_submodule.R")

set_theme(cowplot::theme_cowplot())

plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

# Faster than default "parameters()"
options(modelsummary_get = "broom")

# Mini table of fixed effects and pd for 3 main models. In an ideal world would
# centre the mo() values, etc.

sharepoint_path <- Sys.getenv("REP2SI_DIR")
plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

m <- readRDS(file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v3_flexord", "cd.rds"))


plot_predictions(
  m,
  newdata = datagrid(ego_cent_decile_cat = unique,
                     grid_type = "counterfactual"),
  by = "ego_cent_decile_cat",
  ndraws = 50
) + labs(
  y = "Allocation",
  x = "Donor centrality"
)

ggsave(file.path(plots_dir, "donor_alloc-nomono.pdf"))

plot_predictions(
  m,
  newdata = datagrid(alter_cent_decile_cat = unique,
                     grid_type = "counterfactual"),
  by = "alter_cent_decile_cat",
  ndraws = 50
) + labs(
  y = "Allocation",
  x = "Alter centrality"
)

ggsave(file.path(plots_dir, "alter_alloc-nomono.pdf"))

plot_comparisons(
  m,
  variables = list(Condition = c("anonymous", "revealed")),
  newdata = datagrid(
    ego_cent_decile_cat = unique,
    grid_type = "counterfactual"
  ),
  by = "ego_cent_decile_cat",
  ndraws = 50
) + labs(
  y = "Effect of condition",
  x = "Donor centrality"
)

ggsave(file.path(plots_dir, "donor_effect-nomono.pdf"))

plot_comparisons(
  m,
  variables = list(Condition = c("anonymous", "revealed")),
  newdata = datagrid(
    alter_cent_decile_cat = unique,
    grid_type = "counterfactual"
  ),
  by = "alter_cent_decile_cat",
  ndraws = 50
) + labs(
  y = "Effect of condition",
  x = "Recipient centrality"
)

ggsave(file.path(plots_dir, "alter_effect-nomono.pdf"))
