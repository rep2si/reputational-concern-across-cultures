# Mini table of fixed effects against which main hypotheses are tested

library(tidyverse)
library(brms)
library(broom.mixed)
library(modelsummary)
library(emmeans)
library(posterior)
source("./rep2si-common/source_all_submodule.R")

set_theme(cowplot::theme_cowplot())

# Faster than default "parameters()"
options(modelsummary_get = "broom")

sharepoint_path <- Sys.getenv("REP2SI_DIR")
plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

# Parsable names for parameters
cm <- c(
 "b_Conditionrevealed" = "Allocation revealed",
 "bsp_modist_ordidEQdistmo" = "Distance",
 "bsp_moego_cent_decileidEQegomo" = "Donor Centrality",
 "bsp_moalter_cent_decileidEQaltermo" = "Recipient Centrality",
 "bsp_modist_ordidEQdistmo:Conditionrevealed" = "Revealed x Distance",
 "bsp_moego_cent_decileidEQegomo:Conditionrevealed" = "Revealed x Donor Centrality",
 "bsp_moalter_cent_decileidEQaltermo:Conditionrevealed" = "Revealed x Recipient Centrality",
 "b_anon_first" = "Anonymous First",
 "b_Conditionrevealed:anon_first" = "Revealed x Anonymous First"
)

# Load models
m <- readRDS(file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/v3/cd_imput_combined.rds"))
m_gv50mr <- readRDS(file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/v3_gv50mr/cd_imput_combined.rds"))
m_gv0  <- readRDS(file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/v3_gv0/cd_imput_combined.rds"))

tab <- modelsummary(
  list(
    "M1" = m,
    "M2" = m_gv50mr,
    "M3" = m_gv0
  ),
  # effects = "fixed",
  robust = TRUE,
  conf.method = "quantile",
  estimate = "{estimate} [{conf.low}, {conf.high}]",
  coef_map = cm,
  gof_omit = ".*",
  conf_level = 0.95,
  fmt = "%.2f",
  output = "tinytable"
)

tab %>%
  tinytable::style_tt(fontsize = 0.8) %>%
  tinytable::theme_latex(multipage = FALSE, rowhead = 1) %>%
  print("latex")

# Not shown in article---compares fixed effects for models with and without
# demographic covariates. No changes worth noting.

m_cdagw <- readRDS(file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/v3/cdagw_imput_combined.rds"))
m_gv50mr_cdagw <- readRDS(file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/v3_gv50mr/cdagw_imput_combined.rds"))
m_gv0_cdagw  <- readRDS(file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/v3_gv0/cdagw_imput_combined.rds"))

tab_demog <- modelsummary(
  list(
    "M1" = m,
    "M1 demog" = m_cdagw,
    "M2" = m_gv50mr,
    "M2 demog" = m_gv50mr_cdagw,
    "M3" = m_gv0,
    "M3 demog" = m_gv0_cdagw
  ),
  # effects = "fixed",
  robust = TRUE,
  conf.method = "quantile",
  estimate = "{estimate}",
  statistic = "[{conf.low}, {conf.high}]",
  coef_map = cm,
  gof_omit = ".*",
  conf_level = 0.95,
  fmt = "%.2f",
  output = "tinytable"
)
