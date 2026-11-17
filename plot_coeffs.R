library(tidyverse)
library(brms)
library(tidybayes)
library(posterior)
library(patchwork)
source("./rep2si-common/source_all_submodule.R")

set_theme(cowplot::theme_cowplot())

sharepoint_path <- Sys.getenv("REP2SI_DIR")
model_dir <- file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/v3")
model_dir_gv0 <- file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/v3_gv0")
model_dir_gv50mr <- file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/v3_gv50mr")
plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

m_cd <- readRDS(file.path(model_dir, "cd_imput_combined.rds"))
draws_cd <- as_draws_df(m_cd)

m_cdagw <- readRDS(file.path(model_dir, "cdagw_imput_combined.rds"))
draws_cdagw <- as_draws_df(m_cdagw)

m_cd_gv50mr <- readRDS(file.path(model_dir_gv50mr, "cd_imput_combined.rds"))
draws_cd_gv50mr <- as_draws_df(m_cd_gv50mr)

m_cdagw_gv50mr <- readRDS(file.path(model_dir_gv50mr, "cdagw_imput_combined.rds"))
draws_cdagw_gv50mr <- as_draws_df(m_cdagw_gv50mr)

m_cd_gv0 <- readRDS(file.path(model_dir_gv0, "cd_imput_combined.rds"))
draws_cd_gv0 <- as_draws_df(m_cd_gv0)

m_cdagw_gv0 <- readRDS(file.path(model_dir_gv0, "cdagw_imput_combined.rds"))
draws_cdagw_gv0 <- as_draws_df(m_cdagw_gv0)

params <- c(
  "b_Conditionrevealed" = "Revealed",
  "bsp_modist_ordidEQdistmo" = "Distance",
  "bsp_moego_cent_decileidEQegomo" = "Donor centrality",
  "bsp_moalter_cent_decileidEQaltermo" = "Recipient centrality",
  "bsp_modist_ordidEQdistmo:Conditionrevealed" = "Revealed X Distance",
  "bsp_moego_cent_decileidEQegomo:Conditionrevealed" = "Revealed X Donor Centrality",
  "bsp_moalter_cent_decileidEQaltermo:Conditionrevealed" = "Revealed X Recipient Centrality"
)

demog_params <- c(
  "b_ego_age_demean_div" = "Donor Age",
  "b_alter_age_demean_div" = "Recipient Age",
  "b_ego_genderMale" = "Donor Male",
  "b_alter_genderMale" = "Recipient Male",
  "b_ego_demean_logwealth" = "Donor Wealth",
  "b_alter_demean_logwealth" = "Recipient Wealth",
  "b_Conditionrevealed:ego_age_demean_div" = "Revealed x Donor Age",
  "b_Conditionrevealed:alter_age_demean_div" = "Revealed x Recipient Age",
  "b_Conditionrevealed:ego_genderMale" = "Revealed x Donor Male",
  "b_Conditionrevealed:alter_genderMale" = "Revealed x Recipient Male",
  "b_Conditionrevealed:ego_demean_logwealth" = "Revealed x Donor Wealth",
  "b_Conditionrevealed:alter_demean_logwealth" = "Revealed x Wealth"
)

# S6 ---------------------------------------------------------------------------
# Combo three models
p1 <- draws_cd %>%
  get_all_site_coefs(names(params)) %>%
  plot_coefs_from_draws() +
  scale_y_discrete(limits = rev(names(params)), labels = params) +
  labs(
    title = "Outcome: share of allocation given.",
    subtitle = "Ordered Beta Regression"
  ) +
  theme(
    legend.title = element_blank()
  )

p2 <- draws_cd_gv50mr %>%
  get_all_site_coefs(names(params)) %>%
  plot_coefs_from_draws() +
  scale_y_discrete(limits = rev(names(params)), labels = params) +
  labs(
    title = "Outcome: giving half or more",
    subtitle = "Binary logistic regression"
  ) +
  theme(
    legend.title = element_blank()
  ) +
  scale_x_continuous(limits = c(-1.6, 1.6))

p3 <- draws_cd_gv0 %>%
  get_all_site_coefs(names(params)) %>%
  plot_coefs_from_draws() +
  scale_y_discrete(limits = rev(names(params)), labels = params) +
  labs(
    title = "Outcome: giving exactly 0",
    subtitle = "Binary logistic regression"
  ) +
  theme(
    legend.title = element_blank()
  ) +
  scale_x_continuous(limits = c(-1.6, 1.6))

p1 / p2 / p3 + plot_layout(guides = "collect")


# S7, S8, S9: Pairwise, with and without demogs --------------------------------
p1 <- draws_cd %>%
  get_all_site_coefs(names(params)) %>%
  plot_coefs_from_draws() +
  scale_y_discrete(limits = rev(names(params)), labels = params) +
  labs(
    title = "Outcome: share of allocation given.",
    subtitle = "Network parameters only"
  ) +
  theme(
    legend.title = element_blank()
  )

p2 <- draws_cdagw %>%
  get_all_site_coefs(names(params)) %>%
  plot_coefs_from_draws() +
  scale_y_discrete(limits = rev(names(params)), labels = params) +
  labs(
    title = "Outcome: share of allocation given.",
    subtitle = "With demographic controls"
  ) +
  theme(
    legend.title = element_blank()
  )

p1 / p2 + plot_layout(guides = "collect")

p1 <- draws_cd_gv50mr %>%
  get_all_site_coefs(names(params)) %>%
  plot_coefs_from_draws() +
  scale_y_discrete(limits = rev(names(params)), labels = params) +
  labs(
    title = "Outcome: giving half or more.",
    subtitle = "Network parameters only"
  ) +
  theme(
    legend.title = element_blank()
  )

p2 <- draws_cdagw_gv50mr %>%
  get_all_site_coefs(names(params)) %>%
  plot_coefs_from_draws() +
  scale_y_discrete(limits = rev(names(params)), labels = params) +
  labs(
    title = "Outcome: giving half or more.",
    subtitle = "With demographic controls"
  ) +
  theme(
    legend.title = element_blank()
  )

p1 / p2 + plot_layout(guides = "collect")

p1 <- draws_cd_gv0 %>%
  get_all_site_coefs(names(params)) %>%
  plot_coefs_from_draws() +
  scale_y_discrete(limits = rev(names(params)), labels = params) +
  labs(
    title = "Outcome: giving half or more.",
    subtitle = "Network parameters only"
  ) +
  theme(
    legend.title = element_blank()
  )

p2 <- draws_cdagw_gv0 %>%
  get_all_site_coefs(names(params)) %>%
  plot_coefs_from_draws() +
  scale_y_discrete(limits = rev(names(params)), labels = params) +
  labs(
    title = "Outcome: giving half or more.",
    subtitle = "With demographic controls"
  ) +
  theme(
    legend.title = element_blank()
  )

p1 / p2 + plot_layout(guides = "collect")

# S10 Plot of demographic predictors -------------------------------------------
draws_cdagw %>%
  get_all_site_coefs(names(demog_params)) %>%
  plot_coefs_from_draws() +
  scale_y_discrete(limits = rev(names(demog_params)), labels = demog_params) +
  labs(
    title = "Outcome: allocations",
    subtitle = "Demographic predictors, from model also including network measures"
  ) +
  theme(
    legend.title = element_blank()
  )

# and estimates
coef(m_cdagw, groups = "site")
coef(m_cdagw, groups = "site", ndraws = 100)
