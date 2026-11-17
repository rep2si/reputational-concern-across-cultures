# Generates the plots comparing models that keep vs drop participants who
# failed the comprehension checks, without and with demographic covariates.
# Used in the article to show that both lead to effectively identical
# conclusions.

library(tidyverse)
library(brms)

set_theme(cowplot::theme_cowplot())

sharepoint_path <- Sys.getenv("REP2SI_DIR")
model_dir <- file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models")
plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

params <- c(
  "Conditionrevealed" = "Revealed",
  "modist_ordidEQdistmo" = "Distance",
  "moego_cent_decileidEQegomo" = "Donor centrality",
  "moalter_cent_decileidEQaltermo" = "Recipient centrality",
  "modist_ordidEQdistmo:Conditionrevealed" = "Revealed X Distance",
  "moego_cent_decileidEQegomo:Conditionrevealed" = "Revealed X Donor Centrality",
  "moalter_cent_decileidEQaltermo:Conditionrevealed" = "Revealed X Recipient Centrality",
  "anon_first" = "Anonymous first",
  "Conditionrevealed:anon_first" = "Revealed X Anonymous First"
)

params_demog <- c(
  params,
  "ego_age_demean_div" = "Donor age",
  "alter_age_demean_div" = "Recipient age",
  "ego_genderMale" = "Donor male",
  "alter_genderMale" = "Recipient male",
  "ego_demean_logwealth" = "Donor wealth",
  "alter_demean_logwealth" = "Recipient wealth",
  "Conditionrevealed:ego_age_demean_div" = "Revealed X Donor Age",
  "Conditionrevealed:alter_age_demean_div" = "Revealed X Recipient Age",
  "Conditionrevealed:ego_genderMale" = "Revealed X Donor Male",
  "Conditionrevealed:alter_genderMale" = "Revealed X Recipient Male",
  "Conditionrevealed:ego_demean_logwealth" = "Revealed X Donor Wealth",
  "Conditionrevealed:alter_demean_logwealth" = "Revealed X Recipient Wealth"
)

# Load each model (path relative to model_dir -> legend label) and plot its
# coefficients (median, 95% CI)
plot_comparison <- function(model_paths, params) {
  # one model at a time to keep memory down
  all_coefs <- imap(model_paths, \(label, p) {
    coefs <- fixef(readRDS(file.path(model_dir, p)), robust = TRUE, pars = names(params))
    gc()
    as_tibble(coefs, rownames = "term") %>% mutate(model = label)
  }) %>%
    bind_rows() %>%
    # reversed so the dodge puts the first model on top within each term
    mutate(model = factor(model, levels = rev(model_paths)))

  all_coefs %>%
    ggplot(aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = term, colour = model)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_linerange(orientation = "y", position = position_dodge(width = 0.5)) +
    geom_point(aes(x = Q2.5), shape = "|", size = 3, position = position_dodge(width = 0.5), show.legend = FALSE) +
    geom_point(aes(x = Q97.5), shape = "|", size = 3, position = position_dodge(width = 0.5), show.legend = FALSE) +
    geom_vline(xintercept = 0, linetype = 3) +
    scale_colour_brewer(type = "qual", palette = 2, guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev(names(params)), labels = params) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "vertical"
    ) +
    labs(x = NULL, y = NULL)
}

p <- plot_comparison(
  c(
    "v3/cd.rds" = "All participants",
    "v3_compfilt/cd_comp_lenient.rds" = "Failed comprehension checks dropped"
  ),
  params
)
ggsave(file.path(plots_dir, "coefs-models-compfilt.pdf"), p, width = 10, height = 6)

p_demog <- plot_comparison(
  c(
    "v3/cdagw.rds" = "All participants",
    "v3_compfilt/cdagw_comp_lenient.rds" = "Failed comprehension checks dropped"
  ),
  params_demog
)
ggsave(file.path(plots_dir, "coefs-models-compfilt-demog.pdf"), p_demog, width = 10, height = 12)
