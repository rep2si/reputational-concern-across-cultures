# Generates the plot comparing imputed and non-imputed models, with and
# without demographic covariates. Used in the article to show that all lead to effectively identical
# conclusions. One plot per outcome model (v3, v3_gv50mr, v3_gv0).

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
  "moalter_cent_decileidEQaltermo:Conditionrevealed" = "Revealed X Recipient Centrality"
)

# Load the four variants of one model set (imputed or not, with or without
# demographic covariates) and return their coefficients (median, 95% CI)
get_set_coefs <- function(version) {
  paths <- c(
    cdagw = "cdagw.rds",
    cd = "cd.rds",
    cdagw_imp = "cdagw_imput_combined.rds",
    cd_imp = "cd_imput_combined.rds"
  )

  # one model at a time to keep memory down
  imap(paths, \(p, m) {
    coefs <- fixef(readRDS(file.path(model_dir, version, p)), robust = TRUE, pars = names(params))
    gc()
    as_tibble(coefs, rownames = "term") %>% mutate(model = m)
  }) %>%
    bind_rows() %>%
    mutate(
      dem_covars = str_detect(model, "agw"),
      imputed = str_detect(model, "_imp")
    )
}

plot_set_coefs <- function(all_coefs) {
  all_coefs %>%
    ggplot(aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = term, colour = imputed, group = model, linetype = dem_covars)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_linerange(orientation = "y", position = position_dodge(width = 0.5)) +
    # CI end ticks as points so they stay solid on dashed lines
    geom_point(aes(x = Q2.5), shape = "|", size = 3, position = position_dodge(width = 0.5), show.legend = FALSE) +
    geom_point(aes(x = Q97.5), shape = "|", size = 3, position = position_dodge(width = 0.5), show.legend = FALSE) +
    geom_vline(xintercept = 0, linetype = 3) +
    scale_colour_brewer(type = "qual", palette = 2, labels = c("TRUE" = "Imputed", "FALSE" = "Not Imputed")) +
    scale_linetype_manual(values = c("FALSE" = "22", "TRUE" = "solid"), labels = c("FALSE" = "Without demographic covariates", "TRUE" = "With demographic covariates")) +
    scale_y_discrete(limits = rev(names(params)), labels = params) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "vertical"
    ) +
    labs(x = NULL, y = NULL)
}

# One plot per outcome: model directory -> output file suffix
model_sets <- c(
  v3 = "",
  v3_gv50mr = "-gv50mr",
  v3_gv0 = "-gv0"
)

for (version in names(model_sets)) {
  p <- plot_set_coefs(get_set_coefs(version))
  ggsave(file.path(plots_dir, paste0("coefs-models", model_sets[[version]], ".pdf")), p, width = 10, height = 6)
}
