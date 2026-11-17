library(brms)
library(broom.mixed)
library(modelsummary)
source("./rep2si-common/source_all_submodule.R")

# Faster than default "parameters()"
options(modelsummary_get = "broom")

sharepoint_path <- Sys.getenv("REP2SI_DIR")
plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

cm <- c(
 "b_Intercept" = "Intercept",
 "b_Conditionrevealed" = "Allocation revealed",
 "bsp_moego_cent_decileidEQegomo" = "Donor Centrality",
 "bsp_moalter_cent_decileidEQaltermo" = "Recipient Centrality",
 "bsp_modist_ordidEQdistmo:Conditionrevealed" = "Revealed x Distance",
 "bsp_moego_cent_decileidEQegomo:Conditionrevealed" = "Revealed x Donor Centrality",
 "bsp_moalter_cent_decileidEQaltermo:Conditionrevealed" = "Revealed x Recipient Centrality",
 "b_ego_age_demean_div" = "Donor Age",
 "b_alter_age_demean_div" = "Recipient Age",
 "b_ego_genderMale" = "Donor Male",
 "b_alter_genderMale" = "Recipient Male",
 "b_anon_first" = "Anonymous First",
 "b_ego_demean_logwealth" = "Donor Wealth",
 "b_alter_demean_logwealth" = "Recipient Wealth",
 "b_Conditionrevealed:ego_age_demean_div" = "Revealed x Donor Age",
 "b_Conditionrevealed:alter_age_demean_div" = "Revealed x Recipient Age",
 "b_Conditionrevealed:ego_genderMale" = "Revealed x Donor Male",
 "b_Conditionrevealed:alter_genderMale" = "Revealed x Recipient Male",
 "b_Conditionrevealed:anon_first" = "Revealed x Anonymous First",
 "b_Conditionrevealed:ego_demean_logwealth" = "Revealed x Donor wealth",
 "b_Conditionrevealed:alter_demean_logwealth" = "Revealed x Recipient Wealth",
 "bsp_modist_ordidEQdistmo" = "Distance",
 "b_phi_Intercept" = "Phi"
)

cm_rand <- c(
  # ID
  "sd_ID__Intercept" = "Donor: Intercept",
  "sd_ID__Conditionrevealed" = "Donor: Revealed",
  "sd_ID__alter_age_demean_div" = "Donor: Recipient Age",
  "sd_ID__alter_genderMale" = "Donor: Recipient Gender",
  "sd_ID__alter_demean_logwealth" = "Donor: Recipient Wealth",
  "sd_ID__moalter_cent_decileidEQaltermo" = "Donor: Recipient Centrality",
  "sd_ID__modist_ordidEQdistmo" = "Donor: Distance",
  "sd_ID__Conditionrevealed:alter_age_demean_div" = "Donor: Revealed x Recipient Age",
  "sd_ID__Conditionrevealed:alter_genderMale" = "Donor: Revealed x Recipient Gender",
  "sd_ID__Conditionrevealed:alter_demean_logwealth" = "Donor: Revealed x Recipient Wealth ",
  "sd_ID__moalter_cent_decileidEQaltermo:Conditionrevealed" = "Donor: Revealed x Recipient Centrality",
  # ID_SUID
  "sd_ID_SUID__Intercept" = "Donor Household: Intercept",
  # AID
  "sd_AID__Intercept" = "Recipient: Intercept",
  # ID_SUID
  "sd_AID_SUID__Intercept" = "Recipient Household: Intercept",
  # RID
  "sd_RID__Intercept" = "RA: Intercept",
  # Site
  "sd_site__Intercept" = "Site: Intercept",
  "sd_site__Conditionrevealed" = "Site: Revealed",
  "sd_site__ego_age_demean_div" = "Site: Donor Age",
  "sd_site__alter_age_demean_div" = "Site: Recipient Age",
  "sd_site__ego_genderMale" = "Site: Donor Gender",
  "sd_site__alter_genderMale" = "Site: Recipient Gender",
  "sd_site__ego_demean_logwealth" = "Site: Donor Wealth",
  "sd_site__alter_demean_logwealth" = "Site: Recipient Wealth",
  "sd_site__moego_cent_decileidEQegomo" = "Site: Donor Centrality",
  "sd_site__moalter_cent_decileidEQaltermo" = "Site: Recipient Centrality",
  "sd_site__modist_ordidEQdistmo" = "Site: Distance",
  "sd_site__anon_first" = "Site: Anonymous First",
  "sd_site__Conditionrevealed:ego_age_demean_div" = "Site: Revealed x Donor Age",
  "sd_site__Conditionrevealed:alter_age_demean_div" = "Site: Revealed x Recipient Age",
  "sd_site__Conditionrevealed:ego_genderMale" = "Site: Revealed x Donor Gender",
  "sd_site__Conditionrevealed:alter_genderMale" = "Site: Revealed x Recipient Gender",
  "sd_site__Conditionrevealed:ego_demean_logwealth" = "Site: Revealed x Donor Wealth",
  "sd_site__Conditionrevealed:alter_demean_logwealth" = "Site: Revealed x Recipient Wealth",
  "sd_site__Conditionrevealed:anon_first" = "Site: Revealed x Anonymous First",
  "sd_site__moego_cent_decileidEQegomo:Conditionrevealed" = "Site: Revealed x Donor Centrality",
  "sd_site__moalter_cent_decileidEQaltermo:Conditionrevealed" = "Site: Revealed x Recipient Centrality",
  "sd_site__modist_ordidEQdistmo:Conditionrevealed" = "Site: Revealed x Distance",
  "sd_site__phi_Intercept" = "Site: Phi (Shape of Beta Dist) ",
  # rep dyads
  "sd_dyad_id_reponly__dyad_repeated" = "Repeated Dyad: Intercept"
)

cm_cor <- c(
  "cor_ID__Intercept__Conditionrevealed" = "Donor: Corr Intercept, Revealed",
  "cor_ID__Intercept__alter_age_demean_div" = "Donor: Corr Intercept, Recipient Age",
  "cor_ID__Conditionrevealed__alter_age_demean_div" = "Donor: Corr Revealed, Recipient Age",
  "cor_ID__Intercept__alter_genderMale" = "Donor: Corr Intercept, Recipient Gender",
  "cor_ID__Conditionrevealed__alter_genderMale" = "Donor: Corr Revealed, Recipient Gender",
  "cor_ID__alter_age_demean_div__alter_genderMale" = "Donor: Corr Recipient Age, Recipient Gender",
  "cor_ID__Intercept__alter_demean_logwealth" = "Donor: Corr Intercept, Recipient Wealth",
  "cor_ID__Conditionrevealed__alter_demean_logwealth" = "Donor: Corr Revealed, Recipient Wealth",
  "cor_ID__alter_age_demean_div__alter_demean_logwealth" = "Donor: Corr Recipient Age, Recipient Wealth",
  "cor_ID__alter_genderMale__alter_demean_logwealth" = "Donor: Corr Recipient Gender, Recipient Wealth",
  "cor_ID__Intercept__Conditionrevealed:alter_age_demean_div" = "Donor: Corr Intercept, Revealed x Recipient Age",
  "cor_ID__Conditionrevealed__Conditionrevealed:alter_age_demean_div" = "Donor: Corr Revealed, Revealed x Recipient Age",
  "cor_ID__alter_age_demean_div__Conditionrevealed:alter_age_demean_div" = "Donor: Corr Recipient Age, Revealed x Recipient Age",
  "cor_ID__alter_genderMale__Conditionrevealed:alter_age_demean_div" = "Donor: Corr Recipient Gender, Revealed x Recipient Age",
  "cor_ID__alter_demean_logwealth__Conditionrevealed:alter_age_demean_div" = "Donor: Corr Recipient Wealth, Revealed x Recipient Age",
  "cor_ID__Intercept__Conditionrevealed:alter_genderMale" = "Donor: Corr Intercept, Revealed x Recipient Gender",
  "cor_ID__Conditionrevealed__Conditionrevealed:alter_genderMale" = "Donor: Corr Revealed, Revealed x Recipient Gender",
  "cor_ID__alter_age_demean_div__Conditionrevealed:alter_genderMale" = "Donor: Corr Recipient Age, Revealed x Recipient Gender",
  "cor_ID__alter_genderMale__Conditionrevealed:alter_genderMale" = "Donor: Corr Recipient Gender, Revealed x Recipient Gender",
  "cor_ID__alter_demean_logwealth__Conditionrevealed:alter_genderMale" = "Donor: Corr Recipient Wealth, Revealed x Recipient Gender",
  "cor_ID__Conditionrevealed:alter_age_demean_div__Conditionrevealed:alter_genderMale" = "Donor: Corr Revealed x Recipient Age, Revealed x Recipient Gender",
  "cor_ID__Intercept__Conditionrevealed:alter_demean_logwealth" = "Donor: Corr Intercept, Revealed x Recipient Wealth",
  "cor_ID__Conditionrevealed__Conditionrevealed:alter_demean_logwealth" = "Donor: Corr Revealed, Revealed x Recipient Wealth",
  "cor_ID__alter_age_demean_div__Conditionrevealed:alter_demean_logwealth" = "Donor: Corr Recipient Age, Revealed x Recipient Wealth",
  "cor_ID__alter_genderMale__Conditionrevealed:alter_demean_logwealth" = "Donor: Corr Recipient Gender, Revealed x Recipient Wealth",
  "cor_ID__alter_demean_logwealth__Conditionrevealed:alter_demean_logwealth" = "Donor: Corr Recipient Wealth, Revealed x Recipient Wealth",
  "cor_ID__Conditionrevealed:alter_age_demean_div__Conditionrevealed:alter_demean_logwealth" = "Donor: Corr Revealed x Recipient Age, Revealed x Recipient Wealth",
  "cor_ID__Conditionrevealed:alter_genderMale__Conditionrevealed:alter_demean_logwealth" = "Donor: Corr Revealed x Recipient Gender, Revealed x Recipient Wealth",
  "cor_ID__Intercept__moalter_cent_decileidEQaltermo" = "Donor: Corr Intercept, Recipient Centrality",
  "cor_ID__Conditionrevealed__moalter_cent_decileidEQaltermo" = "Donor: Corr Revealed, Recipient Centrality",
  "cor_ID__alter_age_demean_div__moalter_cent_decileidEQaltermo" = "Donor: Corr Recipient Age, Recipient Centrality",
  "cor_ID__alter_genderMale__moalter_cent_decileidEQaltermo" = "Donor: Corr Recipient Gender, Recipient Centrality",
  "cor_ID__alter_demean_logwealth__moalter_cent_decileidEQaltermo" = "Donor: Corr Recipient Wealth, Recipient Centrality",
  "cor_ID__Conditionrevealed:alter_age_demean_div__moalter_cent_decileidEQaltermo" = "Donor: Corr Revealed x Recipient Age, Recipient Centrality",
  "cor_ID__Conditionrevealed:alter_genderMale__moalter_cent_decileidEQaltermo" = "Donor: Corr Revealed x Recipient Gender, Recipient Centrality",
  "cor_ID__Conditionrevealed:alter_demean_logwealth__moalter_cent_decileidEQaltermo" = "Donor: Corr Revealed x Recipient Wealth, Recipient Centrality",
  "cor_ID__Intercept__modist_ordidEQdistmo" = "Donor: Corr Intercept, Distance",
  "cor_ID__Conditionrevealed__modist_ordidEQdistmo" = "Donor: Corr Revealed, Distance",
  "cor_ID__alter_age_demean_div__modist_ordidEQdistmo" = "Donor: Corr Recipient Age, Distance",
  "cor_ID__alter_genderMale__modist_ordidEQdistmo" = "Donor: Corr Recipient Gender, Distance",
  "cor_ID__alter_demean_logwealth__modist_ordidEQdistmo" = "Donor: Corr Recipient Wealth, Distance",
  "cor_ID__Conditionrevealed:alter_age_demean_div__modist_ordidEQdistmo" = "Donor: Corr Revealed x Recipient Age, Distance",
  "cor_ID__Conditionrevealed:alter_genderMale__modist_ordidEQdistmo" = "Donor: Corr Revealed x Recipient Gender, Distance",
  "cor_ID__Conditionrevealed:alter_demean_logwealth__modist_ordidEQdistmo" = "Donor: Corr Revealed x Recipient Wealth, Distance",
  "cor_ID__moalter_cent_decileidEQaltermo__modist_ordidEQdistmo" = "Donor: Corr Recipient Centrality, Distance",
  "cor_ID__Intercept__moalter_cent_decileidEQaltermo:Conditionrevealed" = "Donor: Corr Intercept, Revealed x Recipient Centrality",
  "cor_ID__Conditionrevealed__moalter_cent_decileidEQaltermo:Conditionrevealed" = "Donor: Corr Revealed, Revealed x Recipient Centrality",
  "cor_ID__alter_age_demean_div__moalter_cent_decileidEQaltermo:Conditionrevealed" = "Donor: Corr Recipient Age, Revealed x Recipient Centrality",
  "cor_ID__alter_genderMale__moalter_cent_decileidEQaltermo:Conditionrevealed" = "Donor: Corr Recipient Gender, Revealed x Recipient Centrality",
  "cor_ID__alter_demean_logwealth__moalter_cent_decileidEQaltermo:Conditionrevealed" = "Donor: Corr Recipient Wealth, Revealed x Recipient Centrality",
  "cor_ID__Conditionrevealed:alter_age_demean_div__moalter_cent_decileidEQaltermo:Conditionrevealed" = "Donor: Corr Revealed x Recipient Age, Revealed x Recipient Centrality",
  "cor_ID__Conditionrevealed:alter_genderMale__moalter_cent_decileidEQaltermo:Conditionrevealed" = "Donor: Corr Revealed x Recipient Gender, Revealed x Recipient Centrality",
  "cor_ID__Conditionrevealed:alter_demean_logwealth__moalter_cent_decileidEQaltermo:Conditionrevealed" = "Donor: Corr Revealed x Recipient Wealth, Revealed x Recipient Centrality",
  "cor_ID__moalter_cent_decileidEQaltermo__moalter_cent_decileidEQaltermo:Conditionrevealed" = "Donor: Corr Recipient Centrality, Revealed x Recipient Centrality",
  "cor_ID__modist_ordidEQdistmo__moalter_cent_decileidEQaltermo:Conditionrevealed" = "Donor: Corr Distance, Revealed x Recipient Centrality"
)

# This will not return correct values for imputed models, because uses the first
# imputed dataset only --- don't forget to update manually
add_r2 <- function(m, nd = 1000, seed = 2025) {
  set.seed(seed)
  r2_marginal    <- brms::bayes_R2(m, re_formula = NA, summary = FALSE, ndraws = nd)
  r2_conditional <- brms::bayes_R2(m, summary = FALSE, ndraws = nd)

  eti_marginal    <- bayestestR::eti(r2_marginal, ci = 0.95)
  eti_conditional <- bayestestR::eti(r2_conditional, ci = 0.95)

  m$r2_cache <- list(
    marginal = list(
      median = median(r2_marginal),
      lower = eti_marginal$CI_low,
      upper = eti_marginal$CI_high
    ),
    conditional = list(
      median = median(r2_conditional),
      lower = eti_conditional$CI_low,
      upper = eti_conditional$CI_high
    )
  )
  m
}

glance_custom.brmsfit <- function(x, ...) {
  g   <- brms::ngrps(x)
  get <- function(nm) if (is.null(g[[nm]])) NA_integer_ else g[[nm]]

  fmt_r2 <- function(r2) {
    if (is.null(r2)) return(NA_character_)
    sprintf("%.2f [%.2f, %.2f]", r2$median, r2$lower, r2$upper)
  }

  data.frame(
    `N donors`     = get("ID"),
    `N recipients` = get("AID"),
    `N RAs`        = get("RID"),
    `N sites`      = get("site"),
    `R2 marginal`    = fmt_r2(x$r2_cache$marginal),
    `R2 conditional` = fmt_r2(x$r2_cache$conditional),
    check.names = FALSE
  )
}

registerS3method("glance_custom", "brmsfit", glance_custom.brmsfit,
                 envir = asNamespace("modelsummary"))


# Could write a function, but am lazy
make_tab <- function(model_dir, caption, out, lab_add) {

  message(caption)

  message("loading models")
  m_0 <- readRDS(file.path(model_dir, "0.rds"))
  m_c <- readRDS(file.path(model_dir, "c.rds"))
  m_d <- readRDS(file.path(model_dir, "d.rds"))
  m_cd <- readRDS(file.path(model_dir, "cd.rds"))
  m_cd_imp <- readRDS(file.path(model_dir, "cd_imput_combined.rds"))
  m_cdagw <- readRDS(file.path(model_dir, "cdagw.rds"))
  m_cdagw_imp <- readRDS(file.path(model_dir, "cdagw_imput_combined.rds"))

  message("Adding R2. Takes a while")
  mods <- lapply(list(m_0, m_c, m_d, m_cd, m_cd_imp, m_cdagw, m_cdagw_imp), add_r2)
  names(mods) <- c("Cond only", "Cent", "Dist", "Cent + Dist",
    "Main: Cent + Dist (with imput)", "Cent + Dist + Demogs", "Main + Demogs (with imput)")


  message("preparing table")
  tab <- modelsummary(
    mods,
    # effects     = c("fixed", "ran_pars"),
    # effects     = "fixed",
    robust      = TRUE,          # posterior median
    conf.method = "quantile",    # ETI; "HPDinterval" for HDI
    estimate    = "{estimate} [{conf.low}, {conf.high}]",
    statistic = NULL,
    conf_level  = 0.95,
    coef_map    = c(cm, cm_rand, cm_cor),
    # metrics     = "R2",
    gof_map = c("nobs", "N donors", "N recipients", "N RAs", "N sites",
      "R2 marginal", "R2 conditional"),
    width       = c(0.3, rep(0.10, 7)),
    fmt         = 2,
    output      = "tinytable",
    title = caption
    # notes = list(
    #   "Posterior medians with 95% equal-tailed credible intervals in brackets.",
    #   "Group-level rows are posterior SDs; correlation rows are group-level correlations for the donor (ID) random effects.",
    #   # "Fitted in brms 4 chains x 2000 iterations (1000 warmup). All R-hat < 1.01; minimum bulk-ESS = X; minimum tail-ESS = Y; Z divergent transitions.",
    #   "Blank cells indicate the term was not in that model."
    # )
  )

  n_fixed <- sum(cm %in% tab@data[[1]])

  message("writing out table")
  tab |>
    tinytable::group_tt(i = setNames(list(1, n_fixed + 1),
      c("Population-level effects", "Group-level effects"))) |>
    tinytable::style_tt(fontsize = 0.6) |>
    tinytable::theme_latex(multipage = TRUE, rowhead = 1, outer = lab_add) |>
    tinytable::save_tt(out, overwrite = TRUE)
}

all_models <- list(
  c(
    path = file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/v3"),
    capt = "Ordered Beta Regression for outcome (1), amount given",
    output = file.path(plots_dir, "coef_table.tex"),
    lab = "label={tab:coefs_main}"
  ),
  c(
    path = file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/v3_gv50mr"),
    capt = "Logistic regression for outcome (2), giving half or more",
    output = file.path(plots_dir, "coef_table_gv50mr.tex"),
    lab = "label={tab:coefs_gv50mr}"
  ),
  c(
    path = file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/v3_gv0"),
    capt = "Logistic regression for outcome (3), giving exactly 0",
    output = file.path(plots_dir, "coef_table_gv0.tex"),
    lab = "label={tab:coefs_gv0}"
  )
)

map(all_models, \(x) make_tab(x[[1]], x[[2]], x[[3]], x[[4]]))
