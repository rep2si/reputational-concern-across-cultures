# Prior sensitivity check: compares fixed effects from models fit with the
# default priors and with wider priors, with and without demographic
# covariates. One plot per outcome model (v3, v3_gv50mr, v3_gv0).

library(tidyverse)
library(brms)

set_theme(cowplot::theme_cowplot())

sharepoint_path <- Sys.getenv("REP2SI_DIR")
model_dir <- file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models")
plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

# Fixed effects shared by cd and cdagw (see cm in coef_table.R). Phi only exists
# in the ordered beta model, and is dropped for the logistic ones.
params <- c(
  "Intercept" = "Intercept",
  "Conditionrevealed" = "Revealed",
  "modist_ordidEQdistmo" = "Distance",
  "moego_cent_decileidEQegomo" = "Donor centrality",
  "moalter_cent_decileidEQaltermo" = "Recipient centrality",
  "modist_ordidEQdistmo:Conditionrevealed" = "Revealed X Distance",
  "moego_cent_decileidEQegomo:Conditionrevealed" = "Revealed X Donor Centrality",
  "moalter_cent_decileidEQaltermo:Conditionrevealed" = "Revealed X Recipient Centrality",
  "anon_first" = "Anonymous First",
  "Conditionrevealed:anon_first" = "Revealed X Anonymous First",
  "phi_Intercept" = "Phi"
)

# Load the four variants of one model set (default or wide prior, with or
# without demographic covariates) and return their coefficients (median, 95% CI)
# and class-level priors
load_set <- function(version) {
  paths <- c(
    cd = "cd.rds",
    cd_wide = "cd-wideprior.rds",
    cdagw = "cdagw.rds",
    cdagw_wide = "cdagw-wideprior.rds"
  )

  # one model at a time to keep memory down
  res <- imap(paths, \(p, m) {
    mod <- readRDS(file.path(model_dir, version, p))
    coefs <- as_tibble(fixef(mod, robust = TRUE), rownames = "term") %>% mutate(model = m)
    priors <- as_tibble(prior_summary(mod)) %>%
      filter(coef == "" | class == "simo", group == "") %>%
      distinct(class, dpar, prior) %>%
      mutate(model = m, version = version)
    rm(mod)
    gc()
    list(coefs = coefs, priors = priors)
  })

  coefs <- map(res, "coefs") %>%
    bind_rows() %>%
    filter(term %in% names(params)) %>%
    mutate(
      dem_covars = str_detect(model, "agw"),
      wide_prior = str_detect(model, "_wide")
    )

  list(coefs = coefs, priors = bind_rows(map(res, "priors")))
}

plot_set_coefs <- function(all_coefs, title) {
  all_coefs %>%
    ggplot(aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = term, colour = wide_prior, group = model, linetype = dem_covars)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_linerange(orientation = "y", position = position_dodge(width = 0.5)) +
    # CI end ticks as points so they stay solid on dashed lines
    geom_point(aes(x = Q2.5), shape = "|", size = 3, position = position_dodge(width = 0.5), show.legend = FALSE) +
    geom_point(aes(x = Q97.5), shape = "|", size = 3, position = position_dodge(width = 0.5), show.legend = FALSE) +
    geom_vline(xintercept = 0, linetype = 3) +
    scale_colour_brewer(type = "qual", palette = 2, labels = c("TRUE" = "Wider prior", "FALSE" = "Weakly regularising prior")) +
    scale_linetype_manual(values = c("FALSE" = "22", "TRUE" = "solid"), labels = c("FALSE" = "Without demographic covariates", "TRUE" = "With demographic covariates")) +
    scale_y_discrete(limits = intersect(rev(names(params)), unique(all_coefs$term)), labels = params) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "vertical"
    ) +
    labs(x = NULL, y = NULL, title = title)
}

# One plot per outcome: model directory -> output file suffix and title
model_sets <- list(
  v3 = c(suffix = "", title = "Share of allocation given (ordered beta)"),
  v3_gv50mr = c(suffix = "-gv50mr", title = "Giving half or more (logistic)"),
  v3_gv0 = c(suffix = "-gv0", title = "Giving exactly 0 (logistic)")
)

all_priors <- list()
for (version in names(model_sets)) {
  ms <- model_sets[[version]]
  set <- load_set(version)
  all_priors[[version]] <- set$priors
  p <- plot_set_coefs(set$coefs, ms[["title"]])
  ggsave(file.path(plots_dir, paste0("prior-sens", ms[["suffix"]], ".pdf")), p, width = 10, height = 6)
}

# Table of priors ---------------------------------------------------------------

# Rows of the table, keyed by brms prior class (phi-specific ones prefixed)
prior_rows <- c(
  "Intercept" = "Intercept",
  "b" = "Coefficients",
  "simo" = "Monotonic effect simplexes",
  "sd" = "Group-level SDs",
  "L" = "Group-level correlations",
  "phi_Intercept" = "Phi: intercept",
  "phi_sd" = "Phi: group-level SDs",
  "cut" = "Cutpoints"
)

# brms prior strings -> something readable in a table
fmt_prior <- function(p) {
  p <- str_replace_all(p, "\\s*,\\s*", ", ")
  case_when(
    p == "" ~ "Flat (improper)",
    str_detect(p, "^induced_dirichlet") ~ paste0("Induced Dirichlet(", str_match(p, "\\[([^]]*)\\]")[, 2], ")"),
    str_detect(p, "^lkj_corr_cholesky") ~ str_replace(p, "^lkj_corr_cholesky", "LKJ"),
    TRUE ~ str_to_sentence(p)
  )
}

priors <- bind_rows(all_priors) %>%
  mutate(
    family = if_else(version == "v3", "ordbeta", "logistic"),
    prior_set = if_else(str_detect(model, "_wide"), "wide", "default"),
    # cutzero and cutone share one induced Dirichlet prior
    key = case_when(
      class %in% c("cutzero", "cutone") ~ "cut",
      dpar == "phi" ~ paste0("phi_", class),
      TRUE ~ class
    ),
    prior = fmt_prior(prior)
  ) %>%
  distinct(family, prior_set, key, prior)

# cd and cdagw, and the two logistic outcomes, should share priors. Stop if not,
# since the table would then hide a difference.
stopifnot(!anyDuplicated(select(priors, family, prior_set, key)))
stopifnot(all(priors$key %in% names(prior_rows)))

prior_tab <- priors %>%
  pivot_wider(names_from = c(family, prior_set), values_from = prior, values_fill = "--") %>%
  mutate(Parameter = prior_rows[key]) %>%
  arrange(match(key, names(prior_rows))) %>%
  select(Parameter, ordbeta_default, ordbeta_wide, logistic_default, logistic_wide)

colnames(prior_tab) <- c("Parameter", "Weakly regularlising", "Wider", "Weakly regularising", "Wider")

prior_tab %>%
  tinytable::tt(caption = "Priors used in the article (weakly regularising) and in the prior sensitivity check (wide)") %>%
  tinytable::group_tt(j = list("Ordered beta (M1)" = 2:3, "Logistic (M2 and M3)" = 4:5)) %>%
  tinytable::theme_latex(outer = "label={tab:priors}") %>%
  tinytable::save_tt(file.path(plots_dir, "prior-sens-priors.tex"), overwrite = TRUE)

# Reweighting the v3 slopes to proper priors -----------------------------------

# ordbetareg ignores coef_prior_SD when phi_reg = "intercepts", so the ordered
# beta models were fitted with flat priors on the slopes (b and bsp). Since the
# fitted prior is flat, the posterior under a Normal(0, s) slope prior is the
# fitted posterior weighted by the Normal(0, s) density of the slope draws. We
# reweight with Pareto-smoothed importance sampling: the default fits to
# Normal(0, 1) and the wide fits to Normal(0, 2), the slope priors that were
# intended for each. Pareto k < 0.7 means the reweighting is reliable.

weighted_quantile <- function(x, w, probs) {
  o <- order(x)
  cw <- cumsum(w[o]) / sum(w)
  map_dbl(probs, \(p) x[o][which(cw >= p)[1]])
}

reweight_slopes <- function(version, file, prior_sd) {
  draws <- as.data.frame(as_draws_df(readRDS(file.path(model_dir, version, file))))
  gc()
  slopes <- str_subset(names(draws), "^(b|bsp)_") %>%
    setdiff(c("b_Intercept", "b_phi_Intercept"))

  log_w <- rowSums(dnorm(as.matrix(draws[, slopes]), 0, prior_sd, log = TRUE))
  ps <- loo::psis(log_w, r_eff = NA)
  w <- as.vector(weights(ps, log = FALSE, normalize = TRUE))

  terms <- intersect(paste0(c("b_", "bsp_"), rep(names(params), each = 2)), slopes)
  map(terms, \(v) {
    x <- draws[[v]]
    tibble(
      term = str_remove(v, "^(b|bsp)_"),
      fitted = quantile(x, c(0.5, 0.025, 0.975), names = FALSE) %>% list(),
      reweighted = weighted_quantile(x, w, c(0.5, 0.025, 0.975)) %>% list()
    )
  }) %>%
    bind_rows() %>%
    mutate(
      file = file,
      prior_sd = prior_sd,
      pareto_k = ps$diagnostics$pareto_k,
      ess = ps$diagnostics$n_eff
    )
}

reweighted <- bind_rows(
  reweight_slopes("v3", "cd.rds", 1),
  reweight_slopes("v3", "cdagw.rds", 1),
  reweight_slopes("v3", "cd-wideprior.rds", 2),
  reweight_slopes("v3", "cdagw-wideprior.rds", 2)
)

fmt_ci <- \(q) sprintf("%.2f [%.2f, %.2f]", q[1], q[2], q[3])
# round before formatting so tiny negative changes do not print as "-0.000"
fmt_change <- \(x) sprintf("%.3f", if_else(abs(x) < 0.0005, 0, x))

rw_tab <- reweighted %>%
  mutate(
    Parameter = params[term],
    `Flat (as fitted)` = map_chr(fitted, fmt_ci),
    `Normal (reweighted)` = map_chr(reweighted, fmt_ci),
    `Change in median` = fmt_change(map_dbl(reweighted, 1) - map_dbl(fitted, 1)),
    group = sprintf(
      "%s, reweighted to Normal(0, %g) (Pareto k = %.2f, ESS = %.0f)",
      c("cd.rds" = "Main model", "cdagw.rds" = "With demographics",
        "cd-wideprior.rds" = "Main model, wide priors",
        "cdagw-wideprior.rds" = "With demographics, wide priors")[file],
      prior_sd, pareto_k, ess
    )
  )

print(select(rw_tab, group, Parameter, `Flat (as fitted)`, `Normal (reweighted)`, `Change in median`), n = Inf, width = Inf)

group_starts <- rw_tab %>%
  mutate(row = row_number()) %>%
  group_by(group) %>%
  summarise(row = min(row)) %>%
  arrange(row)

rw_tab %>%
  select(Parameter, `Flat (as fitted)`, `Normal (reweighted)`, `Change in median`) %>%
  tinytable::tt(caption = "Outcome (1) slopes under the fitted flat prior and reweighted to Normal priors by Pareto-smoothed importance sampling. Posterior medians with 95\\% equal-tailed intervals.") %>%
  tinytable::group_tt(i = setNames(as.list(group_starts$row), group_starts$group)) %>%
  tinytable::style_tt(fontsize = 0.8) %>%
  tinytable::theme_latex(outer = "label={tab:prior_reweight}") %>%
  tinytable::save_tt(file.path(plots_dir, "prior-sens-reweight.tex"), overwrite = TRUE)
