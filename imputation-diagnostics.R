library(tidyverse)
library(brms)
library(ordbetareg)
library(posterior)
library(knitr)
library(kableExtra)

set.seed(4818)

sharepoint_path <- Sys.getenv("REP2SI_DIR")
plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

# Use cmdstanr -- make sure this is installed on slurm
options(brms.backend = "cmdstanr")

# Helper function to extract diagnostics from a single brms model
get_model_diagnostics <- function(model_path) {
  model <- readRDS(model_path)
  np <- nuts_params(model)
  s <- summary(model)

  n_divergent <- sum(np$Value[np$Parameter == "divergent__"])
  n_max_treedepth <- sum(np$Value[np$Parameter == "treedepth__"] >= model$fit@sim$max_td)

  rhats_fixed <- s$fixed$Rhat
  rhats_random <- bind_rows(s$random)$Rhat

  tibble(
    model_name = basename(model_path),
    n_divergent = n_divergent,
    max_rhat_fixed = max(rhats_fixed, na.rm = TRUE),
    max_rhat_random = max(rhats_random, na.rm = TRUE),
  )
}

# Names follow the paper: no suffix -> M1, _gv50mr -> M2, _gv0 -> M3;
# cdagw adds demographic covariates.
models_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models")
dir_list <- c(
  "M1" = file.path(models_dir, "v3", "cd_imput"),
  "M1 (with demog covars)" = file.path(models_dir, "v3", "cdagw_imput"),
  "M2" = file.path(models_dir, "v3_gv50mr", "cd_imput"),
  "M2 (with demog covars)" = file.path(models_dir, "v3_gv50mr", "cdagw_imput"),
  "M3" = file.path(models_dir, "v3_gv0", "cd_imput"),
  "M3 (with demog covars)" = file.path(models_dir, "v3_gv0", "cdagw_imput")
)

models <- map(dir_list, \(x) list.files(x, pattern = "\\.rds$", full.names = TRUE))

diags <- map(models, \(list) map(list, get_model_diagnostics))
diags <- map(diags, bind_rows)
diags <- map(diags, \(d) {
  d %>%
  mutate(model_name = as.numeric(str_extract(model_name, "\\d+"))) %>%
  arrange(model_name)
})

diags

# One table for all model sets, with a subtitle per set.
kbl_out <- diags %>%
  bind_rows() %>%
  kbl(format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      digits = 2,
      col.names = c("Imputation", "Divergent transitions", "Max Rhat (fixed)", "Max Rhat (random)"),
      caption = "Sampler diagnostics for the models fitted on each imputed dataset. The high Rhat values in these models results from the reduced number of warmup iterations (500) and the low numbr of post-warmup posterior draws (100 per chain). The model resulting from combining posterior draws from all imputations is still adequate, as explained above. Note also that all our hypotheses are tested against fixed effects.",
      label = "imputation-diags") %>%
  kable_styling(latex_options = c("striped", "repeat_header")) %>%
  pack_rows(index = map_int(diags, nrow))

writeLines(kbl_out, file.path(plots_dir, "imputation_diags.tex"))
