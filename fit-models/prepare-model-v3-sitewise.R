library(tidyverse)
library(brms)
library(ordbetareg)

# Not necessary, but short and helpful while working on the script
if(Sys.info()["nodename"] == "clustername") {
  d <- readRDS("/home/xxx/rep2si/allocs.rds")
} else {
  # Read from appropriate path when testing locally
  sharepoint_path <- Sys.getenv("REP2SI_DIR")
  d <- readRDS(file.path(sharepoint_path, "rep2si-data/global-data/allocs.rds"))
}

# Formulae ---------------------------------------------------------------------
formulae <- list(
  f_cd = bf(
    percGiven ~
      Condition * mo(dist_ord, id = "distmo") +
        Condition * mo(ego_cent_decile, id = "egomo") +
        Condition * mo(alter_cent_decile, id = "altermo") +
        Condition * anon_first +
        (Condition * mo(alter_cent_decile, id = "altermo") + Condition * mo(dist_ord, id = "distmo") | ID) +
        (1 | AID) +
        (1 | RID) +
        (1 | ID_SUID) +
        (1 | AID_SUID) +
        (0 + dyad_repeated | dyad_id_reponly)
  ),

  f_cdagw = bf(
    percGiven ~
      Condition * mo(dist_ord, id = "distmo") +
        Condition * mo(ego_cent_decile, id = "egomo") +
        Condition * mo(alter_cent_decile, id = "altermo") +
        Condition * ego_age_demean_div +
        Condition * alter_age_demean_div +
        Condition * ego_gender +
        Condition * alter_gender +
        Condition * anon_first +
        Condition * ego_demean_logwealth +
        Condition * alter_demean_logwealth +
        (Condition * mo(alter_cent_decile, id = "altermo") + Condition * mo(dist_ord, id = "distmo") + Condition * alter_age_demean_div + Condition * alter_gender + Condition * alter_demean_logwealth | ID) +
        (1 | AID) +
        (1 | RID) +
        (1 | ID_SUID) +
        (1 | AID_SUID) +
        (0 + dyad_repeated | dyad_id_reponly)
  )
)

# Construct equivalent models for other outcomes
formulae_gv0 <- map(formulae, \(f) {
  update(bf(f$formula), formula = give0 ~ .)
})

formulae_gv50 <- map(formulae, \(f) {
  update(bf(f$formula), formula = give50 ~ .)
})

formulae_gv50mr <- map(formulae, \(f) {
  update(bf(f$formula), formula = give50ormore ~ .)
})

# Create datframes --------------------------------------------------------------

# Strip f_ from names
# Note that we could just not use them in the first place but that throws me off
simple_names <- unlist(
  map(names(formulae), \(name) {str_replace(name, "^f_", "")})
)

all_sites <- c("te", "wl", "wh", "tz", "az", "tm", "np", "up")

models <- tibble(
  outcome = "percGiven",
  model_name = simple_names,
  formula = formulae,
  phireg = "intercepts",
  comprehension = "all",
  base_dir = "v3/sitewise",
)
models <- expand_grid(models, site = all_sites)

models_gv0 <- tibble(
  outcome = "give0",
  model_name = simple_names,
  formula = formulae_gv0,
  phireg = "intercepts",
  comprehension = "all",
  base_dir = "v3_gv0/sitewise",
)
models_gv0 <- expand_grid(models_gv0, site = all_sites)

models_gv50mr <- tibble(
  outcome = "give50ormore",
  model_name = simple_names,
  formula = formulae_gv50mr,
  phireg = "intercepts",
  comprehension = "all",
  base_dir = "v3_gv50mr/sitewise",
)
models_gv50mr <- expand_grid(models_gv50mr, site = all_sites)

v3_models_sitewise <- bind_rows(
  models,
  models_gv0,
  models_gv50mr,
)

# useful when writing sbatch file
v3_models_sitewise %>% print(n = 500)

write_rds(v3_models_sitewise, "/home/xxx/rep2si/v3_models_sitewise.rds")
