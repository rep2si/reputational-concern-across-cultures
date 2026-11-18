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
    f_0 = bf(
      percGiven ~
          Condition * anon_first +
          (Condition | ID) +
          (Condition * anon_first || site) +
          (1 | AID) +
          (1 | RID) +
          (1 | ID_SUID) +
          (1 | AID_SUID) +
          (0 + dyad_repeated | dyad_id_reponly),
      phi ~ (1 | site)
    ),

  f_c = bf(
    percGiven ~
      Condition * mo(ego_cent_decile, id = "egomo") +
        Condition * mo(alter_cent_decile, id = "altermo") +
        Condition * anon_first +
        (Condition * mo(alter_cent_decile, id = "altermo") | ID) +
        (Condition * mo(ego_cent_decile, id = "egomo") + Condition * mo(alter_cent_decile, id = "altermo") + Condition * anon_first || site) +
        (1 | AID) +
        (1 | RID) +
        (1 | ID_SUID) +
        (1 | AID_SUID) +
        (0 + dyad_repeated | dyad_id_reponly),
    phi ~ (1 | site)
  ),

  f_c_Ec.Ac_NI = bf(
    percGiven ~
      Condition * mo(ego_cent_decile, id = "egomo") +
        Condition * mo(alter_cent_decile, id = "altermo") +
        mo(ego_cent_decile, id = "egomo") * mo(alter_cent_decile, id = "altermo") +
        Condition * anon_first +
        (Condition * mo(alter_cent_decile, id = "altermo") | ID) +
        (Condition * mo(ego_cent_decile, id = "egomo") + Condition * mo(alter_cent_decile, id = "altermo") + mo(ego_cent_decile, id = "egomo") * mo(alter_cent_decile, id = "altermo") + Condition * anon_first || site) +
        (1 | AID) +
        (1 | RID) +
        (1 | ID_SUID) +
        (1 | AID_SUID) +
        (0 + dyad_repeated | dyad_id_reponly),
    phi ~ (1 | site)
  ),

  f_d = bf(
    percGiven ~
      Condition * mo(dist_ord, id = "distmo") +
        Condition * anon_first +
        (Condition + mo(dist_ord, id = "distmo") | ID) +
        (Condition * mo(dist_ord, id = "distmo") + Condition * anon_first || site) +
        (1 | AID) +
        (1 | RID) +
        (1 | ID_SUID) +
        (1 | AID_SUID) +
        (0 + dyad_repeated | dyad_id_reponly),
    phi ~ (1 | site)
  ),

  f_cd = bf(
    percGiven ~
      Condition * mo(dist_ord, id = "distmo") +
        Condition * mo(ego_cent_decile, id = "egomo") +
        Condition * mo(alter_cent_decile, id = "altermo") +
        Condition * anon_first +
        (Condition * mo(alter_cent_decile, id = "altermo") + mo(dist_ord, id = "distmo") | ID) +
        (Condition * mo(ego_cent_decile, id = "egomo") + Condition * mo(alter_cent_decile, id = "altermo") + Condition * mo(dist_ord, id = "distmo") + Condition * anon_first || site) +
        (1 | AID) +
        (1 | RID) +
        (1 | ID_SUID) +
        (1 | AID_SUID) +
        (0 + dyad_repeated | dyad_id_reponly),
    phi ~ (1 | site)
  ),

  f_cd_c.d_NI = bf(
    percGiven ~
      # note: have to rearrange order to conform with random effects
      Condition * mo(ego_cent_decile, id = "egomo") +
        Condition * mo(alter_cent_decile, id = "altermo") +
        Condition * mo(dist_ord, id = "distmo") +
        mo(alter_cent_decile, id = "altermo") * mo(dist_ord, id = "distmo") +
        mo(ego_cent_decile, id = "egomo") * mo(dist_ord, id = "distmo") +
        Condition * anon_first +
        (Condition * mo(alter_cent_decile, id = "altermo") + mo(dist_ord, id = "distmo") + mo(alter_cent_decile, id = "altermo") * mo(dist_ord, id = "distmo") | ID) +
        (Condition * mo(ego_cent_decile, id = "egomo") + Condition * mo(alter_cent_decile, id = "altermo") + Condition * mo(dist_ord, id = "distmo") + mo(alter_cent_decile, id = "altermo") * mo(dist_ord, id = "distmo") + mo(ego_cent_decile, id = "egomo") * mo(dist_ord, id = "distmo") + Condition * anon_first  || site) +
        (1 | AID) +
        (1 | RID) +
        (1 | ID_SUID) +
        (1 | AID_SUID) +
        (0 + dyad_repeated | dyad_id_reponly),
    phi ~ (1 | site)
  ),

  f_cd_c.d_NI_Ec.Ac_NI = bf(
    percGiven ~
      Condition * mo(alter_cent_decile, id = "altermo") +
        Condition * mo(ego_cent_decile, id = "egomo") +
        Condition * mo(dist_ord, id = "distmo") +
        mo(alter_cent_decile, id = "altermo") * mo(dist_ord, id = "distmo") +
        mo(ego_cent_decile, id = "egomo") * mo(dist_ord, id = "distmo") +
        mo(alter_cent_decile, id = "altermo") * mo(ego_cent_decile, id = "egomo") +
        Condition * anon_first +
        (Condition * mo(alter_cent_decile, id = "altermo") + mo(dist_ord, id = "distmo") + mo(alter_cent_decile, id = "altermo") * mo(dist_ord, id = "distmo") | ID) +
        (Condition * mo(alter_cent_decile, id = "altermo") + Condition * mo(ego_cent_decile, id = "egomo") + Condition * mo(dist_ord, id = "distmo") + mo(alter_cent_decile, id = "altermo") * mo(dist_ord, id = "distmo") + mo(ego_cent_decile, id = "egomo") * mo(dist_ord, id = "distmo") + mo(alter_cent_decile, id = "altermo") * mo(ego_cent_decile, id = "egomo") + Condition * anon_first || site) +
        (1 | AID) +
        (1 | RID) +
        (1 | ID_SUID) +
        (1 | AID_SUID) +
        (0 + dyad_repeated | dyad_id_reponly),
    phi ~ (1 | site)
  ),

  f_ag_NI= bf(
    percGiven ~
      Condition +
        ego_age_demean_div +
        alter_age_demean_div +
        ego_gender +
        alter_gender +
        Condition * anon_first +
        (alter_age_demean_div + alter_gender | ID) +
        (Condition + ego_age_demean_div + alter_age_demean_div + ego_gender + alter_gender + Condition * anon_first || site) +
        (1 | AID) +
        (1 | RID) +
        (1 | ID_SUID) +
        (1 | AID_SUID) +
        (0 + dyad_repeated | dyad_id_reponly),
    phi ~ (1 | site)
  ),

  f_ag = bf(
    percGiven ~
      Condition * ego_age_demean_div +
        Condition * alter_age_demean_div +
        Condition * ego_gender +
        Condition * alter_gender +
        Condition * anon_first +
        (Condition * alter_age_demean_div + Condition * alter_gender | ID) +
        (Condition * ego_age_demean_div + Condition * alter_age_demean_div + Condition * ego_gender + Condition * alter_gender + Condition * anon_first || site) +
        (1 | AID) +
        (1 | RID) +
        (1 | ID_SUID) +
        (1 | AID_SUID) +
        (0 + dyad_repeated | dyad_id_reponly),
    phi ~ (1 | site)
  ),

  f_agw_NI= bf(
    percGiven ~
      Condition +
        ego_age_demean_div +
        alter_age_demean_div +
        ego_gender +
        alter_gender +
        ego_demean_logwealth +
        alter_demean_logwealth +
        Condition * anon_first +
        (alter_age_demean_div + alter_gender + alter_demean_logwealth | ID) +
        (Condition + ego_age_demean_div + alter_age_demean_div + ego_gender + alter_gender + Condition * anon_first || site) +
        (1 | AID) +
        (1 | RID) +
        (1 | ID_SUID) +
        (1 | AID_SUID) +
        (0 + dyad_repeated | dyad_id_reponly),
    phi ~ (1 | site)
  ),

  f_agw = bf(
    percGiven ~
      Condition * ego_age_demean_div +
        Condition * alter_age_demean_div +
        Condition * ego_gender +
        Condition * alter_gender +
        Condition * ego_demean_logwealth +
        Condition * alter_demean_logwealth +
        Condition * anon_first +
        (Condition * alter_age_demean_div + Condition * alter_gender + Condition * alter_demean_logwealth | ID) +
        (Condition * ego_age_demean_div + Condition * alter_age_demean_div + Condition * ego_gender + Condition * alter_gender + Condition * ego_demean_logwealth + Condition * alter_demean_logwealth + Condition * anon_first || site) +
        (1 | AID) +
        (1 | RID) +
        (1 | ID_SUID) +
        (1 | AID_SUID) +
        (0 + dyad_repeated | dyad_id_reponly),
    phi ~ (1 | site)
  ),

  f_cdag = bf(
    percGiven ~
      Condition * mo(dist_ord, id = "distmo") +
        Condition * mo(ego_cent_decile, id = "egomo") +
        Condition * mo(alter_cent_decile, id = "altermo") +
        Condition * ego_age_demean_div +
        Condition * alter_age_demean_div +
        Condition * ego_gender +
        Condition * alter_gender +
        Condition * anon_first +
        (Condition * mo(alter_cent_decile, id = "altermo") + mo(dist_ord, id = "distmo") + Condition * alter_age_demean_div + Condition * alter_gender | ID) +
        (Condition * mo(ego_cent_decile, id = "egomo") + Condition * mo(alter_cent_decile, id = "altermo") + Condition * mo(dist_ord, id = "distmo") + Condition * ego_age_demean_div + Condition * alter_age_demean_div + Condition * ego_gender + Condition * alter_gender + Condition * anon_first || site) +
        (1 | AID) +
        (1 | RID) +
        (1 | ID_SUID) +
        (1 | AID_SUID) +
        (0 + dyad_repeated | dyad_id_reponly),
    phi ~ (1 | site)
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
        (Condition * mo(alter_cent_decile, id = "altermo") + mo(dist_ord, id = "distmo") + Condition * alter_age_demean_div + Condition * alter_gender + Condition * alter_demean_logwealth | ID) +
        (Condition * mo(ego_cent_decile, id = "egomo") + Condition * mo(alter_cent_decile, id = "altermo") + Condition * mo(dist_ord, id = "distmo") + Condition * ego_age_demean_div + Condition * alter_age_demean_div + Condition * ego_gender + Condition * alter_gender + Condition * ego_demean_logwealth + Condition * alter_demean_logwealth + Condition * anon_first || site) +
        (1 | AID) +
        (1 | RID) +
        (1 | ID_SUID) +
        (1 | AID_SUID) +
        (0 + dyad_repeated | dyad_id_reponly),
    phi ~ (1 | site)
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

models <- tibble(
  site = "total",
  outcome = "percGiven",
  model_name = simple_names,
  formula = formulae,
  phireg = "intercepts",
  comprehension = "all",
  base_dir = "v2",
)

models_compfilt <- tibble(
  site = "total",
  outcome = "percGiven",
  model_name = simple_names,
  formula = formulae,
  phireg = "intercepts",
  comprehension = "lenient",
  base_dir = "v2_compfilt",
)

models_gv0 <- tibble(
  site = "total",
  outcome = "give0",
  model_name = simple_names,
  formula = formulae_gv0,
  phireg = "intercepts",
  comprehension = "all",
  base_dir = "v2_gv0",
)

models_gv50 <- tibble(
  site = "total",
  outcome = "give50",
  model_name = simple_names,
  formula = formulae_gv50,
  phireg = "intercepts",
  comprehension = "all",
  base_dir = "v2_gv50",
)

models_gv50mr <- tibble(
  site = "total",
  outcome = "give50ormore",
  model_name = simple_names,
  formula = formulae_gv50mr,
  phireg = "intercepts",
  comprehension = "all",
  base_dir = "v2_gv50mr",
)

v2_models <- bind_rows(
  models,
  models_compfilt,
  models_gv0,
  models_gv50,
  models_gv50mr,
)

# useful when writing sbatch file
v2_models %>% print(n = 500)

write_rds(v2_models, "/home/xxx/rep2si/v2_models.rds")
