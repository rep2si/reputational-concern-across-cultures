library(tidyverse)
library(brms)
library(ordbetareg)

set.seed(4818)

# Use cmdstanr -- make sure this is installed on slurm
options(brms.backend = "cmdstanr")

if(Sys.info()["nodename"] == "clustername") {
  d <- readRDS("/home/xxx/rep2si/allocs.rds")
} else {
  # Read from appropriate path when testing locally
  sharepoint_path <- Sys.getenv("REP2SI_DIR")
  d <- readRDS(file.path(sharepoint_path, "rep2si-data/global-data/allocs.rds"))
}

# We have reliable 4-level measure in UP and NP only
d_sub <- d %>%
  filter(site %in% c("np", "up")) %>%
  mutate(site = fct_drop(site))

m <- ordbetareg(
  formula = bf(
    percGiven ~
        Condition * anon_first +
        Condition * ego_comp * site +
        (Condition | ID) +
        (1 | AID) +
        (1 | RID),
      phi ~ site
  ),
  data = d_sub,
  iter = 2000,
  true_bounds = c(0, 1),
  # Slightly more informative priors, to help regularise
  dirichlet_prior = c(2, 2, 1), # expect 0s and middle parts, but fewer 100
  intercept_prior_mean = 0,
  intercept_prior_SD = 1.5,       # Separate regularization for the intercept
  coef_prior_mean = 0,
  coef_prior_SD = 1,             # Regularizes all beta coefficients
  phi_reg = "both",
  backend = "cmdstanr",
  extra_prior = c(
    set_prior("exponential(1)", class = "sd"),
    set_prior("lkj(2)", class = "cor"),
    # Expect some precision (exp(2) = 7.4), we know people wont be uniform in allocations
    set_prior("normal(2,1)", class = "Intercept", dpar = "phi")
  ),
  # required to use moment_match for loo
  save_pars = save_pars(all = TRUE),
  # control = list(adapt_delta = 0.9),
  chains = 4,
  cores = 4,
  threads = threading(4),
  file = "/scratch/xxx/rep2si/results/comp/4_levels_componly"
)
