library(tidyverse)
library(brms)
library(ordbetareg)

set.seed(4818)

# Use cmdstanr -- make sure this is installed on slurm
options(brms.backend = "cmdstanr")

# passed from SBATCH file
models_file <- Sys.getenv("MODELS_FILE")
use_subset <- Sys.getenv("USE_SUBSET")
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

if(Sys.info()["nodename"] == "clustername") {
  d <- readRDS("/home/xxx/rep2si/allocs.rds")
  # d_subset <- readRDS("/home/xxx/rep2si/allocs_subsampled.rds")
  models_file <- paste0("/home/xxx/rep2si/", models_file, ".rds")
  all_models <- readRDS(models_file)
} else {
  # Read from appropriate path when testing locally
  sharepoint_path <- Sys.getenv("REP2SI_DIR")
  d <- readRDS(file.path(sharepoint_path, "rep2si-data/global-data/allocs.rds"))
  # d_subset <- readRDS(file.path(sharepoint_path, "rep2si-data/global-data/allocs_subsampled.rds"))
  models_file <- paste0(sharepoint_path, "/rep2si-analysis/global-analysis/alloc_models/", models_file, ".rds")
  all_models <- readRDS(models_file)
}

if(use_subset == 1) {
  d <- d_subset
}

# Curent task
current_site <- all_models$site[task_id]
current_model_name <- all_models$model_name[task_id]
current_formula <- all_models$formula[task_id][[1]]
# Set but removed below
current_phireg <- all_models$phireg[task_id]
current_comprehension <- all_models$comprehension[task_id]
current_base_dir <- all_models$base_dir[task_id]
current_outcome <- all_models$outcome[task_id]

# Filter to site if required
if(current_site %in% c("te", "wl", "wh", "tz", "az", "tm", "np", "up")) {
  d <- d %>% filter(site == current_site)
}

# Special site filter for interacted comprehension models
if (current_site == "comp2lev"){
  d <- d %>% filter(site %in% c("az", "tm", "up", "np"))
}

if (current_site == "comp4lev"){
  d <- d %>% filter(site %in% c("np", "up"))
}

# Filter comprehension if required
if(!(current_comprehension %in% c("all", "lenient", "intermediate", "strict"))) {
  stop("Unknown comprehension level")
} else {
  d <- switch(
    current_comprehension,
    "strict" = d %>% filter(ego_comp_strict == 1),
    "intermediate" = d %>% filter(ego_comp_intermediate == 1),
    "lenient" = d %>% filter(ego_comp_lenient == 1),
    "all" = d
  )
}

# Ensure base_dir exists
dir_path <- paste0(
  "/scratch/xxx/rep2si/results/",
  current_base_dir
)

if(!dir.exists(dir_path)) dir.create(dir_path)

full_path <- if (current_site == "total") {
  file.path(dir_path, current_model_name)
} else {
  file.path(dir_path, paste0(current_model_name,"_", current_site))
}

if (current_comprehension != "all") {
  full_path <- paste0(full_path, "_comp_", current_comprehension)
}

print(paste0("Will save model as ", full_path))

if (current_outcome == "percGiven") {
  m <- ordbetareg(
    formula = current_formula,
    data = d,
    iter = 2000,
    true_bounds = c(0, 1),
    dirichlet_prior = c(2, 2, 1), # expect 0s and middle parts, but fewer 100
    intercept_prior_mean = 0,
    intercept_prior_SD = 1.5,       # Separate regularization for the intercept
    coef_prior_mean = 0,
    coef_prior_SD = 1,              # Regularizes all beta coefficients
    backend = "cmdstanr",
    extra_prior = c(
      set_prior("exponential(1)", class = "sd"),
      set_prior("lkj(2)", class = "cor")
    ),
    # required to use moment_match for loo
    # save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.95),
    chains = 4,
    threads = threading(2),
    cores = 4,
    file = full_path
  )
} else if (current_outcome %in% c("give0", "give50", "give50ormore")) {
  # Run logistic model
  m <- brm(
    formula = current_formula,
    data = d,
    iter = 2000,
    backend = "cmdstanr",
    family = bernoulli(),
    prior = c(
      set_prior("normal(0, 1.5)", class = "Intercept"),
      set_prior("normal(0, 1)", class = "b"),
      set_prior("exponential(1)", class = "sd"),
      set_prior("lkj(2)", class = "cor")
    ),
    # save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.95),
    chains = 4,
    # threads = threading(2),
    cores = 4,
    file = full_path
  )
}
