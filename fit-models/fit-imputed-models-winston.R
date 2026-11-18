library(tidyverse)
library(brms)
library(ordbetareg)

set.seed(4818)

# Use cmdstanr -- make sure this is installed on slurm
options(brms.backend = "cmdstanr")

# passed from SBATCH file
models_file <- Sys.getenv("MODELS_FILE")
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

if(Sys.info()["nodename"] == "clustername") {
  allocs_imput <- readRDS("/home/xxx/rep2si/allocs_imput.rds")
  models_file <- paste0("/home/xxx/rep2si/", models_file, ".rds")
  all_models <- readRDS(models_file)
} else {
  # Read from appropriate path when testing locally
  sharepoint_path <- Sys.getenv("REP2SI_DIR")
  allocs_imput <- readRDS(file.path(sharepoint_path, "rep2si-data/global-data/allocs_imput.rds"))
  models_file <- paste0(sharepoint_path, "/rep2si-analysis/global-analysis/alloc_models/", models_file, ".rds")
  all_models <- readRDS(models_file)
}

d <- allocs_imput[[task_id]]

# Free up memory
rm(allocs_imput)
gc()

# Read in model df

# Curent task
model_id <- as.integer(Sys.getenv("MODEL_ID", "0"))

if(is.na(model_id) | model_id == 0) stop("Model ID env var not read")

current_site <- all_models$site[model_id]
current_model_name <- all_models$model_name[model_id]
current_formula <- all_models$formula[model_id][[1]]
current_phireg <- all_models$phireg[model_id]
current_comprehension <- all_models$comprehension[model_id]
current_base_dir <- all_models$base_dir[model_id]
current_outcome <- all_models$outcome[model_id]

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

dir_path <- if (current_site == "total") {
  paste0(dir_path, "/", current_model_name, "_imput")
} else {
  paste0(dir_path, "/", current_model_name,"_", current_site, "_imput")
}

if(!dir.exists(dir_path)) dir.create(dir_path)

full_path <- if (current_site == "total") {
  file.path(dir_path, paste0(current_model_name, "_imp_", task_id))
} else {
  file.path(dir_path, paste0(current_model_name,"_", current_site, "_imp_", task_id))
}

if (current_comprehension != "all") {
  full_path <- paste0(full_path, "_comp_", current_comprehension)
}

print(paste0("Will save model as ", full_path))

if (current_outcome == "percGiven") {
  m <- ordbetareg(
    formula = current_formula,
    data = d,
    iter = 600,
    warmup = 500,
    true_bounds = c(0, 1),
    # Slightly more informative priors, to help regularise
    dirichlet_prior = c(2, 2, 1), # expect 0s and middle parts, but fewer 100
    intercept_prior_mean = 0,
    intercept_prior_SD = 1.5,       # Separate regularization for the intercept
    coef_prior_mean = 0,
    coef_prior_SD = 1,              # Regularizes all beta coefficients
    phi_reg = current_phireg,
    backend = "cmdstanr",
    extra_prior = c(
      set_prior("exponential(1)", class = "sd"),
      set_prior("lkj(2)", class = "cor"),
      # Expect some precision (exp(2) = 7.4), we know people wont be uniform in allocations
      set_prior("normal(2,1)", class = "Intercept", dpar = "phi"),
      # Regularise the variance in precision, penalise complexity
      set_prior("exponential(1)", class = "sd", dpar = "phi")
    ),
    control = list(adapt_delta = 0.95),
    chains = 4,
    # threads = threading(2),
    cores = 4,
    file = full_path
  )
} else if (current_outcome %in% c("give0", "give50", "give50ormore")) {
  # Run logistic model
  m <- brm(
    formula = current_formula,
    data = d,
    iter = 600,
    warmup = 500,
    backend = "cmdstanr",
    family = bernoulli(),
    prior = c(
      set_prior("normal(0, 1.5)", class = "Intercept"),
      set_prior("normal(0, 1)", class = "b"),
      set_prior("exponential(1)", class = "sd"),
      set_prior("lkj(2)", class = "cor")
    ),
    control = list(adapt_delta = 0.95),
    chains = 4,
    # threads = threading(2),
    cores = 4,
    file = full_path
  )
}
