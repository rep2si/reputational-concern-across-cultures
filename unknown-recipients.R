library(tidyverse)
library(ordbetareg)
library(brms)
library(marginaleffects)
library(kableExtra)
source("./rep2si-common/source_all_submodule.R")

sharepoint_path <- Sys.getenv("REP2SI_DIR")
plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

d <- readRDS(file.path(sharepoint_path, "rep2si-data", "global-data", "allocs_anonalter.rds")) %>%
  # To unify with what I've been doing elsewhere
  mutate(percGiven = percGiven / 100)

# Boxplot, bof
d %>%
  ggplot(aes(y = percGiven, x = factor(Condition))) +
  geom_boxplot() +
  facet_grid(cols = vars(site))

# Paired t-tests. Ultra-clear.
d_pair <- d %>%
  select(site, ID, Condition, percGiven) %>%
  pivot_wider(values_from = percGiven, names_from = Condition)

t.test(d_pair$revealed[d$site == "np"], d_pair$anonymous[d$site == "np"], paired = TRUE)
t.test(d_pair$revealed[d$site == "up"], d_pair$anonymous[d$site == "up"], paired = TRUE)

# Though normality assumption not met, massive 0 inflation
diffs <- d_pair$anonymous[d$site == "np"] - d_pair$revealed[d$site == "np"]
shapiro.test(diffs)
qqnorm(diffs);qqline(diffs)
hist(diffs)

# Model
f <- bf(
  percGiven ~
    Condition * mo(ego_cent_decile, id = "egomo") * site +
      (1 |ID),
  phi ~ site
)

m <- ordbetareg(
  formula = f,
  data = d,
  phi_reg = "both"
)


f2 <- bf(
  percGiven ~
    Condition *  site +
      (1 |ID),
  phi ~ site
)

m2 <- ordbetareg(
  formula = f2,
  data = d,
  phi_reg = "both"
)

fixef(m)

avg_comparisons(m, variables = "Condition", by = "site")

avg_comparisons(m, variables = "ego_cent_decile", by = "site")

cm <- c(
  "b_Intercept" = "Intercept",
  "b_Conditionrevealed" = "Allocation revealed",
  "b_siteup" = "Site UP",
  "b_Conditionrevealed:siteup" = "Revealed x Site UP",
  "bsp_moego_cent_decileidEQegomo" = "Donor Centrality",
  "bsp_moego_cent_decileidEQegomo:Conditionrevealed" = "Revealed x Donor Centrality",
  "bsp_moego_cent_decileidEQegomo:siteup" = "Donor Centrality x Site UP",
  "bsp_moego_cent_decileidEQegomo:Conditionrevealed:siteup" = "Revealed x Donor Centrality x Site UP",
  "b_phi_Intercept" = "Phi",
  "b_phi_siteup" = "Site UP (Phi)",
  "sd_ID__Intercept" = "SD Donor Intercept"
)

# Use brms variable names (b_*, sd_<group>__*) as in coef_table.R
make_coef_tab <- function(m) {
  s <- summary(m)
  bind_rows(
    rownames_to_column(s$fixed, "Parameter") |>
      # Monotonic effects are bsp_* in brms, everything else b_*
      mutate(Parameter = paste0(ifelse(grepl("^mo", Parameter), "bsp_", "b_"), Parameter)),
    imap(s$random, \(t, group) rownames_to_column(t, "Parameter") |>
      mutate(Parameter = paste0("sd_", group, "__", gsub("^sd\\((.*)\\)$", "\\1", Parameter)))) |>
      list_rbind()
  ) |>
    mutate(Parameter = coalesce(cm[Parameter], Parameter))
}

coef_tabs <- list(
  "Without ego centrality parameter" = make_coef_tab(m2),
  "With ego centrality parameter" = make_coef_tab(m)
)
coef_tabs

coef_tab <- bind_rows(coef_tabs)
# Model names as subtitle rows (pack_rows)
model_rows <- map_int(coef_tabs, nrow)

coef_tab |>
  mutate(
    Parameter = gsub("([_%&#$])", "\\\\\\1", Parameter),
    Est = sprintf("%.2f", Estimate),
    SE = sprintf("%.2f", Est.Error),
    CI = sprintf("[%.2f, %.2f]", `l-95% CI`, `u-95% CI`),
    Rhat = sprintf("%.2f", Rhat),
    Bulk = sprintf("%.0f", Bulk_ESS),
    Tail = sprintf("%.0f", Tail_ESS),
    .keep = "none"
  ) |>
  # .keep = "none" leaves modified columns (Rhat) in their original position
  select(Parameter, Est, SE, CI, Rhat, Bulk, Tail) |>
  kable(
    format = "latex",
    booktabs = TRUE,
    align = c("l", rep("r", 6)),
    col.names = c("Parameter", "Est.", "SE", "95\\% CI", "Rhat", "Bulk ESS", "Tail ESS"),
    escape = FALSE,
    caption = "Estimated coefficients and sampler diagnostics. Models fit to allocations made to unknown recipient only.",
    label = "coef-unknown-recipients"
  ) |>
  pack_rows(index = model_rows) |>
  as.character() |>
  writeLines(file.path(plots_dir, "coef_unknown_recipients.tex"))
