# Compares site-specific fixed effects from the pooled v3 cd model (random
# intercepts and slopes by site) with those from the site-by-site models fit
# with no pooling.

library(tidyverse)
library(brms)

sharepoint_path <- Sys.getenv("REP2SI_DIR")
model_dir <- file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/v3/")
plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

m_cd <- readRDS(file.path(model_dir, "cd.rds"))

sitewise_files <- list.files(file.path(model_dir, "sitewise"), pattern = "^cd_.*\\.rds$", full.names = TRUE)
m_sites <- map(sitewise_files, readRDS) %>%
  set_names(str_remove(basename(sitewise_files), "^cd_") %>% str_remove("\\.rds$"))

# Partial pooling: fixef + site ranef, posterior median and 95% CI
coef_pooled <- coef(m_cd, robust = TRUE)$site
pooled <- map(dimnames(coef_pooled)[[3]], function(term) {
  as.data.frame(coef_pooled[, , term]) %>%
    rownames_to_column("site") %>%
    mutate(term = term)
}) %>%
  bind_rows() %>%
  mutate(model = "partial pooling")

# No pooling: fixef of each site model
nopool <- imap(m_sites, function(m, site) {
  # Site-wise models have a scalar phi; log it to match the pooled model's
  # phi_Intercept (log link)
  log_phi <- posterior_summary(log(as.matrix(m, variable = "phi")), robust = TRUE)
  rownames(log_phi) <- "phi_Intercept"

  rbind(fixef(m, robust = TRUE), log_phi) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    mutate(site = site)
}) %>%
  bind_rows() %>%
  mutate(model = "no pooling")

# Global fixed effects from the pooled model, as a reference
global <- fixef(m_cd, robust = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("term")

# Terms present in only one of the two sources (e.g. factor level missing at a site)
only_pooled <- setdiff(unique(pooled$term), unique(nopool$term))
only_nopool <- setdiff(unique(nopool$term), unique(pooled$term))
if (length(only_pooled) > 0) warning("Terms only in pooled model: ", paste(only_pooled, collapse = ", "))
if (length(only_nopool) > 0) warning("Terms only in site-wise models: ", paste(only_nopool, collapse = ", "))

all_coefs <- bind_rows(pooled, nopool) %>%
  select(model, site, term, Estimate, Q2.5, Q97.5)

# Side-by-side table ------------------------------------------------------------

comparison <- all_coefs %>%
  mutate(model = if_else(model == "partial pooling", "pooled", "nopool")) %>%
  pivot_wider(
    id_cols = c(term, site),
    names_from = model,
    values_from = c(Estimate, Q2.5, Q97.5)
  ) %>%
  mutate(diff = Estimate_pooled - Estimate_nopool) %>%
  select(
    term, site,
    Estimate_pooled, Q2.5_pooled, Q97.5_pooled,
    Estimate_nopool, Q2.5_nopool, Q97.5_nopool,
    diff
  ) %>%
  arrange(term, site)

print(comparison, n = Inf)

# Coefficient labels (from coef_table.R, without the b_/bsp_ prefixes) ----------

cm <- c(
  "Intercept" = "Intercept",
  "Conditionrevealed" = "Allocation revealed",
  "moego_cent_decileidEQegomo" = "Donor Centrality",
  "moalter_cent_decileidEQaltermo" = "Recipient Centrality",
  "modist_ordidEQdistmo:Conditionrevealed" = "Revealed x Distance",
  "moego_cent_decileidEQegomo:Conditionrevealed" = "Revealed x Donor Centrality",
  "moalter_cent_decileidEQaltermo:Conditionrevealed" = "Revealed x Recipient Centrality",
  "ego_age_demean_div" = "Donor Age",
  "alter_age_demean_div" = "Recipient Age",
  "ego_genderMale" = "Donor Male",
  "alter_genderMale" = "Recipient Male",
  "anon_first" = "Anonymous First",
  "ego_demean_logwealth" = "Donor Wealth",
  "alter_demean_logwealth" = "Recipient Wealth",
  "Conditionrevealed:ego_age_demean_div" = "Revealed x Donor Age",
  "Conditionrevealed:alter_age_demean_div" = "Revealed x Recipient Age",
  "Conditionrevealed:ego_genderMale" = "Revealed x Donor Male",
  "Conditionrevealed:alter_genderMale" = "Revealed x Recipient Male",
  "Conditionrevealed:anon_first" = "Revealed x Anonymous First",
  "Conditionrevealed:ego_demean_logwealth" = "Revealed x Donor wealth",
  "Conditionrevealed:alter_demean_logwealth" = "Revealed x Recipient Wealth",
  "modist_ordidEQdistmo" = "Distance",
  "phi_Intercept" = "Phi"
)

# Keep unmatched terms under their raw name, and order facets as in the key
rename_terms <- function(df) {
  df %>%
    mutate(
      term = coalesce(cm[term], term),
      term = factor(term, levels = unique(c(cm, term)))
    )
}

plot_coefs <- rename_terms(all_coefs)
plot_global <- rename_terms(global)

# Plot by site ------------------------------------------------------------------

# Intercept and phi are on a much larger scale than the slopes, so they would
# squash everything else within a site panel
by_site_coefs <- plot_coefs %>%
  bind_rows(tidyr::crossing(plot_global, site = unique(all_coefs$site), model = "fixed effect")) %>%
  filter(!term %in% c("Intercept", "Phi")) %>%
  # First level is dodged to the bottom, so fixed effect ends up on top
  mutate(model = factor(model, levels = c("no pooling", "partial pooling", "fixed effect")))

by_site_coefs %>%
  ggplot(aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = fct_rev(fct_drop(term)), colour = model)) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.2) +
  facet_wrap(vars(str_to_upper(site))) +
  scale_colour_manual(
    values = c("fixed effect" = "grey50", "no pooling" = "#F8766D", "partial pooling" = "#00BFC4"),
    breaks = c("fixed effect", "partial pooling", "no pooling")
  ) +
  labs(
    x = "Posterior median and 95% CI",
    y = NULL,
    colour = NULL
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, "pooled-vs-sitewise-coefs-by-site.pdf"), height = 12, width = 14)
