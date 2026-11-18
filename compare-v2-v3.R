# Compares coefficients from v2 and v3 models, where the only difference for the
# latter is that they include a random effect within ID for dist * cond

library(tidyverse)
library(bayestestR)
library(brms)
library(posterior)

sharepoint_path <- Sys.getenv("REP2SI_DIR")
model_dir <- file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/")
plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

m_cd2 <- readRDS(file.path(model_dir, "v2", "cd.rds"))
m_cd2_gv0 <- readRDS(file.path(model_dir, "v2_gv0", "cd.rds"))
m_cd2_gv50mr <- readRDS(file.path(model_dir, "v2_gv50mr", "cd.rds"))

# Version with ID random slopes for dist
m_cd3 <- readRDS(file.path(model_dir, "v3", "cd.rds"))
m_cd3_gv0 <- readRDS(file.path(model_dir, "v3_gv0", "cd.rds"))
m_cd3_gv50mr <- readRDS(file.path(model_dir, "v3_gv50mr", "cd.rds"))

# Plot difference in fixed effects
get_coefs <- function(m) {
  fixef(get(m)) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    mutate(
      model = case_when(
        str_detect(m, "_gv50mr") ~ "M2",
        str_detect(m, "_gv0")   ~ "M3",
        TRUE                    ~ "M1"
      ),
      with_ranslope = factor(
        if_else(str_detect(m, "cd3"), "yes", "no"),
        levels = c("yes", "no"),
        labels = c(
          "Model allowing Condition x Distance to vary by Participant",
          "Model not allowing Condition x Distance to vary by participant"
        )
      )
    )
}

all_coefs <- lapply(
  c("m_cd2", "m_cd2_gv0", "m_cd2_gv50mr", "m_cd3", "m_cd3_gv0", "m_cd3_gv50mr"),
  get_coefs
)

all_coefs <- bind_rows(all_coefs)

all_coefs %>%
  filter(term %in% c("modist_ordidEQdistmo:Conditionrevealed", "modist_ordidEQdistmo")) %>%
  mutate(term = case_when(
    term == "modist_ordidEQdistmo"                    ~ "Distance",
    term == "modist_ordidEQdistmo:Conditionrevealed"  ~ "Revealed x Distance"
  )) %>%
  ggplot(aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, colour = with_ranslope, y = term)) +
  geom_pointrange(position = position_dodge(width = 0.2)) +
  facet_grid(rows = vars(model)) +
  geom_vline(xintercept = 0, linetype = 3) +
  scale_colour_brewer(type = "qual", palette = 2, name = NULL) +
  cowplot::theme_cowplot() +
  labs(y = NULL) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical"
  )

ggsave(file.path(plots_dir, "v2-v3-comp.pdf"), width = 10, height = 7)
