library(tidyverse)
library(brms)
library(patchwork)
source("./rep2si-common/source_all_submodule.R")

plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

set_theme(theme_minimal())

sharepoint_path <- Sys.getenv("REP2SI_DIR")
d <- readRDS(file.path(sharepoint_path, "rep2si-data/global-data/allocs.rds"))
model_dir <- file.path(sharepoint_path, "rep2si-analysis/global-analysis/alloc_models/v2/")

m_cd <- readRDS(file.path(model_dir, "cd.rds"))

# Get only the draws for a model -----------------------------------------------

# Caterpillar plot of giver effects
ran_ef <- ranef(m_cd, robust = TRUE)

# Get effect of condition to compare random effects
cond_eff <- fixef(m_cd, robust = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("param") %>%
  filter(param == "Conditionrevealed")

p1 <- as.data.frame(ran_ef$ID) %>%
  tibble() %>%
  arrange(Estimate.Intercept) %>%
  mutate(group = seq_len(nrow(ran_ef$ID))) %>%
  ggplot(aes(x = group, y = Estimate.Intercept, ymin = Q2.5.Intercept, ymax = Q97.5.Intercept)) +
  geom_pointrange(alpha = 0.3, size = 0.1) +
  geom_hline(yintercept = 0, colour = "blue") +
  scale_y_continuous(limits = c(-6, 6)) +
  labs(
    title = "Donor intercepts",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
  ) +
  geom_pointrange(x = - 0.5, y = cond_eff$Estimate, ymin = cond_eff$Q2.5, ymax = cond_eff$Q97.5, colour = "red", size = 0.1)

p2 <- as.data.frame(ran_ef$AID) %>%
  tibble() %>%
  arrange(Estimate.Intercept) %>%
  mutate(group = seq_len(nrow(ran_ef$AID))) %>%
  ggplot(aes(x = group, y = Estimate.Intercept, ymin = Q2.5.Intercept, ymax = Q97.5.Intercept)) +
  geom_pointrange(alpha = 0.3, size = 0.1)  +
  geom_hline(yintercept = 0, colour = "blue") +
  scale_y_continuous(limits = c(-6, 6)) +
  labs(
    title = "Recipient intercepts",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
  ) +
  geom_pointrange(x = - 0.5, y = cond_eff$Estimate, ymin = cond_eff$Q2.5, ymax = cond_eff$Q97.5, colour = "red", size = 0.1)

p3 <- as.data.frame(ran_ef$RID) %>%
  tibble() %>%
  arrange(Estimate.Intercept) %>%
  mutate(group = seq_len(nrow(ran_ef$RID))) %>%
  ggplot(aes(x = group, y = Estimate.Intercept, ymin = Q2.5.Intercept, ymax = Q97.5.Intercept)) +
  geom_pointrange(alpha = 0.9, size = 0.1) +
  geom_hline(yintercept = 0, colour = "blue") +
  scale_y_continuous(limits = c(-6, 6)) +
  labs(
    title = "Research Assistant intercepts",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
  ) +
  geom_pointrange(x = 0, y = cond_eff$Estimate, ymin = cond_eff$Q2.5, ymax = cond_eff$Q97.5, colour = "red", size = 0.1)


p4 <- as.data.frame(ran_ef$site) %>%
  rownames_to_column("site") %>%
  tibble() %>%
  arrange(Estimate.Intercept) %>%
  mutate(group = seq_len(nrow(ran_ef$site))) %>%
  ggplot(aes(x = group, y = Estimate.Intercept, ymin = Q2.5.Intercept, ymax = Q97.5.Intercept)) +
  geom_pointrange(alpha = 0.9, size = 0.1) +
  geom_hline(yintercept = 0, colour = "blue") +
  geom_text(aes(label = site), hjust = 0, nudge_x = 0.1) +
  # scale_y_continuous(limits = c(-6, 6)) +
  labs(
    title = "Site intercepts",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
  ) +
  scale_y_continuous(limits = c(-6, 6)) +
  geom_pointrange(x = 0.8, y = cond_eff$Estimate, ymin = cond_eff$Q2.5, ymax = cond_eff$Q97.5, colour = "red", size = 0.1)

p1 + p2 + p3 + p4

ggsave(file.path(plots_dir, "ranef.pdf"), height = 10, width = 10)

# Random lsopes for condition by ID
p5 <- as.data.frame(ran_ef$ID) %>%
  tibble() %>%
  arrange(Estimate.Conditionrevealed) %>%
  mutate(group = seq_len(nrow(ran_ef$ID))) %>%
  ggplot(aes(x = group, y = Estimate.Conditionrevealed, ymin = Q2.5.Conditionrevealed, ymax = Q97.5.Conditionrevealed)) +
  geom_pointrange(alpha = 0.3, size = 0.1, colour = "darkgreen") +
  geom_hline(yintercept = 0, colour = "blue") +
  scale_y_continuous(limits = c(-6, 6)) +
  labs(
    title = "ID Random slopes for Condition",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
  ) +
  geom_pointrange(x = - 0.5, y = cond_eff$Estimate, ymin = cond_eff$Q2.5, ymax = cond_eff$Q97.5, colour = "red", size = 0.1)

p6 <- as.data.frame(ran_ef$site) %>%
  rownames_to_column("site") %>%
  tibble() %>%
  arrange(Estimate.Conditionrevealed) %>%
  mutate(group = seq_len(nrow(ran_ef$site))) %>%
  ggplot(aes(x = group, y = Estimate.Conditionrevealed, ymin = Q2.5.Conditionrevealed, ymax = Q97.5.Conditionrevealed)) +
  geom_pointrange(alpha = 0.3, size = 0.1, colour = "darkgreen") +
  geom_hline(yintercept = 0, colour = "blue") +
  geom_text(aes(label = site), hjust = 0, nudge_x = 0.1) +
  scale_y_continuous(limits = c(-6, 6)) +
  labs(
    title = "Site random slopes for Condition",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
  ) +
  geom_pointrange(x = - 0.5, y = cond_eff$Estimate, ymin = cond_eff$Q2.5, ymax = cond_eff$Q97.5, colour = "red", size = 0.1)

p5 + p6

ggsave(file.path(plots_dir, "ranef-slopes.pdf"), height = 5, width = 10)
