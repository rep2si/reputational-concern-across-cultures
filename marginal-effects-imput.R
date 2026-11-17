library(tidyverse)
library(Matrix)
library(tidybayes)
library(brms)
library(posterior)
library(marginaleffects)
library(patchwork)
library(ggrepel)
library(patchwork)
library(cowplot)
library(ggdist)

options(pillar.sigfig = 3)
source("./rep2si-common/source_all_submodule.R")

set_theme(theme_cowplot() + theme(plot.title.position = "plot", plot.subtitle.position = "plot"))

sharepoint_path <- Sys.getenv("REP2SI_DIR")
plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")
site_order <- c("te", "tm", "az", "up", "np", "wh", "tz", "wl", "total")
site_order2 <- c("total", "np", "up", "te", "tm", "wl", "az", "tz", "wh")
cond_cols <- c(anonymous = "#0072B2", revealed = "#D55E00")

# Retrieve names of all models, saved as rds, in a directory
get_model_names <- function (model_dir) {
  files <- list.files(model_dir, full.names = TRUE)
  # fix order of filenames
  nums <- as.numeric(gsub("\\D", "", basename(files)))
  model_names <- files[order(nums)]
  return(model_names)
}

# Generates the various counterfactual grids used to calculated predicted values
gen_ctf <- function(model_dir, ctf_type, extreme_cent_dec_only = FALSE) {

  model_names <- get_model_names(model_dir)
  d_m1 <- readRDS(model_names[1])$data

  # extract features we want in the preds df
  d_basic <- d_m1 %>%
    select(site, ID, AID, Condition, dist_ord, ego_cent_decile, alter_cent_decile, anon_first, ID_SUID, AID_SUID, dyad_repeated, dyad_id_reponly)

  # sanity check
  any_dup <- any(duplicated(d_basic))
  if(any_dup) stop("Duplicated rows in d_basic")

  # Add weights for poststratification
  d_basic <- d_basic %>%
    # pi_aid = probability a given recipient was shown to a participant,
    # estimated as the share of participants in that site x condition who saw them.
    group_by(site, Condition) %>%
    mutate(n_id_site_cond = n_distinct(ID)) %>%
    group_by(site, Condition, AID) %>%
    mutate(n_id_aid = n_distinct(ID)) %>%
    mutate(pi_aid = n_id_aid / n_id_site_cond) %>%
    ungroup() %>%
    # Hajek weights
    # 1/pi_aid undoes the oversampling: a recipient shown twice as often
    # counts half as much per appearance.
    # Dividing by the participant's own sum of 1/pi_aid (the Hajek normalization)
    # does two things: (a) cancels draw-to-draw variation in which recipients a
    # participant happened to see, and (b) makes each participant's weights sum
    # to 1 in each condition, so participants count equally regardless of how
    # many allocations they made.
    group_by(ID, Condition) %>%
    mutate(w = (1 / pi_aid) / sum(1 / pi_aid)) %>%
    ungroup() %>%
    # Site rescale
    group_by(site) %>%
    mutate(n_id_site = n_distinct(ID)) %>%
    ungroup() %>%
    mutate(w_s = w / n_id_site) %>%
    select(-n_id_site_cond, -n_id_aid, -n_id_site)

  # Build counterfactual dataframe for variable(s) whose effect we are isolating
  d_ctf <- switch(
    ctf_type,
    "cond" = bind_rows(
      d_basic %>% mutate(Condition = "anonymous"),
      d_basic %>% mutate(Condition = "revealed"),
    ),
    "dist" = expand_grid(
      d_basic %>% select(-Condition, -dist_ord),
      Condition = c("anonymous", "revealed"),
      dist_ord = levels(d_m1$dist_ord)
    ),
    "alter_cent" = if (extreme_cent_dec_only) {
      expand_grid(
        d_basic %>% select(-Condition, -alter_cent_decile),
        Condition = c("anonymous", "revealed"),
        alter_cent_decile = c("[0,0.1]","(0.9,1]")
      )
    } else {
      expand_grid(
        d_basic %>% select(-Condition, -alter_cent_decile),
        Condition = c("anonymous", "revealed"),
        alter_cent_decile = c("[0,0.1]", "(0.1,0.2]", "(0.3,0.4]", "(0.4,0.5]", "(0.5,0.6]", "(0.6,0.7]", "(0.7,0.8]", "(0.8,0.9]","(0.9,1]")
      )
    },
    "ego_cent" = if (extreme_cent_dec_only) {
      expand_grid(
        d_basic %>% select(-Condition, -ego_cent_decile),
        Condition = c("anonymous", "revealed"),
        ego_cent_decile = c("[0,0.1]","(0.9,1]")
      )
    } else{
      expand_grid(
        d_basic %>% select(-Condition, -ego_cent_decile),
        Condition = c("anonymous", "revealed"),
        ego_cent_decile = c("[0,0.1]", "(0.1,0.2]", "(0.3,0.4]", "(0.4,0.5]", "(0.5,0.6]", "(0.6,0.7]", "(0.7,0.8]", "(0.8,0.9]","(0.9,1]")
      )
    }
  )
  return(d_ctf)
}

# Calculate counterfactual predictions in a memory efficient way (wihtout
# repeating all covariates for each posterior draw).
gen_my_epreds <- function(model_dir, ctf_type, n_dr, extreme_cent_dec_only = FALSE) {
  model_names <- get_model_names(model_dir)
  all_epreds <- list()
  for(i in seq_len(length(model_names))){
    message(paste("Calculating epreds for imputation", i))
    m <- readRDS(model_names[i])
    # Use the imputed dataset for the current model.
    d <- m$data
    # Prepare counterfactual, same as in the resutls df
    d_ctf_tmp <- switch(
      ctf_type,
      "cond" = bind_rows(
        d %>% mutate(Condition = "anonymous"),
        d %>% mutate(Condition = "revealed"),
      ),
      "dist" = expand_grid(
        d %>% select(-Condition, -dist_ord),
        Condition = c("anonymous", "revealed"),
        dist_ord = c("1", "2", "3", "4+")
      ),
      "alter_cent" = if (extreme_cent_dec_only) {
        expand_grid(
          d %>% select(-Condition, -alter_cent_decile),
          Condition = c("anonymous", "revealed"),
          alter_cent_decile = c("[0,0.1]", "(0.9,1]")
        )
      } else {
        expand_grid(
          d %>% select(-Condition, -alter_cent_decile),
          Condition = c("anonymous", "revealed"),
          alter_cent_decile = c("[0,0.1]", "(0.1,0.2]", "(0.3,0.4]", "(0.4,0.5]", "(0.5,0.6]", "(0.6,0.7]", "(0.7,0.8]", "(0.8,0.9]","(0.9,1]")
        )
      },
      "ego_cent" = if (extreme_cent_dec_only) {
        expand_grid(
          d %>% select(-Condition, -ego_cent_decile),
          Condition = c("anonymous", "revealed"),
          ego_cent_decile = c("[0,0.1]","(0.9,1]")
        )
      } else {
        expand_grid(
          d %>% select(-Condition, -ego_cent_decile),
          Condition = c("anonymous", "revealed"),
          ego_cent_decile = c("[0,0.1]", "(0.1,0.2]", "(0.3,0.4]", "(0.4,0.5]", "(0.5,0.6]", "(0.6,0.7]", "(0.7,0.8]", "(0.8,0.9]","(0.9,1]")
        )
      }
    )
    epreds <- posterior_epred(m, newdata = d_ctf_tmp, ndraws = n_dr)
    all_epreds[[i]] <- epreds
  }

  epreds_full <- do.call(rbind, all_epreds)
  return(epreds_full)
}

# Code below computes model predictions and contrasts reported in paper, and
# creates corresponding figures.

# Just Condition ---------------------------------------------------------------
d_ctf <- gen_ctf(
  model_dir = file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v3", "cd_imput"),
  ctf_type = "cond"
)

d_ctf$.epred <- rvar(
  gen_my_epreds(
  model_dir = file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v3", "cd_imput"),
  ctf_type = "cond",
  n_dr = 200
  )
)

# Preds ctf
pred_all_sites <- d_ctf %>%
  group_by(Condition) %>%
  # mean pred per draw and condition
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  # summarise(pred_value = mean(.epred)) %>%
  median_qi(wmean) %>%
  mutate(site = "total")

pred_by_site <- d_ctf %>%
  group_by(site, Condition) %>%
  # mean pred per draw and condition
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  # summarise(pred_value = mean(.epred)) %>%
  median_qi(wmean)

preds <- bind_rows(pred_all_sites, pred_by_site)

preds <- preds %>%
  mutate(site = factor(site, levels = site_order))

# Fig S1
ggplot(preds, aes(x = wmean, xmin = .lower, xmax = .upper, y = site, colour = Condition)) +
  geom_pointrange(position = position_dodge(width = 0.4)) +
  labs(
    title = "Allocations by condition",
    x = "Percentage of total endowment given",
    y = NULL,
    shape = "Condition",
    ) +
  scale_colour_manual(values = cond_cols) +
  scale_x_continuous(labels = scales::percent) +
  guides(colour = "none")

# Preds by *observed* donor centrality---used for later figure
preds_egocent <- d_ctf %>%
  group_by(ego_cent_decile) %>%
  # mean pred per draw and condition
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  # summarise(pred_value = mean(.epred)) %>%
  median_qi(wmean) %>%
  mutate(site = "total")

p_obs_egocent <- preds_egocent %>%
  # mutate(site = factor(site, levels = eff_order)) %>%
  ggplot(aes(x = wmean, xmin = .lower, xmax = .upper, y = ego_cent_decile)) +
  geom_pointrange(position = position_dodge(width = 0.2)) +
  labs(
    title = "Predicted allocation by observed donor centrality",
    x = "Percentage of endowment given",
    y = NULL
  ) +
  scale_x_continuous(labels = scales::percent, limits = c(0.22, 0.33)) +
  theme(legend.position = "none") +
  scale_y_discrete(labels = get_site_labels())

# Effect of condition ----------------------------------------------------------

# Useful to compare effect sizes with and without site weights
eff_all_nositeweight <- d_ctf %>%
  group_by(Condition) %>%
  summarise(wmean = rvar_sum(.epred * w) / sum(w)) %>%
  pivot_wider(names_from=Condition, values_from = wmean) %>%
  mutate(cond_eff = revealed - anonymous) %>%
  median_qi(cond_eff) %>%
  mutate(site = "total")

eff_all <- d_ctf %>%
  group_by(Condition) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from=Condition, values_from = wmean) %>%
  mutate(cond_eff = revealed - anonymous) %>%
  median_qi(cond_eff) %>%
  mutate(site = "total")

eff_site <- d_ctf %>%
  group_by(Condition, site) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from=Condition, values_from = wmean) %>%
  mutate(cond_eff = revealed - anonymous) %>%
  median_qi(cond_eff) %>%
  arrange(cond_eff)

eff <- bind_rows(eff_all, eff_site) %>%
  mutate(site = factor(site, levels = site_order))

# Overall absolute effect sizes reported in article
eff

# Effect sizes relative to amount given in anonymous condition
eff_rel_all <- d_ctf %>%
  group_by(Condition) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from=Condition, values_from = wmean) %>%
  mutate(cond_eff_rel = (revealed / anonymous) - 1) %>%
  median_qi(cond_eff_rel) %>%
  mutate(site = "total")

eff_rel_site <- d_ctf %>%
  group_by(Condition, site) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from=Condition, values_from = wmean) %>%
  mutate(cond_eff_rel = (revealed / anonymous) - 1) %>%
  median_qi(cond_eff_rel) %>%
  arrange(cond_eff_rel)

eff_rel <- bind_rows(eff_rel_all, eff_rel_site) %>%
  mutate(site = factor(site, levels = site_order))

# Relative effect sizes reported in article
eff_rel

eff_order <- eff %>%
  arrange(desc(cond_eff)) %>%
  pull(site)

eff_order <- rev(c("total", setdiff(eff_order, "total")))

# Fig S2
p2 <- eff %>%
  mutate(site = factor(site, levels = eff_order)) %>%
  ggplot(aes(x = cond_eff, xmin = .lower, xmax = .upper, y = site)) +
  geom_pointrange(position = position_dodge(width = 0.2)) +
  labs(
    title = "Effect of reputational exposure",
    x = "Additonal percentage of endowment given when identity is revealed",
    y = NULL
  ) +
  scale_x_continuous(labels = scales::percent) +
  geom_vline(xintercept = 0, linetype = 3) + scale_color_manual(values = get_site_colours()) +
  theme(legend.position = "none") +
  scale_y_discrete(labels = get_site_labels())
p2

# LS-style plot -- not used in article. This is equivalent to the LS-style plot
# using raw descriptives, but with model predictions instead. Done as a sanity
# check, results are very similar
d_ctf %>%
  group_by(Condition, site) %>%
  # mean pred per draw and condition
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  # summarise(pred_value = mean(.epred)) %>%
  median_qi(wmean) %>%
  pivot_wider(
    id_cols = site,
    names_from = Condition,
    values_from = c(wmean, .lower, .upper)
  )  %>%
  ggplot(
    aes(
      x = wmean_anonymous, xmin = .lower_anonymous, xmax = .upper_anonymous,
      y = wmean_revealed, ymin = .lower_revealed, ymax = .upper_revealed,
      colour = site
    )
  ) +
  geom_errorbar(orientation = "x", width = 0, linewidth = 1) +
  geom_errorbar(orientation = "y", width = 0, linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 3) +
  scale_colour_manual(values = get_site_colours()) +
  geom_text_repel(aes(label = str_to_upper(site)), box.padding = 0.7, segment.colour = NA) +
  scale_x_continuous(limits = c(0.04, 0.45), labels = scales::percent) +
  scale_y_continuous(limits = c(0.04, 0.45), labels = scales::percent) +
  labs(
    x = "Anonymous (Mean model prediction and 95% Cred Int)",
    y = "Revealed (Mean model prediction and 95% Cred Int)",
    title = "Weighted model predictions"
  ) +
  theme(legend.position = "none")

# Need lots of memory for next steps!
rm(d_ctf)
gc()

# Distance ---------------------------------------------------------------------
# Generate counterfactual predictions for all combinations of distance and conditon.

d_ctf_dist <- gen_ctf(
  model_dir = file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v3", "cd_imput"),
  ctf_type = "dist"
)

d_ctf_dist$.epred <- rvar(
  gen_my_epreds(
  model_dir = file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v3", "cd_imput"),
  ctf_type = "dist",
  n_dr = 200
  )
)

preds_dist_all <- d_ctf_dist %>%
  group_by(Condition, dist_ord) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  mutate(site = "total")

preds_dist_site <- d_ctf_dist %>%
  group_by(site, Condition, dist_ord) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s))

preds_dist <- bind_rows(preds_dist_all, preds_dist_site)

preds_dist <- preds_dist %>%
  mutate(site = factor(site, levels = site_order))

pred_dist_sum <- preds_dist %>% median_qi(wmean, .width = 0.95)

# Main plot, facet 1
# Helper so both plots share identical encodings and scales
make_plot <- function(d, xvar = dist_ord, line_w = 0.6, pt_size = 2) {
  ggplot(d, aes(x = {{ xvar }}, y = wmean, ymin = .lower, ymax = .upper,
                colour = Condition, fill = Condition, shape = Condition,
                group = Condition)) +
    geom_ribbon(alpha = 0.3, colour = NA) +
    geom_line(linewidth = line_w) +
    geom_point(size = pt_size) +
    scale_colour_manual(values = cond_cols) +
    scale_fill_manual(values = cond_cols) +
    scale_y_continuous(labels = scales::percent, limits = c(0.05, 0.55)) +
    labs(x = NULL, y = NULL) +
    # theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", hjust = 0))
}

# Prominent top panel
p_dist_total <- pred_dist_sum %>%
  filter(site == "total") %>%
  make_plot(line_w = 1.2, pt_size = 3.5) +
  facet_wrap(~site, labeller = labeller(site = get_site_labels())) +
  labs(
    x = "Distance (4+ = four or more)",
    y = "Percentage of total endowment given"
  ) +
  theme(
    strip.text       = element_text(size = 16, face = "bold", hjust = 0),
    # panel.background = element_rect(fill = "grey95", colour = NA),
    panel.border     = element_rect(colour = "grey30", fill = NA, linewidth = 0.8)
  )

# Smaller site panels
p_dist_sites <- pred_dist_sum %>%
  filter(site != "total") %>%
  make_plot() +
  facet_wrap(~site, nrow = 2, labeller = labeller(site = get_site_labels())) +
  labs() +
  guides(shape = "none", colour = "none", fill = "none")

p3b <- (p_dist_total + p_dist_sites) +
   plot_layout(widths = c(1.3, 2), guides = "collect") +
  plot_annotation(
    title    = "Allocations by distance",
    subtitle = "Average counterfactual predictions"
  ) &
  theme(legend.position = "top")
p3b

# Effect of cond by distance
eff_dist_all <- d_ctf_dist %>%
  group_by(Condition, dist_ord) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from=Condition, values_from = wmean) %>%
  mutate(cond_eff = revealed - anonymous) %>%
  median_qi(cond_eff) %>%
  mutate(site = "total")

eff_dist_site <- d_ctf_dist %>%
  group_by(site, Condition, dist_ord) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from=Condition, values_from = wmean) %>%
  mutate(cond_eff = revealed - anonymous) %>%
  median_qi(cond_eff)

eff_dist <- bind_rows(eff_dist_all, eff_dist_site) %>%
  mutate(site = factor(site, levels = rev(site_order)))

# Figure S3
p4 <- eff_dist %>%
  mutate(site = factor(site, levels = site_order2)) %>%
  ggplot(aes(x = cond_eff, xmin = .lower, xmax = .upper, y = dist_ord, colour = site)) +
  geom_pointrange(position = position_dodge(width = 0.2)) +
  labs(
    title = "Effect of reputational exposure by distance",
    x = "Additonal percentage of total endowment given in when identity is revealed",
    y = "Distance"
  ) +
  scale_x_continuous(labels = scales::percent) +
  geom_vline(xintercept = 0, linetype = 3) + scale_color_manual(values = get_site_colours()) +
  facet_grid(rows = vars(site)) +
  theme(
    legend.position = "none"
    # panel.spacing.y = unit(2, "lines")
  )
p4

# Contrasts to distance 1
dist_eff_all <-  d_ctf_dist %>%
  group_by(Condition, dist_ord) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from = dist_ord, values_from = wmean) %>%
  mutate(
    dist_eff_12 = `2` - `1`,
    dist_eff_13 = `3` - `1`,
    dist_eff_14 = `4+` - `1`
  ) %>%
  select(starts_with("dist_eff")) %>%
  pivot_longer(starts_with("dist_eff"), names_to = "contrast", values_to = "value") %>%
  median_qi(value) %>%
  mutate(site = "total")

dist_eff_site <-  d_ctf_dist %>%
  group_by(site, Condition, dist_ord) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from = dist_ord, values_from = wmean) %>%
  mutate(
    dist_eff_12 = `2` - `1`,
    dist_eff_13 = `3` - `1`,
    dist_eff_14 = `4+` - `1`
  ) %>%
  select(starts_with("dist_eff")) %>%
  pivot_longer(starts_with("dist_eff"), names_to = "contrast", values_to = "value") %>%
  median_qi(value)

dist_eff <- bind_rows(dist_eff_all, dist_eff_site) %>%
  mutate(site = factor(site, levels = rev(site_order))) %>%
  mutate(
    contrast = recode_values(
      contrast,
      "dist_eff_12" ~ "2",
      "dist_eff_13" ~ "3",
      "dist_eff_14" ~ "4+",
    )
  )

# Effect of distance, averaged across conditions
dist_eff %>%
  filter(Condition == "revealed" & contrast == "4+") %>%
  arrange(desc(value))


# Condition x distance ---------------------------------------------------------

# Interaction: change in eff of cond between dist 1 and 4+
d_ctf_dist %>%
  group_by(Condition, dist_ord) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from = Condition, values_from = wmean) %>%
  mutate(cond_eff = revealed - anonymous) %>%
  select(dist_ord, cond_eff) %>%
  pivot_wider(names_from = dist_ord, values_from = cond_eff) %>%
  mutate(interaction = `4+` - `1`) %>%
  median_qi(interaction)

# ...by site
d_ctf_dist %>%
  group_by(site, Condition, dist_ord) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from = Condition, values_from = wmean) %>%
  mutate(cond_eff = revealed - anonymous) %>%
  select(dist_ord, cond_eff) %>%
  pivot_wider(names_from = dist_ord, values_from = cond_eff) %>%
  mutate(interaction = `4+` - `1`) %>%
  median_qi(interaction) %>%
  print(n = Inf)

# interaction: size of condition at 4+ relative to 1.
d_ctf_dist %>%
  group_by(Condition, dist_ord) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from = Condition, values_from = wmean) %>%
  mutate(cond_eff = revealed - anonymous) %>%
  select(dist_ord, cond_eff) %>%
  pivot_wider(names_from = dist_ord, values_from = cond_eff) %>%
  mutate(interaction_rel = `4+` / `1`) %>%
  median_qi(interaction_rel)

#...by site
d_ctf_dist %>%
  group_by(site, Condition, dist_ord) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from = Condition, values_from = wmean) %>%
  mutate(cond_eff = revealed - anonymous) %>%
  select(dist_ord, cond_eff) %>%
  pivot_wider(names_from = dist_ord, values_from = cond_eff) %>%
  mutate(interaction_rel = `4+` / `1`) %>%
  median_qi(interaction_rel)

# Alter centrality -------------------------------------------------------------
rm(d_ctf_dist)

d_ctf_alter_cent <- gen_ctf(
  model_dir = file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v3", "cd_imput"),
  ctf_type = "alter_cent"
)

d_ctf_alter_cent$.epred <- rvar(
  gen_my_epreds(
  model_dir = file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v3", "cd_imput"),
  ctf_type = "alter_cent",
  n_dr = 50 # fine for plots
  )
)

preds_alter_all <- d_ctf_alter_cent%>%
  group_by(Condition, alter_cent_decile) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  mutate(site = "total")

preds_alter_site <- d_ctf_alter_cent %>%
  group_by(site, Condition, alter_cent_decile) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s))

preds_alter <- bind_rows(preds_alter_all, preds_alter_site)

preds_alter <- preds_alter %>%
  mutate(site = factor(site, levels = site_order)) %>%
  mutate(
    alter_cent_decile = as.character(alter_cent_decile),
    alter_cent_decile = recode_values(
      alter_cent_decile,
      "[0,0.1]" ~ 1,
      "(0.1,0.2]" ~ 2,
      "(0.2,0.3]" ~ 3,
      "(0.3,0.4]" ~ 4,
      "(0.4,0.5]" ~ 5,
      "(0.5,0.6]" ~ 6,
      "(0.6,0.7]" ~ 7,
      "(0.7,0.8]" ~ 8,
      "(0.8,0.9]" ~ 9,
      "(0.9,1]" ~ 10
    ),
    alter_cent_decile = factor(alter_cent_decile)
  )

pred_alter_sum <- preds_alter %>% median_qi(wmean, .width = 0.95)

# Main plot bottom facet
p_alter_total <- pred_alter_sum %>%
  filter(site == "total") %>%
  make_plot(line_w = 1.2, pt_size = 3.5, xvar = alter_cent_decile) +
  facet_wrap(~site, labeller = labeller(site = get_site_labels())) +
  labs(
    x = "Recipient Centrality (decile)",
    y = "Percentage of total endowment given"
  ) +
  theme(
    strip.text       = element_text(size = 16, face = "bold", hjust = 0),
    # panel.background = element_rect(fill = "grey95", colour = NA),
    panel.border     = element_rect(colour = "grey30", fill = NA, linewidth = 0.8)
  ) +
  guides(shape = "none", colour = "none", fill = "none")

# Smaller site panels
p_alter_sites <- pred_alter_sum %>%
  filter(site != "total") %>%
  make_plot(xvar = alter_cent_decile) +
  facet_wrap(~site, nrow = 2, labeller = labeller(site = get_site_labels())) +
  labs() +
  guides(shape = "none", colour = "none", fill = "none")

p4b <- (p_alter_total + p_alter_sites) +
   plot_layout(widths = c(1.3, 2), guides = "collect") +
  plot_annotation(
    title    = "Allocations by recipient centrality",
    subtitle = "Average counterfactual predictions"
  ) &
  theme(legend.position = "top")
p4b

# We report contrasts for top and bottom decile only. Recompute these with more draws.
d_ctf_alter_cent2 <- gen_ctf(
  model_dir = file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v3", "cd_imput"),
  ctf_type = "alter_cent",
  extreme_cent_dec_only = TRUE
)

d_ctf_alter_cent2$.epred <- rvar(
  gen_my_epreds(
    model_dir = file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v3", "cd_imput"),
    ctf_type = "alter_cent",
    n_dr = 200, # good for reporting values.
    extreme_cent_dec_only = TRUE
  )
)

# alter cent effect
d_ctf_alter_cent2 %>%
  group_by(alter_cent_decile) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from = alter_cent_decile, values_from = wmean) %>%
  mutate(alter_cent_eff = `(0.9,1]` - `[0,0.1]`) %>%
  median_qi(alter_cent_eff)

#...by site
d_ctf_alter_cent2 %>%
  group_by(site, alter_cent_decile) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from = alter_cent_decile, values_from = wmean) %>%
  mutate(alter_cent_eff = `(0.9,1]` - `[0,0.1]`) %>%
  median_qi(alter_cent_eff)

# interaction: cond * alter cent.
d_ctf_alter_cent2 %>%
  group_by(Condition, alter_cent_decile) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from = Condition, values_from = wmean) %>%
  mutate(cond_eff = revealed - anonymous) %>%
  select(alter_cent_decile, cond_eff) %>%
  pivot_wider(names_from = alter_cent_decile, values_from = cond_eff) %>%
  mutate(interaction = `(0.9,1]` - `[0,0.1]`) %>%
  median_qi(interaction)

# ...by site
d_ctf_alter_cent %>%
  group_by(site,Condition, alter_cent_decile) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from = Condition, values_from = wmean) %>%
  mutate(cond_eff = revealed - anonymous) %>%
  select(alter_cent_decile, cond_eff) %>%
  pivot_wider(names_from = alter_cent_decile, values_from = cond_eff) %>%
  mutate(interaction = `(0.9,1]` - `[0,0.1]`) %>%
  median_qi(interaction)

# interaction relative: cond * alter cent.
d_ctf_alter_cent2 %>%
  group_by(Condition, alter_cent_decile) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  pivot_wider(names_from = Condition, values_from = wmean) %>%
  mutate(cond_eff = revealed - anonymous) %>%
  select(alter_cent_decile, cond_eff) %>%
  pivot_wider(names_from = alter_cent_decile, values_from = cond_eff) %>%
  mutate(interaction_rel = `(0.9,1]` / `[0,0.1]`) %>%
  median_qi(interaction_rel)

# Ego cent ---------------------------------------------------------------------
rm(d_ctf_alter_cent)
gc()

# For plots
d_ctf_ego_cent <- gen_ctf(
  model_dir = file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v3", "cd_imput"),
  ctf_type = "ego_cent"
)

d_ctf_ego_cent$.epred <- rvar(
  gen_my_epreds(
  model_dir = file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v3", "cd_imput"),
  ctf_type = "ego_cent",
  n_dr = 50
  )
)

preds_ego_all <- d_ctf_ego_cent%>%
  group_by(Condition, ego_cent_decile) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  mutate(site = "total")

preds_ego_site <- d_ctf_ego_cent %>%
  group_by(site, Condition, ego_cent_decile) %>%
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s))

preds_ego <- bind_rows(preds_ego_all, preds_ego_site)

preds_ego <- preds_ego %>%
  mutate(site = factor(site, levels = site_order)) %>%
  mutate(
    ego_cent_decile = as.character(ego_cent_decile),
    ego_cent_decile = recode_values(
      ego_cent_decile,
      "[0,0.1]" ~ 1,
      "(0.1,0.2]" ~ 2,
      "(0.2,0.3]" ~ 3,
      "(0.3,0.4]" ~ 4,
      "(0.4,0.5]" ~ 5,
      "(0.5,0.6]" ~ 6,
      "(0.6,0.7]" ~ 7,
      "(0.7,0.8]" ~ 8,
      "(0.8,0.9]" ~ 9,
      "(0.9,1]" ~ 10
    ),
    ego_cent_decile = factor(ego_cent_decile)
  )


# Third facet in main plot
pred_ego_sum <- preds_ego %>% median_qi(wmean, .width = 0.95)

p_ego_total <- pred_ego_sum %>%
  filter(site == "total") %>%
  make_plot(line_w = 1.2, pt_size = 3.5, xvar = ego_cent_decile) +
  facet_wrap(~site, labeller = labeller(site = get_site_labels())) +
  labs(
    x = "Donor Centrality (decile)",
    y = "Percentage of total endowment given"
  ) +
  theme(
    strip.text       = element_text(size = 16, face = "bold", hjust = 0),
    # panel.background = element_rect(fill = "grey95", colour = NA),
    panel.border     = element_rect(colour = "grey30", fill = NA, linewidth = 0.8)
  ) +
  guides(shape = "none", colour = "none", fill = "none")

# Smaller site panels
p_ego_sites <- pred_ego_sum %>%
  filter(site != "total") %>%
  make_plot(xvar = ego_cent_decile) +
  facet_wrap(~site, nrow = 2, labeller = labeller(site = get_site_labels())) +
  labs() +
  guides(shape = "none", colour = "none", fill = "none")

p5b <- (p_ego_total + p_ego_sites) +
   plot_layout(widths = c(1.3, 2), guides = "collect") +
  plot_annotation(
    title    = "Allocations by donor centrality",
    subtitle = "Average counterfactual predictions"
  ) &
  theme(legend.position = "top")
p5b

# Preds by *counterfactual* ego centrality
preds_egocent_ctf <- d_ctf_ego_cent %>%
  group_by(ego_cent_decile) %>%
  # mean pred per draw and condition
  summarise(wmean = rvar_sum(.epred * w_s) / sum(w_s)) %>%
  # summarise(pred_value = mean(.epred)) %>%
  median_qi(wmean) %>%
  mutate(site = "total")

p_ctf_egocent <- preds_egocent_ctf %>%
  # mutate(site = factor(site, levels = eff_order)) %>%
  ggplot(aes(x = wmean, xmin = .lower, xmax = .upper, y = factor(ego_cent_decile, levels = levels(d_ctf$ego_cent_decile)))) +
  geom_pointrange(position = position_dodge(width = 0.2)) +
  labs(
    title = "Predicted allocation by donor centrality (counterfactual)",
    x = "Percentage of endowment given",
    y = NULL
  ) +
  scale_x_continuous(labels = scales::percent, limits = c(0.22, 0.33)) +
  theme(legend.position = "none") +
  scale_y_discrete(labels = get_site_labels())

# Combined plots -- S16 in article
p_obs_egocent / p_ctf_egocent



# Combined main plot
p3b / p5b / p4b +
  plot_annotation(
    title = "Allocations by distance, donor centrality, recipient centrality",
    subtitle = "Average counterfactual predictions"
  ) +
  plot_layout(guides = "collect")
