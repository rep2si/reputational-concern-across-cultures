# Functions to extract and plot site-wise coefficients, with 95% CIs.
# Importantly this is designed to work on a dataframe of draws (as_draws_df),
# not the model itself. It is implemented this way so that we can work with
# thinned draws from large imputed models.

require(tidyverse)
require(tidybayes)
require(bayesplot)

# Return a parameter and CI for param, by site
single_coef_from_draws <- function(d_draws, param, site = "total") {
  # drop the "Intercept" var on brms' inner scale
  if (site == "total") {
  d_draws <- d_draws %>%
    select(!Intercept) %>%
    select(param) %>%
    rename(estimate = {{ param }}) %>%
    median_qi() %>%
    mutate(
      term = param,
      site = "total"
    ) %>%
    relocate(term, site)
  } else{
    d_draws <- d_draws %>%
      select(!Intercept) %>%
      select(param, paste0("r_site[", site, ",", str_remove(param, "^b_|^bsp_"),"]"))
    if(ncol(d_draws) != 2) stop("Regex borked panic now")
    d_draws[,"estimate"] <- d_draws[,1] + d_draws[,2]
    d_draws <- d_draws %>%
      select("estimate") %>%
    median_qi() %>%
      mutate(
        term = param,
        site = site
      ) %>%
    relocate(term, site)
  }
  return(d_draws)
}

multi_coefs_from_draws <- function(d_draws, params, site = "total") {
  d <- map(params, \(p) single_coef_from_draws(d_draws, p, site))
  d <- bind_rows(d)
  return(d)
}

get_all_site_coefs <- function(d_draws, params, sites = c("total", "az", "np", "te", "tm", "tz", "up", "wh", "wl")) {
  d <- map(sites, \(s) multi_coefs_from_draws(d_draws, params, s))
  d <- bind_rows(d)
  d <- d %>% arrange(factor(term, levels = params))
  return(d)
}

plot_coefs_from_draws <- function(coefs) {
  site_levs <- unique(coefs$site)
  total_included <- "total" %in% site_levs
  p <- coefs %>%
    filter(site != "total") %>%
    ggplot(aes(y = term, x = estimate, xmin = .lower, xmax = .upper, colour = site)) +
    geom_pointrange(position = position_dodge(width = 0.5)) +
    geom_vline(xintercept = 0, linetype = 3) +
    labs(title = "Params", x = NULL, y = NULL)

  if (total_included) p <- p +
    geom_pointrange(data = coefs %>% filter(site == "total"), position = position_nudge(y = 0.35)) +
    scale_colour_manual(
      values = get_site_colours(),
      labels = get_site_labels(),
      breaks = names(get_site_labels())
    )

  # Note: need to add labels for new effects as I move on
  p <- p +
    scale_y_discrete(
      labels = c(
        "b_Intercept" = "Intercept",
        "b_Conditionrevealed" = "Revealed",
        "b_ego_age_div" = "Ego age",
        "b_alter_age_div" = "Alter age",
        "b_ego_genderMale" = "Ego male",
        "b_alter_genderMale" = "Alter male",
        "b_anon_first" = "Anon first",
        "bsp_moego_cent_decile" = "Ego centrality",
        "bsp_moalter_cent_decile" = "Alter centrality",
        "bsp_modist_ord" = "Distance",
        "bsp_moego_cent_decileConditionrevealed" = "Ego centrality : Revealed",
        "moalter_cent_decile:Conditionrevealed" = "Alter centrality : Revealed",
        "modist_ordid:Conditionrevealed" = "Distance : Revealed"                                                                                                                                                                                                                     ),
    )

  return(p)
}

# Example use
# get_all_site_params(m5_draws, c("b_Conditionrevealed", "bsp_modist_ord", "bsp_moego_cent_decile")) %>%
  # plot_coeffs
