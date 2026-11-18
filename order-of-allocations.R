# Exploratory plots showing that time may well have effects that we are not
# accounting for in our models.

library(tidyverse)
library(cowplot)
library(feasts)
library(tsibble)
source("./rep2si-common/source_all_submodule.R")

set_theme(theme_cowplot())

sharepoint_path <- Sys.getenv("REP2SI_DIR")
d <- readRDS(file.path(sharepoint_path, "rep2si-data/global-data/allocs.rds"))

rep2si_sharepoint <- Sys.getenv("REP2SI_DIR")
model_dir <- file.path(rep2si_sharepoint, "rep2si-analysis/global-analysis/alloc_models")
plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

# Ordering descriptives ---------------------------------------------------------

# Descriptive to check Condition * anon_first--clearly obtains in tz, less clear
# for others
d %>%
  group_by(ID, Condition, anon_first, site) %>%
  summarise(mean_alloc = mean(percGiven, na.rm = TRUE)) %>%
  pivot_wider(names_from = Condition, values_from = mean_alloc) %>%
  mutate(delta = revealed - anonymous) %>%
  group_by(site, anon_first) %>%
  summarise(mean_delta = mean(delta, na.rm = TRUE)) %>%
  pivot_wider(names_from = anon_first, values_from = mean_delta) %>%
  mutate(diff_mean_delta = `0` - `1`)

# Any overall drop in giving? DEFINITELY in TZ, which are the
# two sites with repeated measures.
d %>%
  ggplot(aes(x = q_order, y = percGiven, colour = site)) +
  geom_smooth() +
  labs(
    x = "Allocation #",
    y = "Percentage of endowment given"
  ) +
  scale_y_continuous(labels = scales::percent)

ggsave(file.path(plots_dir, "alloc_over_time.pdf"))


# Within condition drop?
# I first need a within-block clock. Doing this the naive way
within_blockify <- function(q_order, cutoff) {
  ifelse (q_order > cutoff, q_order - cutoff, q_order)
}

d <- d %>%
  mutate(
    q_order_within = case_when(
      site %in% c("te", "wl", "az", "tm") ~ within_blockify(q_order, 14),
      site %in% c("wh") ~ within_blockify(q_order, 12),
      site %in% c("tz", "np", "up") ~ within_blockify(q_order, 15)
    )
  )

d %>%
  ggplot(aes(x = q_order_within, y = percGiven, colour = site)) +
  geom_smooth() +
  labs(
    x = "Allocation # (within condition)",
    y = "Percentage of endowment given"
  ) +
  scale_y_continuous(labels = scales::percent)
ggsave(file.path(plots_dir, "alloc_over_time_within.pdf"))


# Within-subject heteroskedasticity
sliding_sd <- function(allocs, win_size) {
  if(win_size < 3) stop("Can't do an sd with fewer than 3 mesures")
  n_allocs <- length(allocs)
  num_windows <- n_allocs - win_size + 1
  num_windows <- ifelse(num_windows > 0, num_windows, 0)
  results <- seq_len(num_windows)
  results <- map(results, \(x){
    sd(allocs[x:(x + win_size - 1)], na.rm = TRUE)
  })
  return(unlist(results))
}

get_ordered_allocs <- function(id, df) {
  df %>% filter(ID == id) %>% arrange(q_order) %>% pull(percGiven)
}

slsd <- map(unique(d$ID), \(x) {
  ord_a <- get_ordered_allocs(x, d)
  result <- sliding_sd(ord_a, win_size = 5)
  ll <- length(result)
  tibble(
    ID = rep(x, ll),
    win_num = seq_len(ll),
    slsd = result
  )
})

d_slsd <- bind_rows(slsd)

d_slsd <- left_join(
  d %>% select(ID, site) %>% distinct,
  d_slsd)

d_slsd %>%
  ggplot(aes(x = win_num, y = slsd, colour = site)) +
  geom_smooth(alpha = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Moving window of 5 allocations",
    y = "Percentage of amount allocated, SD in window"
  )
ggsave(file.path(plots_dir, "sd_over_time.pdf"))




d_slsd %>%
  ggplot(aes(x = win_num, y = slsd, colour = site)) +
  geom_smooth(alpha = 0.2, method = "lm") +
  labs(
    x = "Moving window of 5 allocations",
    y = "SD in window"
  )


# Autocorrelation
dt <- d %>%
  as_tsibble(index = q_order, key = ID)

# Plot the ACF
d_acf <- dt %>%
  fill_gaps() %>%
  ACF(percGiven, lag_max = 12) %>%
  left_join(d %>% select(ID, site) %>% distinct, .)

d_mean_acf <- d_acf %>%
  group_by(site, lag) %>%
  # reframe handles the multiple columns returned by mean_cl_boot
  reframe(mean_cl_boot(acf)) %>%
  rename(
    mean_acf = y,
    conf.low = ymin,
    conf.high = ymax
  )

d_acf %>%
  ggplot(aes(x = factor(lag))) +
  geom_line(aes(y = acf, group = ID), colour = "black", alpha = 0.05) +
  geom_line(data = d_mean_acf, aes(x = factor(lag), y = mean_acf, group = NA), colour = "red") +
  geom_hline(yintercept = -1/28, linetype = 3) +
  facet_wrap(facets = vars(site))

d_mean_acf %>%
  ggplot(aes(x = lag, y = mean_acf, ymin = conf.low, ymax = conf.high)) +
  # geom_col(position = position_dodge(width = 1.1), width = 0.5) +
  # geom_col(position = position_dodge2(width = 1), width = 0.5) +
  # geom_errorbar(position = position_dodge2(width = 1), width = 0.5) +
  geom_lineribbon() +
  geom_hline(yintercept = 0, linetype = 3) +
  scale_x_continuous(breaks = 1:12) +
  facet_wrap(facets = vars(site))

ggsave(file.path(plots_dir, "mean_acf.pdf"))

# Examples from Nepal
d_acf %>%
  filter(site == "np") %>%
  ggplot(aes(x = factor(lag), y = acf)) +
  geom_line(aes(group = ID), colour = "black", alpha = 0.05) +
  geom_hline(yintercept = -1/28, linetype = 3) +
  geom_line(data = filter(d_acf, ID %in% c("np10301", "np01101", "np01904")),
            aes(group = ID, color = ID),
            linewidth = 1)


d %>%
  filter(ID %in% c("np10301", "np01101", "np01904")) %>%
  ggplot(aes(x = q_order, y = percGiven, colour = ID)) +
  geom_point() +
  geom_line() +
  facet_grid(rows = vars(ID))
