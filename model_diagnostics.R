library(tidyverse)
library(knitr)
library(kableExtra)
library(brms)
library(ordbetareg)
library(posterior)

sharepoint_path <- Sys.getenv("REP2SI_DIR")
plots_dir <- file.path(sharepoint_path, "rep2si-analysis", "global-analysis")

m2 <- readRDS(file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v2", "cd.rds"))
m2_gv0 <- readRDS(file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v2_gv0", "cd.rds"))
m2_gv50mr <- readRDS(file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v2_gv50mr", "cd.rds"))

m3 <- readRDS(file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v3", "cd.rds"))
m3_gv0 <- readRDS(file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v3_gv0", "cd.rds"))
m3_gv50mr <- readRDS(file.path(sharepoint_path, "rep2si-analysis", "global-analysis", "alloc_models", "v3_gv50mr", "cd.rds"))

# Order of tables in the output: M1, M2, M3, each v3 then v2
models <- list(
  m3 = m3,
  m2 = m2,
  m3_gv50mr = m3_gv50mr,
  m2_gv50mr = m2_gv50mr,
  m3_gv0 = m3_gv0,
  m2_gv0 = m2_gv0
)

escape_latex <- function(x) gsub("([_%&#$])", "\\\\\\1", x)

# Display names: base models (v2, v3) -> M1, _gv50mr -> M2, _gv0 -> M3.
# v2 models lack the condition * distance | ID slopes; v3 gets the plain name.
display_name <- function(x) {
  name <- case_when(
    str_detect(x, "_gv50mr") ~ "M2",
    str_detect(x, "_gv0") ~ "M3",
    TRUE ~ "M1"
  )
  if_else(
    str_detect(x, "^m2"),
    paste(name, "(without condition * distance \\textbar{} ID slopes)"),
    name
  )
}

# Labels follow cm_rand / cm_cor in coef_table.R
re_group_labels <- c(
  "ID" = "Donor",
  "ID_SUID" = "Donor Household",
  "AID" = "Recipient",
  "AID_SUID" = "Recipient Household",
  "RID" = "RA",
  "site" = "Site",
  "dyad_id_reponly" = "Repeated Dyad"
)

re_term_labels <- c(
  "Intercept" = "Intercept",
  "Conditionrevealed" = "Revealed",
  "ego_age_demean_div" = "Donor Age",
  "alter_age_demean_div" = "Recipient Age",
  "ego_genderMale" = "Donor Gender",
  "alter_genderMale" = "Recipient Gender",
  "ego_demean_logwealth" = "Donor Wealth",
  "alter_demean_logwealth" = "Recipient Wealth",
  "moego_cent_decileidEQegomo" = "Donor Centrality",
  "moalter_cent_decileidEQaltermo" = "Recipient Centrality",
  "modist_ordidEQdistmo" = "Distance",
  "anon_first" = "Anonymous First",
  "phi_Intercept" = "Phi (Shape of Beta Dist)",
  "dyad_repeated" = "Intercept"
)

lookup <- function(x, labels) if_else(x %in% names(labels), labels[x], x)

# "Conditionrevealed:alter_genderMale" -> "Revealed x Recipient Gender"
term_label <- function(term) {
  parts <- lookup(str_split_1(term, ":"), re_term_labels)
  paste(c(parts[parts == "Revealed"], parts[parts != "Revealed"]), collapse = " x ")
}

# "sd(a)" -> "SD a"; "cor(a,b)" -> "Corr a, b"
re_row_label <- function(rn) {
  map_chr(rn, function(r) {
    terms <- str_split_1(str_match(r, "^\\w+\\((.*)\\)$")[, 2], ",")
    labs <- map_chr(terms, term_label)
    if (str_starts(r, "cor")) paste("Corr", paste(labs, collapse = ", ")) else paste("SD", labs)
  })
}

# One LaTeX table per model: fixed effects first, then random effects grouped
# by grouping level. longtable, since donor-level correlations span pages.
model_tables <- imap(models, function(m, model_name) {
  fe <- summary(m)$fixed
  fe_tab <- tibble(
    Term = escape_latex(map_chr(rownames(fe), term_label)),
    Est = sprintf("%.2f", fe$Estimate),
    SE = sprintf("%.2f", fe$Est.Error),
    CI = sprintf("[%.2f, %.2f]", fe$`l-95% CI`, fe$`u-95% CI`),
    Rhat = sprintf("%.2f", fe$Rhat),
    Bulk = sprintf("%.0f", fe$Bulk_ESS),
    Tail = sprintf("%.0f", fe$Tail_ESS)
  )

  re <- summary(m)$random
  re <- re[order(match(names(re), names(re_group_labels)))]
  re_tab <- imap(re, function(t, group) {
    tibble(
      Term = escape_latex(re_row_label(rownames(t))),
      Est = sprintf("%.2f", t$Estimate),
      SE = sprintf("%.2f", t$Est.Error),
      CI = sprintf("[%.2f, %.2f]", t$`l-95% CI`, t$`u-95% CI`),
      Rhat = sprintf("%.2f", t$Rhat),
      Bulk = sprintf("%.0f", t$Bulk_ESS),
      Tail = sprintf("%.0f", t$Tail_ESS)
    )
  }) |> list_rbind()
  re_group_sizes <- set_names(map_int(re, nrow), lookup(names(re), re_group_labels))

  index <- c("Fixed Effects" = nrow(fe_tab), re_group_sizes)

  bind_rows(fe_tab, re_tab) |>
    kable(
      format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      align = c("l", rep("r", 6)),
      col.names = c("", "Est.", "SE", "95\\% CI", "Rhat", "Bulk ESS", "Tail ESS"),
      escape = FALSE,
      caption = sprintf("Estimated coefficients and sampler diagnostics for model %s, complete-case analysis.", display_name(model_name)),
      label = gsub("[^A-Za-z0-9]+", "-", paste("coef", model_name))
    ) |>
    kable_styling(latex_options = "repeat_header") |>
    pack_rows(index = index)
})

writeLines(
  paste(
    sprintf("{\\relscale{0.7}\\setlength{\\tabcolsep}{3pt}\n%s\n}", map_chr(model_tables, as.character)),
    collapse = "\n\\clearpage\n\n"
  ),
  file.path(plots_dir, "fit-diags.tex")
)
