# reputational-concern-across-cultures

This repository contains code used for the manuscript "Reputational Concern Across Cultures". Note that the data cannot be provided for reasons explained in the article. The code here if provided for review purposes only.

The full workflow is complex. The code here focuses on the models we fit and results we report from them. If the used to fit ERGMS to the social network data is also required (this is only used for imputation purposes), please contact the author.

## Fitting models

Code used to fit the models is contained in the `fit-models` subdirectory. This assumes a SLURM installation. Running this code locally is prohibitively expensive from a computational perspective, except for the most trivial models---fit time is typically 4 to 5 days per model.

The most important code is in the `prepare-model-XXX.R` scripts; this contains the model specifications. Reviewers may also want to look at the `fit-models-XXX.R` scripts, which specify priors. The rest of the code is simply engineered to run mulitple models in parallel on a SLURM cluster.

## Outputs (tables and plots)

- `dg-descriptives.R` generates the histograms, Ligon-Schehter style plot, table 2

- `marginal-effects-imput.R` computes marginal effects, produces main plot and values reported in the article, and some of the supplemetary plots

- `mini-coef-table.R` produces Table 1

- `coef-table.R` produces the full coefficient tables in SI B

- `plot-coefs.R` plots figures S6 to S10

- `compare-pooled-sitewise.R` is used for SI C (sensitivity to partial pooling)

- `no-monotonic.R` is used for SI D (models without monotonicity constraint)

- `plot-random-effects.R` generates the plots for SI F (Unexplained heterogeneity)

- `compare-imp-noimp.R` produces S26 -- S28

- `model-diagnostics.R` procduces tables S11 -- S16

- `imputation-diagnostics.R` produces table S18

- `prior-sens.R` is used for prior sensitivity checks (S29 -- S31)

- `compare-comp-drop` is used for comparing models with and without comprehension checks (S32-S33)

- `comprehension-model.R` shows results for models with enhanced comprehension checks (NP and UP, S35)

- `order-of-allocations.R` produces the exploratory plots for SI N

- `unknown-recipients.R` generates table S19 for SI O

- `comp-v2-v3.R` is used for SI P (comparing models with and without one participant-level random slope)

- `strand-models.R` is used for SI Q (STRAND models as robustness check)
