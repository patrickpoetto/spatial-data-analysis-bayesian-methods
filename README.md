# Spatial Data Analysis — Bayesian Methods

University coursework from the **Spatial Data Analysis** course,
covering the probabilistic/Bayesian-inference foundations (Sessions 4-5
assignment) that underpin the spatial modelling techniques taught later
in the module.

**Authors:** Corrà (5206191), Poetto (5201214), Shaboian (5217162), Bottino (5205740)
**Date:** 2024-03-03 (Bayesian Methods) / 2024-03-07 (Applied Models)

## Contents

| File | Description |
|---|---|
| `01_bayesian_methods_sda.Rmd` | Full reconstruction — Bayesian inference assignment (Sessions 4-5) |
| `02_applied_models_sda.Rmd` | Placeholder — see note below, source content was not recoverable |

## 1. Bayesian Methods for SDA (Sessions 4-5)

### Methodology and techniques applied

- **Conjugate Bayesian inference** for a Poisson likelihood: three
  competing Gamma priors (an informative prior built from English-text
  letter-frequency assumptions; a second, more diffuse prior informed by
  an external historical dataset; and a non-informative **Jeffreys'
  prior** derived analytically from the Poisson Fisher information).
- **Posterior summaries**: posterior mean, mode, and median (the latter
  approximated numerically via root-finding on the Gamma CDF, since the
  Gamma median has no closed form), plus 95% credible intervals and a
  posterior-probability hypothesis contrast (`P(mean < 25)`).
- **Prior/posterior predictive distributions**: analytically derived as
  **Negative Binomial** distributions (by marginalising the Poisson
  likelihood over the Gamma prior/posterior), used to answer predictive
  questions about a page not yet observed.
- **Generalised linear mixed models (GLMMs)**: a Poisson GLMM for
  repeated seizure counts (patient- and observation-level random
  effects) and a second-order random-walk smoother for a binomial
  rainfall time series — both specified as full Bayesian graphical
  models (founder parameters / deterministic and stochastic nodes /
  observed data).
- **Approximate Bayesian inference with INLA** (Integrated Nested
  Laplace Approximation): a random-intercept model on the `pulp`
  dataset (`faraway` package) comparing operators as a random vs. fixed
  effect via DIC/WAIC, and a Poisson regression with an exposure offset
  on a classic 10-pump failure/lifetime dataset.
- **Bayesian hierarchical modelling**: pump-specific failure rates
  modelled as exchangeable draws from a shared Gamma population
  distribution, with full conditional distributions derived analytically
  for a Gibbs sampler (and a Metropolis-Hastings step for the one
  hyperparameter without a closed-form conditional).
- **Survival analysis theory**: hazard/survivor functions, the Cox
  proportional-hazards model, and two extensions for *recurrent* event
  data — the Andersen-Gill proportional-intensity model and the
  Prentice-Williams-Peterson total-time model — discussed as
  alternatives to plain Poisson regression for modelling repeated pump
  failures.

### Datasets

- 16 hand-collected counts of the letter "u" per page, sampled from the
  Prefaces of the first four *A Game of Thrones* books.
- `GGOT.csv` — 31 additional "u"-per-page counts collected at the 2nd
  Valencia International Bayesian Summer School (July 2018); external
  file, not redistributed in this repository.
- `pulp` dataset (R package `faraway`) — paper brightness by shift
  operator, used for the INLA random-effects example.
- A classic 10-pump failure-count / lifetime dataset used for the
  Poisson/hierarchical/INLA pump-failure analysis.

### Key findings

- The Poisson MLE for "u"s-per-page (35) matched the sample mean, and
  all three priors (informative, diffuse, Jeffreys) converged to nearly
  identical posteriors (means 34.7-35.0) once updated with 16
  observations — the data dominated the choice of prior.
- The posterior probability that the true mean is below 25 was
  numerically negligible under all three priors (~1e-14 to 1e-19).
- For the operator-brightness data, the random-effect precision was
  orders of magnitude larger than the residual precision, suggesting
  minimal variation between operators — yet a fixed-effects model with
  dummy variables achieved lower DIC/WAIC, highlighting the tension
  between the exchangeability assumption and raw model fit.
- For pump failures, the Bayesian hierarchical model and the INLA
  Poisson GLMM with an exposure offset both pointed to a negative,
  credibly-nonzero intercept (`P(beta_0 < 0) = 0.993`) with pronounced
  between-plant variability, motivating the move from a single pooled
  rate parameter to a hierarchical/mixed-effects structure.

## 2. Applied Models (Assignment 3) — content not recoverable

The second deliverable, `02_applied_models_sda.Rmd`, was supplied for
this migration only as a ~1 MB raw-text extraction of a knitted R
Markdown HTML report. That extraction was read in full (100%, 988,622
characters / 4,179 lines) and found to contain **only the HTML page's
embedded front-end libraries** (jQuery, Bootstrap CSS, highlight.js /
pagedtable, Leaflet + Leaflet.markercluster) — truncated before the
report's actual `<body>` content (narrative, R code, tables, figures)
was ever reached. No analysis text, R code, or results of any kind were
present in the source; only two `<meta>` tags survived (`title`:
"Assignment 3"; `date`: "2024-03-07"), confirming the same author group
as above.

Rather than inventing an analysis that cannot be verified against any
source text, `02_applied_models_sda.Rmd` documents this limitation
explicitly and flags the one available clue: the presence of the
Leaflet/Leaflet.markercluster libraries in the page head suggests the
original report rendered an interactive map, consistent with a
"Spatial Data Analysis" applied-models deliverable. This file should be
replaced once the report can be re-exported from a readable copy of the
original HTML or `.Rmd` source.

## Tools

R, `stringr`, `faraway`, `INLA` (Integrated Nested Laplace Approximation),
base `stats` (`dpois`, `dgamma`, `pgamma`, `qgamma`, `pnbinom`, `uniroot`).

## License

MIT — see [LICENSE](LICENSE).
