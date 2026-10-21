# QRJM — Bayesian Quantile Joint Model for Multivariate Longitudinal and Time-to-Event Data

R/JAGS code accompanying the paper

> Kundu, D., Krishnan, S., Gogoi, M. P. and Das, K. (2024). **A Bayesian quantile joint modeling of multivariate longitudinal and time-to-event data.** *Lifetime Data Analysis*, 30(3), 680–699. https://doi.org/10.1007/s10985-024-09622-1

The repository implements a **quantile regression joint model (QRJM)** that links three longitudinal biomarkers to a time-to-event outcome, and allows **each biomarker to be modelled at its own quantile level**. It was developed for, and applied to, data from the Acute Lymphoblastic Leukemia (ALL) maintenance study conducted at Tata Medical Center, Kolkata.

---

## Motivation

Joint models for longitudinal and survival data usually describe the longitudinal process with a linear mixed model, i.e. they model the *mean* trajectory under a Gaussian assumption. Haematological biomarkers measured during ALL maintenance therapy are skewed and heterogeneous, and clinically the interest is often in the tails — e.g. how persistently *low* neutrophil counts relate to relapse — rather than the mean.

QRJM replaces the Gaussian longitudinal sub-model with a **linear quantile mixed model**, so that one can ask:

- How do covariate effects (including time-varying drug doses) change from one quantile of a biomarker to another?
- How does the association between the biomarkers and the hazard of relapse change across quantile levels?
- What happens when different biomarkers are examined at *different* quantiles simultaneously (e.g. a low quantile of ANC together with a high quantile of platelets)?

## The ALL maintenance application

Children with ALL received two standard maintenance drugs, **6-mercaptopurine (6MP)** and **methotrexate (MTX)**, for roughly two years, during which three blood counts were recorded repeatedly:

| Biomarker | Variable in code | Notes |
|---|---|---|
| Lymphocyte count | `WBC` | Derived as WBC − ANC on the original scale, then log-transformed |
| Absolute neutrophil count | `ANC` | Log scale |
| Platelet count | `Platelet` | Log scale |

All three are median-centred. Patients were then followed for about three more years and the time to relapse (if any) was recorded.

Covariates:

- **Time-varying:** 6MP and MTX doses (normalised by their medians), carried forward as step functions between visits.
- **Baseline:** age, gender, lineage, presenting WBC, NCI risk, bulky disease, CNS status, risk group, Day 8 / Day 35 response, morphological CR1, MRD and status.
- **Time:** linear and quadratic terms in (scaled) week.

Subjects with fewer than five longitudinal observations are excluded.

## Model

For subject $`i`$, biomarker $`k = 1, 2, 3`$ (lymphocyte, ANC, platelet) and quantile level $`\tau_k \in (0, 1)`$:

### Longitudinal sub-model

```math
y_{ik}(t) = \mathbf{x}_i(t)^{\top} \boldsymbol{\beta}_k(\tau_k) + U_{ik}(t) + \varepsilon_{ik}(t),
\qquad
\varepsilon_{ik}(t) \sim \mathrm{ALD}(0, \sigma_k, \tau_k)
```

so that, given the latent process, the $`\tau_k`$-th quantile of $`y_{ik}(t)`$ is $`\mathbf{x}_i(t)^{\top} \boldsymbol{\beta}_k(\tau_k) + U_{ik}(t)`$. Each biomarker can have its own $`\tau_k`$.

The asymmetric Laplace error is written as a normal–exponential mixture,

```math
\varepsilon_{ik}(t) = \kappa_{1k}\, e_{ik}(t) + \sqrt{\kappa_{2k}\, \sigma_k\, e_{ik}(t)}\; z_{ik}(t),
\qquad
e_{ik}(t) \sim \mathrm{Exp}(\text{mean } \sigma_k), \quad z_{ik}(t) \sim N(0, 1),
```

```math
\kappa_{1k} = \frac{1 - 2\tau_k}{\tau_k (1 - \tau_k)},
\qquad
\kappa_{2k} = \frac{2}{\tau_k (1 - \tau_k)},
```

which makes the likelihood conditionally Gaussian and gives straightforward Gibbs-type updates in JAGS.

The subject-specific term $`\mathbf{U}_i(t) = \left( U_{i1}(t), U_{i2}(t), U_{i3}(t) \right)^{\top}`$ is a trivariate Brownian-motion process,

```math
\mathbf{U}_i(t) - \mathbf{U}_i(s) \sim N_3\left( \mathbf{0}, (t - s) D \right), \qquad s < t,
```

with a Wishart prior on $`D^{-1}`$. It captures individual deviations over time and the correlation between the three biomarkers. In the code the process is evaluated on a subject-specific grid built from the Gauss–Kronrod nodes, so the same latent path enters both sub-models.

### Survival sub-model

```math
h_i(t) = h_0(t) \exp\left( \mathbf{w}_i^{\top} \boldsymbol{\theta} + \sum_{k=1}^{3} \psi_k\, m_{ik}(t) \right),
\qquad
m_{ik}(t) = \mathbf{x}_i(t)^{\top} \boldsymbol{\beta}_k(\tau_k) + U_{ik}(t)
```

- $`\mathbf{w}_i`$ are the baseline covariates, with effects $`\boldsymbol{\theta}`$.
- $`\psi_k`$ measures the association between the $`\tau_k`$-th quantile trajectory of biomarker $`k`$ and the hazard of relapse.
- The baseline log-hazard is a penalised cubic B-spline (P-spline), $`\log h_0(t) = \sum_{l} \gamma_l B_l(t)`$, with a second-order difference penalty and a gamma prior on the smoothing parameter.
- The cumulative hazard $`\int_0^{T_i} h_i(s)\, ds`$ is approximated with 15-point Gauss–Kronrod quadrature.
- The joint likelihood is coded in JAGS using the zeros trick.

### Priors

Vague normal priors (variance $`10^3`$) for $`\boldsymbol{\beta}_k`$, $`\boldsymbol{\theta}`$ and $`\boldsymbol{\psi}`$; a $`\mathrm{Gamma}(0.001, 0.001)`$ prior on each $`1/\sigma_k`$; and a Wishart prior on $`D^{-1}`$.

## Key findings from the paper

- Higher lymphocyte counts are associated with a shorter time to relapse.
- Higher neutrophil and platelet counts, taken together, are associated with a lower risk of relapse.
- 6MP lowers lymphocyte counts and MTX raises neutrophil counts at almost every quantile level examined.

The paper also includes simulation studies assessing the performance of the method.

## Repository contents

```
QRJM/
├── QRJM_main.R        # Data preparation, initialisation, QRJM fit in JAGS, posterior summaries
├── Visualization.R    # Comparisons across quantile combinations and figures for the paper
├── Data_Viz.R         # Exploratory plots of the biomarker data (Q-Q plot, density contours)
├── Plots/             # Figures (EPS/JPEG) produced for the manuscript
└── LICENSE            # MIT
```

### `QRJM_main.R`

The single parameter to set at the top is the vector of quantile levels, one per biomarker (lymphocyte, ANC, platelet):

```r
tau <- c(0.10, 0.10, 0.10)
```

The script then:

1. **Prepares the data** — scales time, transforms responses and covariates, drops subjects with < 5 visits, and saves `Transformed_DATA.rds`.
2. **Builds design matrices** — longitudinal design `X`, the design at the event time `XT`, the design at quadrature points `Xs`, and the B-spline bases `W2T` / `W2s` for the baseline hazard.
3. **Initialises in three stages**
   - separate univariate quantile mixed models with `lqmm::lqmm()` for each biomarker;
   - a longitudinal-only multivariate quantile model in JAGS (`QR.model.txt`);
   - a time-dependent Cox model on the fitted quantile trajectories to start $`\boldsymbol{\psi}`$ and $`\boldsymbol{\theta}`$.
4. **Fits the full QRJM** in JAGS (`QRJM.model.txt`) — by default 2 chains, 20,000 iterations, 10,000 burn-in, thinning 50, run in parallel.
5. **Summarises the posterior** — means, SDs and 95% credible intervals for $`\boldsymbol{\beta}`$, $`\sigma`$, $`D`$, $`\boldsymbol{\psi}`$ and $`\boldsymbol{\theta}`$, written to `QRJM(<tau>).csv`, plus DIC and $`p_D`$.

Each run writes to its own folder, `QRJM_<tau1>_<tau2>_<tau3>/`, containing the JAGS model files, the initial fit (`init_QR_fit_*.rds`), the final fit (`QRJM_fit_*.rds`) and the summary CSV.

### `Visualization.R`

Reads fitted models for a set of quantile combinations — by default

| Lymphocyte | ANC | Platelet |
|---|---|---|
| 0.25 | 0.25 | 0.25 |
| 0.25 | 0.75 | 0.75 |
| 0.50 | 0.50 | 0.50 |
| 0.75 | 0.25 | 0.25 |
| 0.75 | 0.75 | 0.75 |

and produces:

- overlap of credible intervals between pairs of quantile combinations (how different the fitted models are);
- association parameters $`\psi_k`$ across quantile combinations;
- 6MP and MTX effects on each biomarker across quantile combinations;
- significance maps for baseline covariates in the longitudinal and survival sub-models;
- $`\hat{R}`$ convergence diagnostics, trace plots and posterior densities;
- elements of the Brownian-motion covariance matrix $`D`$ across quantile levels.

Run `QRJM_main.R` once for each row of `T_set` before running this script.

### `Data_Viz.R`

Exploratory graphics of subject-level mean ANC and platelet counts, including a multivariate normal Q-Q plot and 2-D density contours split by relapse status, labelled at the empirical quantiles.

## Requirements

- **R** (≥ 4.0 recommended)
- **JAGS** (≥ 4.3) — install separately from <https://mcmc-jags.sourceforge.io/>
- R packages:

```r
install.packages(c(
  "rjags", "R2jags", "jagsUI", "lqmm", "survival", "splines", "lme4",
  "MASS", "dplyr", "data.table",
  "ggplot2", "gridExtra", "ggpubr", "RColorBrewer", "reshape2", "lattice",
  "readxl", "stringr"
))
```

The `jags()` call in `QRJM_main.R` uses the `jagsUI` interface (`parallel = TRUE`, `n.adapt`); load `jagsUI` after `R2jags`, as in the script, so its version of `jags()` is the one used.

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/Damitri411/QRJM.git
   ```
2. Edit the paths at the top of each script (`fp`, `path`, `path_out`) to point to your data and output folders. They are currently hard-coded to the author's machine.
3. Place the input data as `data_0.rds` in the data folder (see *Data availability* below for the expected structure).
4. Set `tau` in `QRJM_main.R` and run it. Repeat for each quantile combination of interest.
5. Run `Visualization.R` to compare fits across quantile combinations.

**Runtime.** The full joint model is computationally heavy (three latent processes on a fine time grid, quadrature for the cumulative hazard). Expect a run to take hours rather than minutes per quantile combination; parallel chains help.

### Adapting to other data

To use the code on another study, the input data frame needs, in long format, one row per subject-visit with:

- a subject identifier (`Child.Id`, `id_new`) and visit index (`time_id`);
- visit time (`Week`);
- the three longitudinal outcomes;
- time-varying covariates (here `X6_mp`, `MTX`);
- baseline covariates (repeated on every row);
- event time (`Time_to_event`), event indicator (`Relapsed`) and last follow-up (`Last_follow_up_date`).

The transformation block, the covariate lists in the `coxph()` and model formulas, and the parameter names in `Stats_pos()` would need to be changed to match.

## Data availability

The ALL maintenance data from Tata Medical Center, Kolkata are **not included** in this repository, because they are confidential patient-level clinical data. The code is provided so that the methodology can be inspected, reproduced on comparable data, or adapted to other studies.

## Citation

If you use this code, please cite:

```bibtex
@article{Kundu2024QRJM,
  author  = {Kundu, Damitri and Krishnan, Shekhar and Gogoi, Manash Pratim and Das, Kiranmoy},
  title   = {A {B}ayesian quantile joint modeling of multivariate longitudinal and time-to-event data},
  journal = {Lifetime Data Analysis},
  year    = {2024},
  volume  = {30},
  number  = {3},
  pages   = {680--699},
  doi     = {10.1007/s10985-024-09622-1}
}
```

## Related work

- Kundu, D. and Das, K. (2024). A quantile-regression approach to bivariate longitudinal joint modeling. *Journal of Statistical Research*, 58(1), 111–130.
- Kundu, D. (2024). *Bayesian joint modeling of multivariate longitudinal and event-time outcomes with applications to ALL maintenance studies.* PhD thesis, Indian Statistical Institute, Kolkata.

## License

Released under the [MIT License](LICENSE).

## Contact

Damitri Kundu — [@Damitri411](https://github.com/Damitri411)
