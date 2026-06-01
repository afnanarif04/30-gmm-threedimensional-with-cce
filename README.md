# Dynamic Three-Dimensional Common Correlated Effects Estimation

## Research Overview

This replication package accompanies a paper that proposes a **dynamic
three-dimensional common correlated effects (CCE) mean group estimator**
for heterogeneous panel data in which outcomes are observed for ordered
pairs of units (such as bilateral trade flows, cross-border lending, or
migration) over time.

### The Estimation Problem

Three-dimensional panel data carry cross-sectional dependence that
operates simultaneously at three levels: global shocks that affect all
units, shocks specific to each origin unit (the sending side), and shocks
specific to each destination unit (the receiving side). Existing estimators
handle only the global level (Pesaran, 2006) or are limited to static
specifications (Kapetanios, Serlenga and Shin, 2021). This paper extends
the framework to the dynamic case, where the outcome variable depends on
its own lagged value as well as on weakly exogenous regressors.

### The Estimator

The proposed estimator augments each unit-pair regression with the current
and lagged values of three cross-sectional averages:

- the **grand average** over all pairs (absorbs global factors)
- the **origin-margin average** over all destinations for a given origin
  (absorbs origin-specific factors)
- the **destination-margin average** over all origins for a given
  destination (absorbs destination-specific factors)

The mean group estimator averages the pair-specific coefficients. The
procedure is closed form at the pair level, requires no iteration, and
does not require the number of factors to be known in advance.

### Key Theoretical Results

- The estimator is consistent and asymptotically normal under joint
  expansion of the origin dimension N, the destination dimension M,
  and the time dimension T.
- The convergence rate is governed by the **smaller of N and M**, not
  the pair count N × M, because the dominant source of sampling variation
  is heterogeneity along the two cross-sectional margins.
- Valid inference requires a **margin-based variance estimator** that
  accumulates dispersion across the two margin means. The standard
  pair-level variance estimator is inconsistent in this setting.

---

## Monte Carlo Simulation

The simulation assesses finite-sample bias, root mean squared error, and
empirical test size across four sample configurations. The data-generating
process follows the proposed hierarchical factor structure with
heterogeneous slopes, dynamic dependence through the lagged outcome, and
factor loadings that vary across pairs.

### Sample configurations

| Configuration | N | M | T |
|---|---|---|---|
| C1 (small) | 20 | 20 | 50 |
| C2 | 25 | 25 | 100 |
| C3 | 40 | 40 | 100 |
| C4 (large) | 50 | 50 | 150 |

### Competitors

The simulation compares four estimators:
1. **Proposed** — dynamic three-dimensional CCE mean group (this paper)
2. **Static 3D** — Kapetanios, Serlenga and Shin (2021) — omits dynamics
3. **Naive 2D dynamic** — Chudik and Pesaran (2015) — ignores the
   destination factor level
4. **Pooled FE** — imposes slope homogeneity and additive fixed effects

### Summary of simulation findings

The short-run mean group estimator shows negligible bias and declining root
mean squared error as the sample dimensions grow. Test size for the
proposed estimator is near-nominal at moderate dimensions. All three
misspecified competitors show rejection rates rising toward one as sample
dimensions expand, confirming that omitting dynamics, omitting a factor
level, or imposing additive structure produces growing distortion.

### Scripts (`code/simulation/`)

| Script | What it does |
|--------|-------------|
| `01_dgp_functions.R` | Generates three-dimensional panel data under the hierarchical factor structure; builds N×M×T arrays |
| `02_estimators.R` | Implements all four estimators: proposed dynamic 3D-CCE, static 3D-CCE, naive 2D dynamic CCE, and pooled fixed effects |
| `03_simulation.R` | Runs the parallel simulation loop across all four sample configurations (2,000 replications each) |
| `04_summarise.R` | Computes bias, RMSE, size, and coverage; exports Tables 1 and 2 as CSV |

**Estimated runtime:** approximately 45 minutes on a 16-core machine.
Single-core runtime will be longer in proportion.

---

## Empirical Application

The empirical illustration applies the estimator to a dynamic gravity
model of bilateral export flows. The three-dimensional structure arises
naturally: each flow is indexed by an exporting economy (origin), an
importing economy (destination), and a year.

### Data

- **Bilateral exports:** CEPII BACI database (HS92 classification for
  1995–2016; HS17 classification for 2017–2024). Flows are aggregated
  over all product codes to produce bilateral totals.
- **Income (GDP):** World Bank World Development Indicators
  (series NY.GDP.MKTP.CD, current USD). Chinese Taipei (Taiwan) GDP
  is sourced from the IMF World Economic Outlook.

### Panel dimensions

| | |
|---|---|
| Origin economies (N) | 18 |
| Destination economies (M) | 19 |
| Time period | 1995–2024 (T = 30 years) |
| Balanced pairs | 322 |
| Total observations | 9,660 |

### Key empirical findings

- The raw pair-level first-order autocorrelation of log bilateral exports
  exceeds 0.93, indicating strong apparent persistence.
- After absorbing the three-level factor structure with the proposed
  estimator, the factor-adjusted persistence coefficient is approximately
  0.095 and is statistically significant at the one per cent level. This
  means that roughly ninety per cent of apparent bilateral persistence
  reflects common shocks rather than genuine bilateral relationship
  inertia — a finding that a static estimator cannot recover.
- The cross-sectional dependence statistics confirm the hierarchical
  structure: adding only the grand average reduces but does not eliminate
  dependence; the full three-margin augmentation drives both origin-margin
  and destination-margin statistics to near zero.
- Short-run and long-run income elasticities are imprecisely estimated
  after removing the factor structure, consistent with the structural
  gravity interpretation that GDP variation is largely absorbed by the
  margin-level factor proxies.

### Scripts (`code/empirical/`)

| Script | What it does |
|--------|-------------|
| `E01_data_from_excel.R` | Reads `data/Empirical_Data_3DCCE_Paper.xlsx`; constructs N×M×T arrays for the outcome and regressors |
| `E02_cd_diagnostics.R` | Three-step cross-sectional dependence diagnostic (Pesaran 2015 CD statistic and BKP 2016 exponent at each augmentation stage) |
| `E03_estimation.R` | Dynamic 3D-CCE mean group estimation with margin-based standard errors and delta-method long-run standard errors |
| `E04_results_excel.R` | Consolidates all results into `output/EMPIRICAL_RESULTS.xlsx` |

**Estimated runtime:** approximately 10 minutes on a standard desktop.

---

## Software Requirements

- **R** version 4.3.0 or later (tested on R 4.5.2)
- **RStudio** (recommended)

### Required R Packages (all on CRAN)

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "openxlsx",
  "parallel"
))
```

---

## Repository Structure

```
.
├── README.md
├── LICENSE
├── DATA_AVAILABILITY_STATEMENT.md
├── code/
│   ├── simulation/
│   │   ├── 01_dgp_functions.R
│   │   ├── 02_estimators.R
│   │   ├── 03_simulation.R
│   │   └── 04_summarise.R
│   └── empirical/
│       ├── E01_data_from_excel.R
│       ├── E02_cd_diagnostics.R
│       ├── E03_estimation.R
│       └── E04_results_excel.R
├── data/
│   ├── README_data.md
│   └── Empirical_Data_3DCCE_Paper.xlsx   ← processed bilateral panel
└── output/                               ← populated when code runs
    └── .gitkeep
```

---

## Replication Instructions

### Step 1 — Set the working directory

Every script sets its working directory automatically when run
interactively in RStudio. If running non-interactively, set the working
directory to the folder containing the script before running.

### Step 2 — Run the Monte Carlo simulation

Run the simulation scripts in order from the `code/simulation/` folder:

```r
source("code/simulation/01_dgp_functions.R")
source("code/simulation/02_estimators.R")
source("code/simulation/03_simulation.R")   # runs in parallel
source("code/simulation/04_summarise.R")
```

Output is saved to `output/mc_table1_bias_rmse.csv` and
`output/mc_table2_size_coverage.csv`.

### Step 3 — Run the empirical estimation

The processed data file `data/Empirical_Data_3DCCE_Paper.xlsx` must
be present. Run the empirical scripts in order:

```r
source("code/empirical/E01_data_from_excel.R")
source("code/empirical/E02_cd_diagnostics.R")
source("code/empirical/E03_estimation.R")
source("code/empirical/E04_results_excel.R")
```

Output is saved to `output/cd_diagnostics.csv`,
`output/estimation_results.csv`, and `output/EMPIRICAL_RESULTS.xlsx`.

### Step 4 — Check output

| Output file | Contents |
|-------------|----------|
| `output/mc_table1_bias_rmse.csv` | Bias and RMSE by configuration and estimator |
| `output/mc_table2_size_coverage.csv` | Size and coverage by configuration and estimator |
| `output/cd_diagnostics.csv` | CD statistics and BKP exponents at each augmentation stage |
| `output/estimation_results.csv` | Mean group estimates, standard errors, t-statistics |
| `output/EMPIRICAL_RESULTS.xlsx` | Full consolidated results workbook (four sheets) |

---

## Data Availability

The processed empirical dataset is included in this repository as
`data/Empirical_Data_3DCCE_Paper.xlsx`. This file contains the bilateral
export panel and the GDP series, assembled from publicly available sources.
Raw bilateral trade data from CEPII BACI are not redistributable; see
`data/README_data.md` for download instructions.

---

## License

MIT License. See `LICENSE` for details.

---

## Notes on Reproducibility

- Random seeds are set at the top of `03_simulation.R`.
- Results may vary slightly across platforms due to floating-point
  arithmetic and parallel random number generation.
- The lag-truncation order `p_T` is determined automatically by the
  `get_pT()` function in `E03_estimation.R`. It evaluates to zero for
  the current panel (T = 30, N = 18, M = 19) and can be inspected or
  overridden manually if the panel dimensions change.

---

*Author name, institutional affiliation, paper title, and journal name
will be added after acceptance.*
