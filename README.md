# CQRM — Composite Quantile Regression Models for ncRNA Differential Expression

## Overview

CQRM is an R package implementing quantile-regression-based differential
expression (DE) tests for non-coding RNA (ncRNA) count data. ncRNA expression
is characterized by low counts, high variability, and multi-modal distributions
that violate the assumptions of conventional negative-binomial-based DE tools.
CQRM addresses these challenges by aggregating rank-score statistics across
multiple combinations of quantiles and selecting the quantile set that best
separates two conditions, while controlling the false discovery rate via a
permutation-calibrated MinP procedure.

---

## Installation

```r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("xzc0907/CQRM")

# Or install dependencies first:
install.packages(c("MASS", "lqmm", "quantreg"))
```

---

## Dependencies

| Package | Version | Used for |
|---------|---------|----------|
| `MASS` | ≥ 7.3 | `mvrnorm()` for null distribution sampling |
| `lqmm` | ≥ 1.5 | `lqm.fit.gs()`, `addnoise()`, `lqmControl()` (jitter-based fitting) |
| `quantreg` | ≥ 5.9 | `rq.wfit()` (single-pass LP solver used by `CQRM.approx`) |

---

## Quick Start

```r
library(CQRM)

set.seed(1)
G  <- 200          # genes
n  <- 20           # samples per group

# Simulate count data (NB distribution)
count1 <- matrix(rnbinom(G * n, mu = 5,  size = 0.5), G, n)
count2 <- matrix(rnbinom(G * n, mu = 10, size = 0.5), G, n)

# Fast approximate version (no jittering)
pv_approx <- CQRM.approx(count1, count2, tau_vec = c(0.2, 0.4, 0.6, 0.8), B = 1000)
fdr_approx <- p.adjust(pv_approx, method = "BH")

# Full CQRM (slower but more accurate)
pv_cqrm <- CQRM(count1, count2, tau_vec = c(0.2, 0.4, 0.6, 0.8), B = 1000, M = 50)
fdr_cqrm <- p.adjust(pv_cqrm, method = "BH")
```

---

## Functions

### `CQRM(count1, count2, tau_vec, B, seed, M)`
**Method: CQRM**

Full composite quantile regression test. Tests all \(2^K - 1\) non-empty
subsets of `tau_vec` quantiles and takes the minimum p-value (MinP). The MinP
is calibrated against draws from the asymptotic multivariate normal null
distribution of the joint score vector. `solve(subsigma)` is precomputed once
per weight pattern outside the B-loop.

| Argument | Default | Description |
|----------|---------|-------------|
| `count1` | — | Gene × sample count matrix, group 1 |
| `count2` | — | Gene × sample count matrix, group 2 |
| `tau_vec` | `c(0.2, 0.4, 0.6, 0.8)` | Quantile levels to aggregate |
| `B` | `10000` | Number of null draws for calibration |
| `seed` | `12345` | Random seed for reproducibility |
| `M` | `50` | Number of jitter replicates per gene per quantile |

**Returns:** Numeric vector of length G. Apply `p.adjust(..., method = "BH")` for FDR control.

---

### `CQRM.best(count1, count2, tau_vec, B, seed, M)`
**Method: CQRM-Best**

Tests each quantile in `tau_vec` independently (one chi-squared test per
quantile, scalar variance) and takes the MinP across quantiles. The null
distribution is generated analytically from independent normals (no matrix
solve), making this faster than `CQRM`.

Same arguments as `CQRM`.

---

### `CQRM.region(count1, count2, tau_vec, k, B, seed, M)`
**Method: CQRM-Region**

Tests all contiguous windows of `k` consecutive quantiles from `tau_vec` and
takes the MinP across windows. Intermediate between `CQRM.best` (`k = 1`) and
`CQRM` (all subsets). `solve(subsigma)` is precomputed once per window outside
the B-loop.

| Additional argument | Default | Description |
|--------------------|---------|-------------|
| `k` | `2` | Window width (number of consecutive quantiles per test) |

---

### `CQRM.approx(count1, count2, tau_vec, B, seed, zeta)`
**Method: CQRM-Approx**

Fast approximation of `CQRM`. Replaces the M-replicate jitter loop with a
single `rq.wfit` call (C-compiled LP solver from `quantreg`) on the
un-jittered Machado–Santos-Silva transformed response. Each quantile uses its
own tau-specific transformation. The permutation null and MinP calibration
proceed identically to `CQRM`. `solve(subsigma)` is precomputed once per
weight pattern outside the B-loop.

| Additional argument | Default | Description |
|--------------------|---------|-------------|
| `zeta` | `1e-5` | Minimum value for the log transformation (numerical stability) |

No `M` argument (no jittering). Requires `quantreg` instead of `lqmm`.

---

## Simulation & Comparison Code

The `data_code_paper/` folder contains reproducible code for the simulation
study in the accompanying paper:

| File | Contents |
|------|----------|
| `generate_data.R` | Generates a synthetic kidney-inspired dataset from the SimSeq `kidney` dataset |
| `comparison_methods.R` | Clean wrappers for all comparison methods (DESeq2, edgeR, limma-voom, Wilcoxon, NOISeq, ZIAQ, all CQRM variants) |
| `run_simulation.R` | Runs one complete simulation replicate and saves results as a formatted table |
| `CQRM_tutorial.Rmd` | R Markdown tutorial: data generation → method application → results display |

---

## Citation

If you use CQRM in your research, please cite:

> Xu Z. et al. (2026). Composite Quantile Regression Models for Differential
> Expression of Non-Coding RNA Count Data. *[Journal]*.

---

## License

MIT © 2026 Zhangchi Xu
