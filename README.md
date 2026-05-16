# Reproduction code for "Evidence Aggregation with *e*-Values"

This repository contains R code that reproduces the two simulation tables in

> *Evidence Aggregation with e-Values: Dependence-Robust and Anytime-Valid Inference.*

There are no external data dependencies. Each script is self-contained and writes its results to standard output.

## Contents

- `simulation_fixed.R` &mdash; reproduces Table 3 (fixed-sample size and power under three dependence structures and three signal regimes; M = 10<sup>5</sup> replicates).
- `simulation_sequential.R` &mdash; reproduces Table 4 (rejection rates under T = 10 interim looks without alpha-spending; M = 5 &times; 10<sup>4</sup> replicates).

## Requirements

- R 4.0 or later.
- The `mvtnorm` package (`install.packages("mvtnorm")`).

The scripts use only base R plus `mvtnorm`. No tidyverse or other heavy dependencies.

## Usage

From a shell, in the repository directory:

```
Rscript simulation_fixed.R
Rscript simulation_sequential.R
```

Each script prints a block of empirical rejection rates per scenario, with Monte Carlo standard errors. The fixed-sample script takes a few seconds; the sequential script takes roughly half a minute on a modern laptop.

## Reproducibility

Both scripts call `set.seed()` at the top:

- `simulation_fixed.R` uses seed `202605`.
- `simulation_sequential.R` uses seed `202606`.

Running on the same R version and `mvtnorm` version will produce identical numerical output. Across versions, the numbers will match to within Monte Carlo error (about 0.001 for the null sizes and 0.002 for the alternatives at the simulation sizes used here).

## Notes on the methods compared

Each script documents its methods in the header comment. In brief:

- Fisher's combination uses &minus;2 &sum; log p<sub>k</sub> on &chi;<sup>2</sup><sub>2K</sub>.
- Stouffer's combination uses (&sum; z<sub>k</sub>) / sqrt(K) for the independence variant and (&sum; z<sub>k</sub>) / sqrt(sum of &Sigma;) for the variant with known &Sigma;.
- The Cauchy combination test (Liu and Xie, 2020) uses the tan-based statistic with equal weights and the standard Cauchy null calibration.
- Mean *e*-value tests use the power calibrator f<sub>&kappa;</sub>(u) = &kappa; u<sup>&kappa;&minus;1</sup> with uniform weights; the outer-mixture variant averages over &kappa; &isin; {0.1, 0.3, 0.5, 0.7, 0.9}.
- The sequential mean *e*-process is built from one likelihood-ratio test martingale per stream against the point alternative &mu;<sub>alt</sub> = 0.20; the analyst rejects on first crossing of 1/&alpha; = 20.

## License

Released under the MIT License. See `LICENSE` for details.
