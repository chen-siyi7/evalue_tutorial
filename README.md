# Evidence Aggregation with e-Values

R code for *Evidence Aggregation with e-Values: Dependence-Robust and Anytime-Valid Inference*, by Siyi Chen.

The `revision/` directory reproduces the revised numerical examples, Gaussian, Bernoulli and Poisson simulations, nine table fragments, and two ggplot2 figures. No external datasets are required. The original `simulation_fixed.R` and `simulation_sequential.R` scripts are retained unchanged at the repository root.

## Run the revised analyses

The workflow was tested with R 4.5.1, mvtnorm 1.4-2, ggplot2 4.0.3, and xml2 1.6.0. Install the three R packages if needed:

```r
install.packages(c("mvtnorm", "ggplot2", "xml2"))
```

From a terminal in the repository folder, run:

```sh
Rscript revision/run_all.R
```

The workflow reruns both original scripts, generates the revised analyses, builds tables and figures, and validates the outputs. It uses only R; neither Python nor a LaTeX installation is required. R must have Cairo graphics support (`capabilities("cairo")`) for SVG and PNG export. Helvetica is requested for figure lettering; font substitution on other systems can change appearance without changing the plotted data.

The full run takes several minutes on the test machine. Outputs are written under `revision/`. Existing output files are replaced when their corresponding stage runs. The manuscript source is not required or modified.

To check the supplied results without repeating the simulations:

```sh
Rscript revision/validate_revision.R
```

## Scripts

| Script in `revision/` | Purpose |
| --- | --- |
| `run_all.R` | Runs the complete workflow in order. |
| `revision_simulations.R` | Gaussian examples, fixed-sample and sequential comparisons, pilot/evaluation stopping designs, extended monitoring, the SPRT benchmark, twenty-component examples, continuation across studies, and log-growth calculations. |
| `non_gaussian_simulations.R` | Bernoulli and Poisson designs, exact spending boundaries, simulated moments, and sampling effort. |
| `build_tables.R` | Generates all nine table fragments from result files and static templates. |
| `build_figures.R` | Creates both ggplot2 figures and their SVG, PNG, RDS, and LaTeX exports. |
| `embed_svg.R` | Converts the plotted SVG paths into optional LaTeX vector fragments. |
| `validate_revision.R` | Checks probabilities, counts, Monte Carlo errors, boundaries, moments, quadrature, baseline reproduction, and the saved plotting data. |

Run any individual script with `Rscript revision/<script>.R`; its inputs must already be present. Table and figure builders can be run directly against the supplied CSVs.

## Results and manuscript correspondence

| Manuscript item | Result file(s) in `revision/results/` |
| --- | --- |
| Table 1: method comparison | Descriptive template in `revision/table_templates/methods.tex`; no simulation. |
| Table 2: worked examples | `worked_examples.csv`, `example_7_1_components.csv` |
| Table 3: fixed-sample Gaussian comparisons | `fixed_results.csv` |
| Table 4: planned sequential comparisons | `sequential_results.csv`, `group_sequential_boundaries.csv` |
| Table 5 and Figure 1: capped sampling effort | `stopping_design.csv`, `stopping_results.csv`, `stopping_plot_data.csv` |
| Table 6: unplanned extension | `extension_results.csv` |
| Table 7: optional continuation across studies | `continuation_results.csv` |
| Figure 2: expected log growth | `growth_rates.csv`, `growth_curves.csv`, `growth_quadrature_R.csv` |
| Table 8: twenty-component design | `fixed_results_K20.csv` |
| Table 9: Bernoulli and Poisson comparisons | `non_gaussian_results.csv`, `non_gaussian_design.csv`, `non_gaussian_boundaries.csv`, `non_gaussian_moments.csv` |
| Appendix SPRT benchmark | `sprt_results.csv` |
| Additional computational example | `workflow_example.csv` |

Full-precision CSVs contain rates, rejection counts, Monte Carlo standard errors, and sampling-effort summaries as applicable. Table fragments are in `revision/tables/`; their cross-references and mathematical macros refer to the manuscript. Table 1 is descriptive; the other eight tables are generated from the calculations. Figures are in `revision/figures/` as SVG and 600-dpi PNG, with saved ggplot objects and optional LaTeX fragments. LaTeX fragments are supplied for reuse, not compiled by this workflow.

## Reproducibility

Random seeds, replication counts, planning alternatives, dependence structures, thresholds, and weights are specified in the R scripts. The cap-selection pilot and evaluation samples use separate seeds. The discrete designs also record scenario seeds and parameters in `non_gaussian_design.csv`. Run logs and R/package versions are saved in `revision/results/`.

The supplied numerical results use the full simulation sizes, including 100,000 replications in the main fixed-sample and cross-study designs, and 50,000 in the main sequential, stopping, and non-Gaussian designs. Floating-point arithmetic, R/package versions, and graphics devices can affect last digits or rendering across systems; the supplied CSVs are the reference results for this revision. The validator includes independent calculations of group-sequential boundaries, discrete spending probabilities, distribution moments, and log-growth integrals.

Methods are compared under the assumptions and sampling rules stated in the manuscript. Unadjusted repeated p-value tests are included as diagnostic comparisons. A high rejection rate from such a method should not be interpreted as power at the nominal error level.

## Verification of this upload

A full run from a clean folder completed on October 2, 2026 and passed 1,303 checks. All twenty numerical CSVs match the manuscript's reference results, all nine table fragments match the manuscript tables, and both PNG figure exports are unchanged. The test also used an absolute script path containing spaces. Details and software versions are recorded in `revision/results/release_verification.json` and `reproduction_run.txt`.

## Updating the GitHub repository

Upload the contents of this folder to `chen-siyi7/evalue_tutorial`, preserving the directory structure. The root README describes the revised workflow; the two original root scripts retain their existing contents. The `revision/` folder contains the new code and reference outputs. No data download or credentials are needed to run the analyses.

## License

MIT; see `LICENSE`. The original scripts' license notice is preserved verbatim in `original_license/LICENSE`.
