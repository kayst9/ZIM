# Cohort examples

The two additional [case-control examples](CASE_CONTROL.md) cover one- and
two-population settings under the same compact design.

Each script loads a saved dataset and runs one ZIM-GMM analysis with
unadjusted external summaries. Its `--generate` mode creates a new
internal/external dataset pair, saves the analysis inputs, and runs the fit. These are
compact, runnable illustrative examples.

| Setting | One population | Two populations |
| --- | --- | --- |
| Script | `cohort_one_population.R` | `cohort_two_populations.R` |
| Source study | Bundled illustrative design | Bundled illustrative design |
| Saved dataset seed | 1 | 1 |
| Internal cohort | 2,000 | 2,000 |
| External cohort | 4,000 | 4,000 |
| Independent reference sample | None | 500 external-population predictor observations |
| SNPs / blocks | 100 / one | 100 / one |
| Latent AR(1) correlation | 0.90 in both studies | 0.90 internal; 0.95 external and reference |
| Summary SNPs | All 100 | All 100 |
| External summaries | Unadjusted | Unadjusted |

## Run

After installing and loading `ZIM`, run an example from an R session:

```r
source(zim_example_path("cohort_one_population.R"))
run_cohort_one_population()

source(zim_example_path("cohort_two_populations.R"))
run_cohort_two_populations()
```

Keep `cohort_example_helpers.R` beside the two scripts. Only R's standard
libraries are needed. The examples contain their fixed calibration constants;
the original simulation folders and calibration files are not required. Keep the
installed `extdata/` directory available with the examples.

By default, each example writes to
`zim-example-output/<example_name>/` in the current working directory. To
choose a different destination when using the scripts through `Rscript`:

```sh
Rscript inst/examples/cohort_one_population.R --output-dir=/path/to/onepop-output
Rscript inst/examples/cohort_two_populations.R --output-dir=/path/to/twopop-output
```

The examples use the cohort implementation bundled with `ZIM`. An explicit
`--zim-source=/path/to/program.R` overrides this selection. The exact loaded
source path and MD5 checksum are saved with results.

To generate a new dataset and analyze it:

```sh
Rscript inst/examples/cohort_one_population.R --generate --data-file=/path/to/new-onepop.rds
Rscript inst/examples/cohort_two_populations.R --generate --data-file=/path/to/new-twopop.rds
```

To analyze an existing dataset, supply `--data-file=PATH` without `--generate`.
Generation uses the caller's RNG state; the scripts contain no seed-setting
statements or comments. The supplied datasets were generated with seed 1 set
outside the distributed scripts. `--generate` without `--data-file` writes a
new data file in the selected output directory.

## Data and models

The fixed design file `extdata/fixed_snp_design.csv` contains the original SNP
weights from Uniform(-0.03, 0.07) and allele frequencies from Uniform(0.10, 0.40).
These vectors were originally drawn using design seed 2026; they are read from
the file without resetting the RNG. Genotypes are dosages 0/1/2
obtained from correlated latent normal draws using the original simulation's
allele-frequency thresholds. The AR(1) correlation in the table describes the
latent normal variables. The score `S` is the weighted sum of the first 100 SNPs.

Both studies use the target model:

```text
logit Pr(Y = 1 | S) = -0.014 + 0.015*S
```

The intercept is calibrated to give approximately 50% outcome prevalence.
Internal data are generated before external data. In the two-population example,
the independent external-population reference is generated next and contains
predictors only.

Each external model is `Y ~ xj`; only the SNP coefficient is passed to ZIM. The
internal target is `Y ~ S`. Both examples set:

```r
beta_method = "GMM"
inverse_method = "hbb_pca"
pc_var = 0.999
pc_max_condition = 500
pca_rule = "condition_or_variance"
pca_scale = "covariance"
ridge = 0
Sigma_beta_beta = NULL
n_external = 4000
```

The selection rule retains the smaller number of components allowed by the
99.9% variance target and the condition-number limit of 500. The retained
variance can therefore be below 99.9%.

Covariance is estimated by the cohort program from internal data for one
population, or internal plus reference data for two populations. External
individual-level data are used only to obtain coefficient estimates. No
external covariance matrix is computed or supplied.

The `settings` list near the top of each script exposes the design and analysis
parameters. The current `pc_var`, `pc_max_condition`, and `pca_rule` settings
apply to both saved-data and generation runs. Other settings are loaded from
the saved dataset during a default run. Thus, changing the PCA controls does
not require regenerating or modifying the data.

To demonstrate a smaller summary set, change `m` to a value below 100 and use
`--generate`. This selects the first `m` SNPs, while `S` still uses all 100 SNPs.
Each script fits a single setting to one dataset.

## Output

- `extdata/<example_name>.rds`: supplied analysis inputs, described in
  [the dataset README](../extdata/README.md). Generation mode writes a new data
  file in the output directory, or
  the path specified by `--data-file`, before fitting.
- `estimates.csv`: internal MLE and ZIM-GMM estimates, standard errors, 95% Wald
  confidence intervals, and two-sided p-values for the score.
- `result.rds`: design settings, sample summaries, external SNP coefficients,
  estimates, retained-PC and convergence diagnostics, and source provenance.
- `provenance.csv`: program path/checksum, example/helper checksums, R version,
  covariance source, and elapsed time.

These are single-dataset demonstrations. Their estimates and p-values are not
simulation bias, coverage, or power summaries. Running a script again with the
same output directory replaces those three output files.

## Verification

All four examples are checked with the canonical ZIM programs: one generation
run followed by a fit loaded from the saved inputs. `verification.csv` records
the cohort checks and resource use. Expected seed-1 estimates are provided under
`expected-results/`. Runtime and memory depend on the machine and R configuration.

The original verification CSVs record dataset-generation checks.
`pca_verification.csv` records saved-data replay with the current PCA settings.
All four examples reproduce the existing expected estimates without modifying
the supplied datasets.
