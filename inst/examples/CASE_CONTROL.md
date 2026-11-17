# Case-control examples

Each script loads one saved dataset from the bundled illustrative design and
runs one ZIM-GMM fit with unadjusted external summaries. Its
`--generate` mode generates and saves a new dataset before running the fit.

| Setting | One population | Two populations |
| --- | --- | --- |
| Script | `case_control_one_population.R` | `case_control_two_populations.R` |
| Source study | Bundled illustrative design | Bundled illustrative design |
| Saved dataset seed | 1 | 1 |
| Internal study | 1,000 cases + 1,000 controls | 1,000 cases + 1,000 controls |
| External study | 2,000 cases + 2,000 controls | 2,000 cases + 2,000 controls |
| Independent reference | None | 500 external-population controls |
| SNPs / blocks | 100 / one | 100 / one |
| Latent AR(1) correlation | 0.90 in both populations | 0.90 internal; 0.95 external/reference |
| Summary SNPs | All 100 | All 100 |
| Summary adjustment | None | None |

## Run

After installing and loading `ZIM`, run an example from an R session:

```r
source(zim_example_path("case_control_one_population.R"))
run_case_control_one_population()

source(zim_example_path("case_control_two_populations.R"))
run_case_control_two_populations()
```

Keep `case_control_example_helpers.R` beside the scripts. Only R's standard
libraries are required. The fixed calibration constants are included; neither
example needs the original simulation folders or calibration RDS files. Keep the
installed `extdata/` directory available with the examples.

By default, outputs are saved under `zim-example-output/<example_name>/` in the
current working directory. Optional
arguments are `--output-dir=/path/to/output` and `--zim-source=/path/to/program.R`.
The examples load the case-control implementation bundled with `ZIM`. The source
path and checksum are recorded for each run.

To generate and analyze a new dataset:

```sh
Rscript inst/examples/case_control_one_population.R --generate --data-file=/path/to/new-onepop.rds
Rscript inst/examples/case_control_two_populations.R --generate --data-file=/path/to/new-twopop.rds
```

To analyze an existing dataset, supply `--data-file=PATH` without `--generate`.
Generation uses the caller's RNG state; the scripts contain no seed-setting
statements or comments. The supplied datasets were generated with seed 1 set
outside the distributed scripts. `--generate` without `--data-file` writes a
new data file in the selected output directory.

## Data and models

The fixed design file `extdata/fixed_snp_design.csv` supplies SNP weights from
Uniform(-0.03, 0.07) and allele frequencies from Uniform(0.10, 0.40), originally
drawn with design seed 2026. The scripts read these vectors without resetting
the RNG. Both populations share these vectors.
Latent Gaussian AR(1) draws are transformed to genotype dosages 0/1/2 using the
original allele-frequency thresholds. The score `S` uses the first 100 SNPs.

The source-population disease model is:

```text
logit Pr(Y = 1 | S) = -0.014 + 0.015*S
```

The intercept is calibrated to give source-population disease prevalence of
approximately 50%.

Each study is sampled by generating source individuals, then retaining the
required cases and controls. Source batches contain 5,000 individuals for the
internal and external studies. Each completed sample is
shuffled. The two-population reference is sampled independently from external
controls, conditional on `Y = 0`, using source batches of 2,000; its outcome
column is removed before fitting ZIM. Generation order is internal, external,
then reference. The study case fractions are 50%, distinct from the source
disease prevalence.

Each external summary is the SNP slope from logistic `Y ~ xj`. The target model
is `Y ~ S`. Both scripts use:

```r
beta_method = "GMM"
inverse_method = "hbb_pca"
pc_var = 0.99
pc_max_condition = 200
pca_rule = "condition_or_variance"
pca_scale = "covariance"
ridge = 0
fixed_case_control = TRUE
Sigma_beta_beta = NULL
external_n_case = 2000
external_n_control = 2000
```

Covariance is estimated from internal data for one population, or internal data
and external-population reference controls for two populations. The generated
external individual-level data supply coefficient estimates only; no external
covariance matrix is computed or provided.

The scripts retain the source studies' GLM controls and GMM starting methods.
All scientific and analysis settings are exposed in each script's `settings`
list for generation mode; default analysis uses the settings stored in the
dataset. Each example fits only one setting to one dataset.

## Outputs

- `extdata/<example_name>.rds`: supplied analysis inputs, described in
  [the dataset README](../extdata/README.md). Generation mode writes a new data
  file in the output directory, or
  the path specified by `--data-file`, before fitting.
- `estimates.csv`: internal MLE and ZIM-GMM estimates for `S`, standard errors,
  95% confidence intervals, and two-sided p-values.
- `result.rds`: settings, sample-generation summaries, external coefficients,
  estimates, PC/convergence diagnostics, and provenance.
- `provenance.csv`: source path/checksum, example/helper checksums, R version,
  covariance source, and runtime.

These are single-dataset examples. Simulation bias, coverage, and power require
many datasets and are not estimated here. Rerunning in the same output directory
replaces the output files.

## Verification

The full-size case-control examples are checked with the canonical ZIM program:
one generation run followed by a fit loaded from each saved dataset.
`case_control_verification.csv` records the checks and resource use. Expected
seed-1 estimates are provided under `expected-results/`. Runtime and memory
depend on the machine and R configuration.
