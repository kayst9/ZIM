# Saved example datasets

Each file contains one simulated example generated with dataset seed **1**:

| RDS file | Internal observations | External summary sample | Reference observations |
| --- | --- | --- | --- |
| `cohort_one_population.rds` | 2,000 | 4,000 | None |
| `cohort_two_populations.rds` | 2,000 | 4,000 | 500 external-population predictor observations |
| `case_control_one_population.rds` | 1,000 cases + 1,000 controls | 2,000 cases + 2,000 controls | None |
| `case_control_two_populations.rds` | 1,000 cases + 1,000 controls | 2,000 cases + 2,000 controls | 500 external-population controls |

All four use 100 SNPs in one block, unadjusted external single-SNP summary
coefficients, and a target model containing the polygenic risk score `S` only.
The two-population reference samples are independent of both study samples.

## Contents

Load a file with `readRDS()`. It contains:

- `internal_data`: `Y` (binary outcome), `S` (weighted score), and genotype
  dosages `x1` through `x100`.
- `reference_data`: `S` and the genotype dosages without `Y`, or `NULL` for one
  population.
- `beta_tilde`: the ordered vector of external unadjusted SNP coefficients.
- `settings`: target-model parameters, sample sizes, correlation settings, and
  ZIM options.
- `generation`: sample-generation summaries.
- `design`: fixed score weights (`omega`) and allele frequencies (`fre`).

External individual-level observations are generated to estimate `beta_tilde`,
then discarded. They are not needed to run ZIM. No external covariance matrix is
stored or supplied: the ZIM program estimates covariance from the internal data
and, for two populations, the reference data on each analysis run.

The compressed RDS files are the complete analysis inputs. The matching scripts
in `examples/` load them by default. Use `--generate --data-file=PATH` to create a
new dataset and analyze it. Simulation scripts contain no comments, seed-setting
statements, or fixed dataset seeds. New generation follows the caller's RNG state;
it does not automatically reproduce the provided seed-1 datasets.

## Generation metadata

`datasets.json` records the dataset seed, RNG kind, R version, file sizes and
checksums, and source program and generator checksums for the supplied files.
The original fixed SNP design was drawn using design seed 2026; its exact values
are provided in `fixed_snp_design.csv`. This preserves the score calibration while
allowing all four dataset realizations to use seed 1.

The seed was set only in the private preparation wrapper. Seed metadata remain
available outside the distributed simulation scripts. All data are simulated.
