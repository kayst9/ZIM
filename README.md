# ZIM

`ZIM` implements the Zero-Augmented Integration Method for integrating
individual-level internal data with external regression coefficient estimates.
It provides separate implementations for cohort and case-control studies and
supports common and different source populations.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("kayst9/ZIM")
```

Then load the package:

```r
library(ZIM)
```

## Main functions

- `zim_cohort()` fits ZIM under a cohort design.
- `zim_case_control()` fits ZIM under a case-control design.
- `zim_cohort_wkm()` and `zim_case_control_wkm()` use the corresponding
  working-model-based internal estimate.
- `zim_case_control_crossfit()` provides the cross-fitted case-control
  procedure.

The main functions take internal data, a target-model formula, working-model
specifications, and external coefficient estimates. See their help pages and
the bundled examples for the required argument structure.

## Examples and supplied data

The installed package includes four reproducible examples and their simulated
analysis inputs. Locate an example with:

```r
zim_example_path("cohort_one_population.R")
zim_example_path("case_control_one_population.R")
```

The examples may also be run from a clone of this repository. They use the
bundled simulation inputs by default and write results to a
`zim-example-output/` directory in the current working directory. See
`inst/examples/README.md` and `inst/examples/CASE_CONTROL.md` for details.

## Optional dependency

The default solvers use base R functionality. The optional
`gmm_solver = "gmm"` setting requires the [`gmm`](https://cran.r-project.org/package=gmm)
package.

## Verification

The four simulated example datasets are documented in
`inst/extdata/datasets.json`, which records their checksums, dimensions,
target model, external-summary specification, and successful saved-data replay.
Expected results and verification outputs are provided under `inst/examples/`.
The source scripts are loaded into separate internal environments so that the
cohort and case-control implementations retain their distinct helper functions.

## License

This software is available under the [MIT License](LICENSE.md).
