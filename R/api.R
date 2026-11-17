#' Fit a Zero-Augmented Integration Method model for cohort data
#'
#' This wrapper calls the cohort implementation bundled with the package.
#' Supply the arguments described in the package README and the example
#' scripts. The return value contains the ZIM estimate, its standard error,
#' the internal estimate, covariance information, and diagnostics.
#'
#' @param ... Arguments passed to the cohort ZIM implementation.
#' @return A list containing fitted estimates and diagnostics.
#' @export
zim_cohort <- function(...) {
  .zim_impl_env("cohort")$zim_cohort(...)
}

#' Fit a working-model-based cohort ZIM model
#'
#' This convenience wrapper calls [zim_cohort()] with `beta_method = "WKM"`.
#'
#' @param ... Arguments passed to [zim_cohort()].
#' @return A list containing fitted estimates and diagnostics.
#' @export
zim_cohort_wkm <- function(...) {
  .zim_impl_env("cohort")$zim_cohort_wkm(...)
}

#' Fit a Zero-Augmented Integration Method model for case-control data
#'
#' This wrapper calls the case-control implementation bundled with the
#' package. Supply the arguments described in the package README and the
#' example scripts.
#'
#' @param ... Arguments passed to the case-control ZIM implementation.
#' @return A list containing fitted estimates and diagnostics.
#' @export
zim_case_control <- function(...) {
  .zim_impl_env("case_control")$zim_case_control(...)
}

#' Fit a working-model-based case-control ZIM model
#'
#' This convenience wrapper calls [zim_case_control()] with
#' `beta_method = "WKM"`.
#'
#' @param ... Arguments passed to [zim_case_control()].
#' @return A list containing fitted estimates and diagnostics.
#' @export
zim_case_control_wkm <- function(...) {
  .zim_impl_env("case_control")$zim_case_control_wkm(...)
}

#' Fit a cross-fitted case-control ZIM model
#'
#' This wrapper calls the cross-fitted case-control implementation bundled with
#' the package.
#'
#' @param ... Arguments passed to the cross-fitted case-control implementation.
#' @return A list containing fitted estimates and diagnostics.
#' @export
zim_case_control_crossfit <- function(...) {
  .zim_impl_env("case_control")$zim_case_control_crossfit(...)
}

#' Locate a bundled ZIM example
#'
#' @param ... Path components below the installed `examples` directory.
#' @param mustWork Whether to signal an error when the requested file is not
#'   present.
#' @return A normalized path to the requested example file.
#' @export
zim_example_path <- function(..., mustWork = TRUE) {
  path <- system.file("examples", ..., package = "ZIM")
  if (mustWork && !nzchar(path)) {
    stop("The requested bundled example was not found.", call. = FALSE)
  }
  if (nzchar(path)) normalizePath(path, mustWork = mustWork) else path
}
