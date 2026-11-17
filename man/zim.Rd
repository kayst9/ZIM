\name{zim}
\alias{zim_cohort}
\alias{zim_cohort_wkm}
\alias{zim_case_control}
\alias{zim_case_control_wkm}
\alias{zim_case_control_crossfit}
\title{Fit Zero-Augmented Integration Method Models}
\usage{
zim_cohort(...)

zim_cohort_wkm(...)

zim_case_control(...)

zim_case_control_wkm(...)

zim_case_control_crossfit(...)
}
\arguments{
  \item{...}{Arguments passed to the corresponding ZIM implementation. The
  installed \code{README.md} and bundled examples describe the required
  internal data, target model, working models, external coefficient estimates,
  and study-specific sample-size information.}
}
\value{A list containing the ZIM estimate, its estimated standard error, the
internal estimate, covariance information, and diagnostics.}
\description{Fits ZIM models for cohort or case-control data using
individual-level internal data and external regression coefficient estimates.}
\details{The cohort and case-control implementations use separate private
environments because they have distinct internal helper functions. The
working-model-based convenience functions set \code{beta_method = "WKM"}.
The case-control cross-fitting function is provided for the corresponding
case-control workflow.}
\examples{
\dontrun{
library(ZIM)
?zim_cohort
?zim_case_control
}
}
