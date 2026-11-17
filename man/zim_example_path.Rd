\name{zim_example_path}
\alias{zim_example_path}
\title{Locate Bundled ZIM Examples}
\usage{
zim_example_path(..., mustWork = TRUE)
}
\arguments{
  \item{...}{Path components below the installed \code{examples} directory.}
  \item{mustWork}{Whether to signal an error if the requested file is absent.}
}
\value{A normalized path to a bundled example file.}
\description{Returns the installed path to an example script or supporting
file bundled with \pkg{ZIM}.}
\examples{
zim_example_path("cohort_one_population.R")
}
