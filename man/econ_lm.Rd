\name{econ_lm}
\alias{econ_lm}
\title{Fit an OLS or WLS model}

\usage{
econ_lm(
  formula,
  data,
  weights = NULL,
  subset = NULL,
  na.action = stats::na.omit,
  vcov = c(
    "classical",
    "HC0",
    "HC1",
    "HC2",
    "HC3",
    "CR0",
    "CR2",
    "stata"
  ),
  cluster = NULL
)
}

\arguments{
\item{formula}{Model formula.}
\item{data}{Data frame.}
\item{weights}{Optional non-negative observation weights.}
\item{subset}{Optional subset.}
\item{na.action}{Missing-data action.}
\item{vcov}{Covariance estimator.}
\item{cluster}{Cluster IDs or column name.}
}

\description{
Fits ordinary or weighted least squares using model-specific missing-data handling.
}