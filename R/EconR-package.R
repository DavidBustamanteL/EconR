###### 1. EconR Package Documentation ####

#' EconR: Econometric Models with Consistent Estimation Semantics
#'
#' @description
#' EconR provides a unified framework for estimating econometric models
#' in R, with consistent interfaces for estimation, covariance estimation,
#' statistical inference, and model reporting.
#'
#' @details
#' The package currently supports:
#'
#' \itemize{
#'   \item \code{econ_lm()}: Ordinary least squares (OLS) and
#'         weighted least squares (WLS).
#'   \item \code{econ_binary()}: Binary logit and probit models.
#'   \item \code{econ_ordered()}: Ordered logit and probit models.
#' }
#'
#' EconR supports classical, heteroskedasticity-robust, and
#' cluster-robust covariance estimators, depending on the model.
#'
#' The package integrates with standard R methods, including
#' \code{summary()}, \code{coef()}, \code{vcov()},
#' \code{predict()}, \code{confint()}, and \code{logLik()}.
#'
#' It also provides compatibility with packages such as
#' \code{broom}, \code{sandwich}, and \code{modelsummary}.
#'
#' @section Model Estimation:
#' \describe{
#'   \item{Linear models}{\code{\link{econ_lm}}}
#'   \item{Binary response models}{\code{\link{econ_binary}}}
#'   \item{Ordered response models}{\code{\link{econ_ordered}}}
#' }
#'
#' @section Goodness of Fit:
#' Binary and ordered response models provide McFadden's
#' pseudo-R-squared alongside standard likelihood-based
#' goodness-of-fit statistics.
#'
#' @author
#' David Bustamante Lazo
#'
#' @seealso
#' \code{\link{econ_lm}},
#' \code{\link{econ_binary}},
#' \code{\link{econ_ordered}}
#'
#' @docType package
#' @name EconR-package
#' @aliases EconR
NULL