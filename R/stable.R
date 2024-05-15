#' Wrapper to libstable4u
#'
#' @param x argument of PDF
#' @param pars parameters of stable distribution, see libstable4u
#'
#' @description {A wrapper to \code{libstable4u} functions to match
#' naming conventions of basic pdfs/cdfs/generators.}
#' @return numeric
#' @export dstable
dstable <- function(x, pars)
{
  return(libstable4u::stable_pdf(x, pars))
}

#' Wrapper to libstable4u
#'
#' @param n number of variates to simulate
#' @param pars parameters of stable distribution, see libstable4u
#'
#' @description {A wrapper to \code{libstable4u} functions to match
#' naming conventions of basic pdfs/cdfs/generators.}
#' @return numeric
#' @export rstable
rstable <- function(n, pars)
{
  return(libstable4u::stable_rnd(n, pars))
}

#' Wrapper to libstable4u
#'
#' @param x argument of CDF
#' @param pars parameters of stable distribution, see libstable4u
#'
#' @description {A wrapper to \code{libstable4u} functions to match
#' naming conventions of basic pdfs/cdfs/generators.}
#' @return numeric
pstable <- function(x, pars)
{
  return(libstable4u::stable_cdf(x, pars))
}

#' Wrapper to libstable4u
#'
#' @param p argument of inverse CDF
#' @param pars parameters of stable distribution, see libstable4u
#'
#' @description {A wrapper to \code{libstable4u} functions to match
#' naming conventions of basic pdfs/cdfs/generators.}
#' @return numeric
#' @export qstable
qstable <- function(p, pars)
{
  return(libstable4u::stable_q(p, pars))
}

#' Wrapper to libstable4u
#'
#' @param x argument of inverse CDF
#' @param pars parameters of stable distribution, see libstable4u
#'
#' @description {A wrapper to \code{libstable4u} functions to match
#' naming conventions of basic pdfs/cdfs/generators.}
#' @return numeric
#' @export pstable
pstable <- function(x, pars)
{
  return(libstable4u::stable_cdf(x, pars))
}

#' Value at Risk for stable distribution
#'
#' @param p level of confidence
#' @param pars parameters of stable distribution, see libstable4u
#'
#' @description {Computes the inverse CDF at \code{1-p}.}
#' @return numeric
#' @export stableVAR
stableVAR <- function(p, pars)
{

  loss <- qstable(1-p, pars)
  return(loss)
}

