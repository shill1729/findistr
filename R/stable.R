#' Wrapper to libstableR
#'
#' @param x argument of PDF
#' @param pars parameters of stable distribution, see libstableR
#'
#' @description {A wrapper to \code{libstableR} functions to match
#' naming conventions of basic pdfs/cdfs/generators.}
#' @return numeric
#' @export dstable
dstable <- function(x, pars)
{
  return(libstableR::stable_pdf(x, pars))
}

#' Wrapper to libstableR
#'
#' @param n number of variates to simulate
#' @param pars parameters of stable distribution, see libstableR
#'
#' @description {A wrapper to \code{libstableR} functions to match
#' naming conventions of basic pdfs/cdfs/generators.}
#' @return numeric
#' @export rstable
rstable <- function(n, pars)
{
  return(libstableR::stable_rnd(n, pars))
}

#' Wrapper to libstableR
#'
#' @param x argument of CDF
#' @param pars parameters of stable distribution, see libstableR
#'
#' @description {A wrapper to \code{libstableR} functions to match
#' naming conventions of basic pdfs/cdfs/generators.}
#' @return numeric
pstable <- function(x, pars)
{
  return(libstableR::stable_cdf(x, pars))
}

#' Wrapper to libstableR
#'
#' @param p argument of inverse CDF
#' @param pars parameters of stable distribution, see libstableR
#'
#' @description {A wrapper to \code{libstableR} functions to match
#' naming conventions of basic pdfs/cdfs/generators.}
#' @return numeric
#' @export qstable
qstable <- function(p, pars)
{
  return(libstableR::stable_q(p, pars))
}

#' Wrapper to libstableR
#'
#' @param x argument of inverse CDF
#' @param pars parameters of stable distribution, see libstableR
#'
#' @description {A wrapper to \code{libstableR} functions to match
#' naming conventions of basic pdfs/cdfs/generators.}
#' @return numeric
#' @export pstable
pstable <- function(x, pars)
{
  return(libstableR::stable_cdf(x, pars))
}

#' Value at Risk for stable distribution
#'
#' @param p level of confidence
#' @param pars parameters of stable distribution, see libstableR
#'
#' @description {Computes the inverse CDF at \code{1-p}.}
#' @return numeric
#' @export stableVAR
stableVAR <- function(p, pars)
{

  loss <- qstable(1-p, pars)
  return(loss)
}

