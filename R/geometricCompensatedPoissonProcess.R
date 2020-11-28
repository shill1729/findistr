#' Cumulative distribution function for geometric compensated Poisson processes
#'
#' @param s price argument to CDF
#' @param t time argument to CDF
#' @param spot initial spot price
#' @param a the coefficient of the Poisson RV
#' @param b the coefficient of the compensated drift
#' @param lambda the mean rate of arrivals
#'
#' @description {The CDF of a geometrically compensated Poisson process at a given
#' time and initial spot price.}
#' @return numeric
#' @export pgcpp
pgcpp <- function(s, t, spot, a, b, lambda)
{
  z <- stats::ppois((log(s/spot)+b*t)/a, lambda*t)
  return(z)
}

#' Simulate variates under geometric compensated Poisson distributions at a given time.
#'
#' @param n number of variates to simulate
#' @param t time argument to CDF
#' @param spot initial spot price
#' @param a the coefficient of the Poisson RV
#' @param b the coefficient of the compensated drift
#' @param lambda the mean rate of arrivals
#'
#' @description {Simulate variates that follow a geometric compensated Poisson distribution
#' of a given time and initial spot price.}
#' @return numeric
#' @export rgcpp
rgcpp <- function(n, t, spot, a, b, lambda)
{
  x <- a*stats::rpois(n, lambda*t)-b*t
  return(spot*exp(x))
}
