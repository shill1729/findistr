#' Cumulative distribution function for Black-Scholes distribution
#'
#' @param s real-numeric input to PDF representing price level
#' @param t positive time value
#' @param spot the initial spot price
#' @param rate the risk-free rate of return
#' @param volat the volatility level
#'
#' @description {The cumulative distribution function for Black-Scholes distributed
#' variables. }
#' @return numeric
#' @export pblackscholes
pblackscholes <- function(s, t, spot, rate, volat)
{
  z <- stats::pnorm(log(s/spot), (rate-0.5*volat^2)*t, sqrt(t)*volat)
  return(z)
}

#' #' Probability density function for Black-Scholes distribution
#'
#' @param s real-numeric input to PDF representing price level
#' @param t positive time value
#' @param spot the initial spot price
#' @param rate the risk-free rate of return
#' @param volat the volatility level
#'
#' @description {The probability density function for Black-Scholes distributed
#' variables.}
#' @return numeric
#' @export dblackscholes
dblackscholes <- function(s, t, spot, rate, volat)
{
  z <- stats::dnorm(log(s/spot), (rate-0.5*volat^2)*t, sqrt(t)*volat)/s
  return(z)
}

#' Generate variates under Black-Scholes distribution
#'
#' @param n number of variates to generate
#' @param t positive time value
#' @param spot the initial spot price
#' @param rate the risk-free rate of return
#' @param volat the volatility level
#'
#' @description {Generate random variates with Black-Scholes distribution.}
#' @return numeric
#' @export rblackscholes
rblackscholes <- function(n, t, spot, rate, volat)
{
  x <- stats::rnorm(n, (rate-0.5*volat^2)*t, sqrt(t)*volat)
  return(spot*exp(x))
}
