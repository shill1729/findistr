#' Log-likelihood function for GBM with time-dependent coefficient functions.
#'
#' @param theta parameter to optimize
#' @param x observations of realized variates
#' @param h time-scale
#' @param idynamics the integrated mean and volatility functions, as a list
#'
#' @description {The log-likelihood of the daily log returns under a time-dependent
#' GBM model. Used for numerical MLE routines.}
#' @details {Here the coefficients in the SDE \eqn{dS_t=\mu(t) S_t dt +\sigma(t) S_t dB_t} are assumed
#' to be functions of time.}
#' @return numeric
#' @export logLikeGBM
logLikeGBM <- function(theta, x, h, idynamics)
{
  M <- idynamics[[1]]
  V <- idynamics[[2]]
  n <- length(x)
  ll <- -n*log(sqrt(2*pi*V(h, theta)^2))-(1/(2*V(h, theta)^2))*sum((x-(M(h, theta)-0.5*V(h, theta)^2))^2)
  return(ll)
}

#' Mean square error between empirical PDF and time-dependent Gaussian
#'
#' @param theta parameter to optimize
#' @param epdf object returned from \code{density()}, the empirical probability density function.
#' @param h time-scale
#' @param idynamics the integrated mean and volatility functions, as a list
#'
#' @description {The mean-square-error between the empirical PDF of daily log returns against a time-dependent
#' GBM model. Used for numerical calibration to the empirical density.}
#' @details {Here the coefficients in the SDE \eqn{dS_t=\mu(t) S_t dt +\sigma(t) S_t dB_t} are assumed
#' to be functions of time.}
#' @return numeric
#' @export mseGBM
mseGBM <- function(theta, epdf, h, idynamics)
{
  M <- idynamics[[1]]
  V <- idynamics[[2]]
  # Compute model-density
  mf <- stats::dnorm(epdf$x, mean = M(h, theta)-0.5*V(h, theta)^2, sd = V(h, theta))
  z <- (mf-epdf$y)^2
  return(mean(z))
}
