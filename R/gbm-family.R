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
#' @export pgbm
pgbm <- function(s, t, spot, rate, volat)
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
#' @export dgbm
dgbm <- function(s, t, spot, rate, volat)
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
#' @export rgbm
rgbm <- function(n, t, spot, rate, volat)
{
  x <- stats::rnorm(n, (rate-0.5*volat^2)*t, sqrt(t)*volat)
  return(spot*exp(x))
}

#' Multivariate normal probability density function
#'
#' @param x vector of reals
#' @param mu vector means
#' @param Sigma covariance matrix
#'
#' @description {Everybody's favorite distribution, now in multiple dimensions!
#' (Why isn't this in base R anyway?)}
#' @return numeric
#' @export dmvtnorm
dmvtnorm <- function(x, mu, Sigma)
{
  k <- length(x)
  detSigma <- det(Sigma)
  Siginv <- solve(Sigma)
  z <- exp(-0.5*t(x-mu)%*%Siginv%*%(x-mu))/sqrt(detSigma*(2*pi)^k)
  return(as.numeric(z))
}

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

#' Generate a vector of two correlated standard Gaussian RVs
#'
#' @param n number of vectors to simulate
#' @param rho the correlation factor
#'
#' @description {Uses the Cholesky decomposition on the 2x2 matrix
#' to simulate a vector of two correlated standard normals.}
#' @return numeric vector/matrix
#' @export rcornorm
rcornorm <- function(n, rho)
{
  w <- matrix(0, n, ncol = 2)
  for(i in 1:n)
  {
    z <- stats::rnorm(2)
    ch <- rbind(c(1, 0),
                c(rho, sqrt(1-rho^2))
    )
    w[i, ] <- t(ch%*%z)
  }
  return(w)
}
