#' The probability density function under Merton's jump-diffusion model for log-price dynamics
#'
#' @param x input value for the CDF
#' @param t the time input, for daily use 1/252
#' @param param a vector of parameters defining Merton's jump dynamics. See details.
#' @description {The PDF of the log-returns in a time interval under Merton's jump diffusion model. Includes
#' compensator in the drift.}
#' @details {The argument \code{param} must be a vector of five real numbers, \eqn{(\mu, \sigma, \lambda, \alpha, \beta)},
#' representing the mean drift, volatility, mean-rate of jumps, and the mean
#' log-jump size and the standard deviation of the log-jump size, in this order. The volatility,
#' mean rate of jumps, and standard deviation of log-jump size must all be positive.}
#' @return numeric
dmerton1 <- function(x, t, param)
{
  if(param[2] <= 0)
  {
    stop("Volatility must be positive.")
  }
  if(param[3] <= 0)
  {
    stop("The mean rate of jumps must be positive.")
  }
  if(param[5] <= 0)
  {
    stop("The standard deviation of log-jump sizes must be positive.")
  }
  drift <- param[1]
  volat <- param[2]
  lambda <- param[3]
  a <- param[4]
  b <- param[5]
  # Generate enough summands
  N <- stats::qpois(0.99, lambda = lambda*t)
  n <- 0:N
  p <- stats::dpois(n, lambda = lambda*t)
  eta1 <- exp(a+0.5*b^2)-1
  phi <- stats::dnorm(x, (drift-0.5*volat^2-lambda*eta1)*t+n*a, sqrt(t*volat^2+n*b^2))
  return(sum(p*phi))
}

#' The probability density function under Merton's jump-diffusion model for log-price dynamics
#'
#' @param x input value for the CDF
#' @param t the time input, for daily use 1/252
#' @param param a vector of parameters defining Merton's jump dynamics. See details.
#' @description {The PDF of the log-returns in a time interval under Merton's jump diffusion model. Includes
#' compensator in the drift.}
#' @details {The argument \code{param} must be a vector of five real numbers, \eqn{(\mu, \sigma, \lambda, \alpha, \beta)},
#' representing the mean drift, volatility, mean-rate of jumps, and the mean
#' log-jump size and the standard deviation of the log-jump size, in this order. The volatility,
#' mean rate of jumps, and standard deviation of log-jump size must all be positive.}
#' @return numeric or vector depending on the length of the input \code{x}.
#' @export dmerton
dmerton <- function(x, t, param)
{
  if(length(x)==1)
  {
    return(dmerton1(x, t, param))
  } else if(length(x) > 1)
  {
    v <- matrix(0, nrow = length(x))
    for(i in 1:length(x))
    {
      v[i] <- dmerton1(x[i], t, param)
    }
    # mapply(function(X){
    #   dmerton1(X, t, param)
    # }, X = x)
    return(v)
  } else if(length(x) == 0)
  {
    stop("Bad input for x")
  }
}

#' The cumulative distribution function under Merton's jump-diffusion model for log-price dynamics
#'
#' @param x input value for the CDF
#' @param t the time input, for daily use 1/252
#' @param param a vector of parameters defining Merton's jump dynamics. See details.
#' @description {The PDF of the log-returns in a time interval under Merton's jump diffusion model. Includes
#' compensator in the drift.}
#' @details {The argument \code{param} must be a vector of five real numbers, \eqn{(\mu, \sigma, \lambda, \alpha, \beta)},
#' representing the mean drift, volatility, mean-rate of jumps, and the mean
#' log-jump size and the standard deviation of the log-jump size, in this order. The volatility,
#' mean rate of jumps, and standard deviation of log-jump size must all be positive.}
#' @return numeric or vector depending on the length of the input \code{x}.
#' @export pmerton
pmerton <- function(x, t, param)
{
  if(length(x) == 1)
  {
    w <- stats::integrate(dmerton, lower = -Inf, upper = x, t = t, param = param)
    return(w$value)
  } else if(length(x) > 1)
  {
    w <- matrix(0, nrow = length(x))
    for(i in 1:length(x))
    {
      w[i] <- stats::integrate(dmerton, lower = -Inf, upper = x[i], t = t, param = param)$value

    }
    return(w)
  } else
  {
    stop("bad input for the CDF")
  }
}

#' Simulate the terminal RV of a log-Merton jump diffusion
#'
#' @param n number of variates to simulate
#' @param t the terminal time-horizon
#' @param param the parameter set defining the jump diffusion dynamics
#'
#' @description {Simulates the terminal value of a log jump diffusion
#' under Merton dynamics. This should not be used to simulate sample paths over time,
#' but only the terminal log variate.}
#'
#' @return numeric or vector
#' @export rmerton
rmerton <- function(n, t, param)
{
  mu <- param[1]
  volat <- param[2]
  lambda <- param[3]
  alpha <- param[4]
  beta <- param[5]
  eta <- exp(alpha+0.5*beta^2)-1
  variates <- matrix(0, nrow = n)
  drift <- (mu-0.5*volat^2-lambda*eta)
  jumpCounts <- stats::rpois(n, lambda = lambda*t)
  z <- stats::rnorm(n, mean = 0, sd = 1)
  for(i in 1:n)
  {
    jumps <- sum(stats::rnorm(jumpCounts[i], alpha, beta))
    variates[i] <- drift*t+volat*sqrt(t)*z[i]+jumps
  }
  return(variates)
}


