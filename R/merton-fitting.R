#' Log-likelihood for Merton jump diffusion log-increments
#'
#' @param param a vector of parameters defining Merton's jump dynamics. See details.
#' @param x input value for the CDF
#' @param t the time input, for daily use 1/252
#' @description {The log of the PDF of the log-returns in a time interval under Merton's jump diffusion model. Includes
#' compensator in the drift.}
#' @details {The argument \code{param} must be a vector of five real numbers, \eqn{(\mu, \sigma, \lambda, \alpha, \beta)},
#' representing the mean drift, volatility, mean-rate of jumps, and the mean
#' log-jump size and the standard deviation of the log-jump size, in this order. The volatility,
#' mean rate of jumps, and standard deviation of log-jump size must all be positive.}
#' @return numeric or vector depending on the length of the input \code{x}.
#' @export logLikeMerton
logLikeMerton <- function(param, x, t)
{
  f <- dmerton(x, t, param)
  return(sum(log(f[f>0])))
}

#' MSE between empirical and model PDF under Merton's jump-diffusion
#'
#' @param param a vector of parameters defining Merton's jump dynamics. See details.
#' @param epdf the empirical pdf returned from \code{density}
#' @param t the time input, for daily use 1/252
#' @description {MSE between empirical model PDF under Merton's jump-diffusion
#' for a given time-series of log-returns.}
#' @details {The argument \code{param} must be a vector of five real numbers, \eqn{(\mu, \sigma, \lambda, \alpha, \beta)},
#' representing the mean drift, volatility, mean-rate of jumps, and the mean
#' log-jump size and the standard deviation of the log-jump size, in this order. The volatility,
#' mean rate of jumps, and standard deviation of log-jump size must all be positive.}
#' @return numeric or vector depending on the length of the input \code{x}.
#' @export mseMerton
mseMerton <- function(param, epdf, t)
{
  mpdf <- dmerton(epdf$x, t, param)
  m <- mean((mpdf-epdf$y)^2)
  return(m)
}

#' Fit the symmetric Merton jump diffusion model to log-return data.
#'
#' @param x the data-set of daily log returns
#' @param k time scale, defaults to 1/252 for daily data
#' @param numDev number of daily-standard deviations to use in initial guess for jump sizes.
#' @return output from optim
#' @description {Use base R's \code{optim} to maximize the log-likelihood in
#' the Merton jump diffusion model with symmetric jumps.}
#' @export fitSymmetricMerton
fitSymmetricMerton <- function(x, k = 1/252, numDev = 1)
{
  volat <- stats::sd(x[abs(x) < numDev*stats::sd(x)])/sqrt(k)
  mu <- mean(x[abs(x) < numDev*stats::sd(x)])/k+0.5*volat^2
  beta <- stats::sd(x[abs(x) >= numDev*stats::sd(x)])
  lambda <- (sum(abs(x>=numDev*stats::sd(x)))/length(x))/k
  initialGuess <- c(mu, volat, lambda, beta)
  print(initialGuess)
  logLikeSymMert <- function(param, x, t)
  {
    p <- c(param[1:3], 0, param[4])
    logLikeMerton(p, x, t)
  }

  # MLE numerical rroutine
  mle <- stats::optim(par = initialGuess,
                      fn = logLikeSymMert, x = x, t = k,
                      method = "L-BFGS-B",
                      lower = c(-Inf, 0.01, 0.01, 0.01),
                      control = list(fnscale = -1, trace = 6)
  )
  print(mle)
  return(mle)
}
