#' Sampling algorithm for Gaussian mixtures
#'
#' @param n number of variates to simulate
#' @param p the probability vector for component chances
#' @param mus the vector of component means
#' @param sigmas the vector of volatilities.
#' @description {The arguments can be the component vectors of a mixture model}
#' @return vector
#' @import mclust
#' @export rgmm
rgmm <- function(n, p, mus, sigmas)
{
  modelName <- "V"
  parameters <- list(pro = p, mean = mus, variance = list(sigmasq = sigmas^2))
  mclust::sim(modelName, parameters, n = n)[,2]
}

#' Probability density function for Gaussian mixtures
#'
#' @param x the real variable argument of the PDF
#' @param p component probabilities
#' @param mus component means
#' @param sigmas component standard deviations
#' @description {A mixture model of PDFs is a probabilistic weighted sum of individual PDFs}
#' @return numeric
dgmm1 <- function(x, p, mus, sigmas)
{
  sum(p*stats::dnorm(x, mean = mus, sd = sigmas))
}
#' Probability density function for Gaussian mixtures
#'
#' @param x the real variable argument of the PDF
#' @param p component probabilities
#' @param mus component means
#' @param sigmas component standard deviations
#' @description {A mixture model of PDFs is a probabilistic weighted sum of individual PDFs}
#' @return numeric
#' @export dgmm
dgmm <- Vectorize(dgmm1, vectorize.args = "x")


#' Cumulative distribution function for Gaussian mixtures
#'
#' @param x the real variable argument of the PDF
#' @param p component probabilities
#' @param mus component means
#' @param sigmas component standard deviations
#' @description {A mixture model of PDFs is a probabilistic weighted sum of individual PDFs}
#' @return numeric
pgmm1 <- function(x, p, mus, sigmas)
{
  sum(p*stats::pnorm(x, mean = mus, sd = sigmas))
}
#' Cumulative distribution function for Gaussian mixtures
#'
#' @param x the real variable argument of the PDF
#' @param p component probabilities
#' @param mus component means
#' @param sigmas component standard deviations
#' @description {A mixture model of PDFs is a probabilistic weighted sum of individual PDFs}
#' @return numeric
#' @export pgmm
pgmm <- Vectorize(pgmm1, vectorize.args = "x")

#' Convenience to extract parameters
#'
#' @param mcfit object returned from \code{Mclust}
#' @param scale to scale returns (use only for log-returns data)
#' @param continuous boolean for discrete time or continuoust time
#'
#' @return matrix
#' @export extract_mixture
extract_mixture <- function(mcfit, scale = 252, continuous = FALSE)
{
  # Extract parameters
  probs <- mcfit$parameters$pro
  if(!continuous)
  {
    mus <- mcfit$parameters$mean
    sigmas <- sqrt(mcfit$parameters$variance$sigmasq)
  } else if(continuous)
  {
    sigmas <- sqrt(mcfit$parameters$variance$sigmasq*scale)
    mus <- mcfit$parameters$mean*scale+0.5*sigmas^2
  }

  parameters <- rbind(probs, mus, sigmas)
  return(parameters)
}

#' Fit log-normal mixture diffusion to daily-time series of log-returns
#'
#' @param log_returns multi-dimensional time-series of daily log-returns
#' @param timeScale time-scale to convert by.
#'
#' @description {wrapper to \code{mclust} EM algorithm for the parameters of a Gaussian-mixture diffusion.}
#'
#' @return matrix of component parameters, three rows: probabilities, drifts, volatilities
#' @export fitMixtureDiffusion
fitMixtureDiffusion <- function(log_returns, timeScale = 1/252)
{
  if(is.null(log_returns))
  {
    stop("Must pass 'log_returns' as an argument")
  }
  if(ncol(log_returns) > 1)
  {
    stop("'fitMixtureDiffusion' is for univariate log-returns only. Use 'fitGBMs' instead.")
  }
  mcfit <- mclust::Mclust(data = log_returns)
  mix_param <- extract_mixture(mcfit, 1/timeScale, TRUE)
  return(mix_param)
}


#' The diffusion coefficient function for a log-normal mixture diffusion
#'
#' @param t the current time
#' @param x the current log-price
#' @param param the matrix of component parameters by row, see details.
#' @param truncation the truncation of time
#'
#' @description {The mixture diffusion coefficient is a weighted average of variances.}
#' @details {The argument \code{param} must be a square matrix with rows representing,
#' in order, the probabilities, drifts, and volatilities of each component of the mixture.}
#' @return numeric matrix
#' @export mixture_vol
mixture_vol <- function(t, x, param, truncation = 0.5/252)
{
  probs <- param[1, ]
  mus <- param[2, ]
  vols <- param[3, ]

  n <- length(t)
  m <- length(x)
  u <- matrix(0, nrow = n, ncol = m)
  if(n > 1)
  {
    t[which(t == 0)] <- truncation
  } else
  {
    if(t == 0)
    {
      t <- truncation
    }
  }

  for(i in 1:n)
  {
    for(j in 1:m)
    {
      p_i <- findistr::dgbm(exp(x[j]), t[i], 1, mus, vols)
      p <- sum(probs*p_i)
      numerator <- sum(probs*(vols^2)*p_i)
      u[i, j] <- sqrt(numerator/p)
    }
  }
  return(u)
}

#' The drift coefficient function for a log-normal mixture diffusion
#'
#' @param t the current time
#' @param x the current log-price
#' @param param the matrix of component parameters by row, see details.
#' @param truncation the truncation of time
#'
#' @description {The mixture drift coefficient is a weighted average of variances.}
#' @details {The argument \code{param} must be a square matrix with rows representing,
#' in order, the probabilities, drifts, and volatilities of each component of the mixture.}
#' @return numeric matrix
#' @export mixture_drift
mixture_drift <- function(t, x, param, truncation = 0.5/252)
{
  probs <- param[1, ]
  mus <- param[2, ]
  vols <- param[3, ]
  n <- length(t)
  m <- length(x)
  u <- matrix(0, nrow = n, ncol = m)
  if(n > 1)
  {
    t[which(t == 0)] <- truncation
  } else
  {
    if(t == 0)
    {
      t <- truncation
    }
  }
  for(i in 1:n)
  {
    for(j in 1:m)
    {
      p_i <- findistr::dgbm(exp(x[j]), t[i], 1, mus, vols)
      p <- sum(probs*p_i)
      numerator <- sum(probs*mus*p_i)
      u[i, j] <- numerator/p
    }
  }
  return(u)
}

