#' Fit a geometric Brownian motion to a daily log-returns time-series
#'
#' @param log_returns uni-variate time-series of daily log-returns
#' @param timeScale time-scale to convert by.
#'
#' @description {MLE estimates for the parameters of a GBM, the mean
#' drift rate and the volatility coefficient.}
#'
#' @return data.frame containing \code{drift} and \code{volat} point-estimates.
#' @export fitGBM
fitGBM <- function(log_returns, timeScale = 1/252)
{
  if(is.null(log_returns))
  {
    stop("Must pass 'log_returns' as an argument")
  }
  if(ncol(log_returns) > 1)
  {
    stop("'fitGBM' is for univariate log-returns only. Use 'fitGBMs' instead.")
  }
  volat <- stats::sd(log_returns)/sqrt(timeScale)
  drift <- mean(log_returns)/timeScale+0.5*volat^2
  param <- c(drift = drift, volat = volat)
  return(param)
}

#' Fit correlated geometric Brownian motions to a matrix of daily log-returns time-series
#'
#' @param log_returns multi-dimensional time-series of daily log-returns
#' @param timeScale time-scale to convert by.
#'
#' @description {MLE estimates for the parameters of a vector of correlated GBMs, the mean
#' drift rate and the volatility coefficient.}
#'
#' @return list containing the vector \code{drift} and matrix \code{Sigma} of point-estimates.
#' @export fitGBMs
fitGBMs <- function(log_returns, timeScale = 1/252)
{
  if(is.null(log_returns))
  {
    stop("Must pass 'log_returns' as an argument")
  }
  if(ncol(log_returns) == 1)
  {
    stop("Use 'fitGBM' for univariate log-returns")
  }
  Sigma <- stats::cov(log_returns)/timeScale
  drift <- apply(log_returns, 2, mean)/timeScale+0.5*diag(Sigma)
  return(list(drift = drift, Sigma = Sigma))
}


#' Exponentially moving average
#'
#' @param x data-set
#' @param lambda weight parameter (near 0 corresponds to sample mean)
#' @param h timescale
#'
#' @description {Exponentially weighted sample average. Past data has exponentially
#' decaying weights.
#' }
#' @return numeric
#' @export ema
ema <- function(x, lambda = 0.94, h = 1/252)
{
  n <- dim(x)[1]
  weights <- (1-lambda)^(0:(n-1))
  weights <- weights*(lambda/(1-(1-lambda)^n))
  mu <- rev(t(weights))%*%x
  return(mu/h)
}

#' Exponentially moving covariance
#'
#' @param X data-set
#' @param lambda weight parameter (near 0 corresponds to sample mean)
#' @param h timescale
#'
#' @description {Exponentially weighted sample covariance Past data has exponentially
#' decaying weights.
#' }
#' @return matrix
#' @export ewmc
ewmc <- function(X, lambda = 0.94, h = 1/252)
{
  N <- nrow(X)
  # Center data
  a1 <- X-apply(X, 2, mean)
  # Compute weights
  ws <- (1-lambda)^(0:(N-1))
  ws <- ws*(lambda/(1-(1-lambda)^N))
  # Compute weighted sample covariance in matrix form
  Sigma <- t(rev(ws)*a1)%*%a1
  # Tidy and return
  colnames(Sigma) <- colnames(X)
  rownames(Sigma) <- colnames(X)
  return(Sigma/h)
}

#' Fit GBM model where parameters are estimated using EMAs
#'
#' @param X data-set
#' @param lambda weight parameter (near 0 corresponds to sample mean)
#' @param h timescale
#'
#' @description {Assuming log-returns. This estimates drift vector and covariance
#' matrix but with an EMA filter. A log-adjustment is made to the drift per
#' the GBM dynamics.
#' }
#' @return list of drift and covariance matrix Sigma
#' @export fit_ema_gbm
fit_ema_gbm <- function(X, lambda = 0.94, h = 1/252)
{
  if(is.null(X))
  {
    stop("Must pass 'log_returns' as an argument")
  }
  if(ncol(X) == 1)
  {
    Sigma <- as.numeric(ewmc(X, lambda, h))
    drift <- as.numeric(ema(X, lambda, h)+0.5*diag(Sigma))
    return(list(drift = drift, Sigma = Sigma))
  }
  Sigma <- ewmc(X, lambda, h)
  drift <- ema(X, lambda, h)+0.5*diag(Sigma)
  return(list(drift = drift, Sigma = Sigma))
}

