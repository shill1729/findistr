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
