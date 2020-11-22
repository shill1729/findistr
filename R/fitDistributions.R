#' Fit returns distribution to daily-arithmetic returns
#'
#' @param x time series of daily arithmetic returns
#' @param distr distribution name: "norm" or "unif"
#'
#' @return vector
#' @import mclust
#' @export fitDTFM
fitDTFM <- function(x, distr = "norm")
{
  if(distr == "unif")
  {
    return(c(min = min(x), max = max(x)))
  } else if(distr == "norm")
  {
    return(c(mean = mean(x), sd = stats::sd(x)))
  } else if(distr == "gmm")
  {
    mcfit <- mclust::Mclust(data = x)
    mix_param <- extract_mixture(mcfit, 1)
    return(mix_param)
  } else if(distr == "stable")
  {
    # Fit stable distribution
    stable_param <- libstableR::stable_fit_mle(rnd = x)
    return(stable_param)
  }
}

#' Probability density function of fitted model
#'
#' @param r region of returns
#' @param distr distribution name: "norm" or "unif"
#' @param param parameters of distribution as a vector
#'
#' @return numeric
#' @export ddtfm
ddtfm <- function(r, distr = "norm", param)
{
  ddistr <- paste("d", distr, sep = "")
  ddistr <- get(ddistr)
  if(distr == "norm" || distr == "unif")
  {
    args <- c(list(r), as.list(param))
  } else if(distr == "gmm")
  {
    args <- list(c(r), param[1, ], param[2, ], param[3, ])
  } else if(distr == "stable")
  {
    args <- list(c(r), param)
  }

  return(do.call(ddistr, args))
}


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

