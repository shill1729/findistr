% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbm-fitting.R
\name{fitGBM}
\alias{fitGBM}
\title{Fit a geometric Brownian motion to a daily log-returns time-series}
\usage{
fitGBM(log_returns, timeScale = 1/252)
}
\arguments{
\item{log_returns}{uni-variate time-series of daily log-returns}

\item{timeScale}{time-scale to convert by.}
}
\value{
data.frame containing \code{drift} and \code{volat} point-estimates.
}
\description{
{MLE estimates for the parameters of a GBM, the mean
drift rate and the volatility coefficient.}
}
