% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbm-family.R
\name{mseGBM}
\alias{mseGBM}
\title{Mean square error between empirical PDF and time-dependent Gaussian}
\usage{
mseGBM(theta, epdf, h, idynamics)
}
\arguments{
\item{theta}{parameter to optimize}

\item{epdf}{object returned from \code{density()}, the empirical probability density function.}

\item{h}{time-scale}

\item{idynamics}{the integrated mean and volatility functions, as a list}
}
\value{
numeric
}
\description{
{The mean-square-error between the empirical PDF of daily log returns against a time-dependent
GBM model. Used for numerical calibration to the empirical density.}
}
\details{
{Here the coefficients in the SDE \eqn{dS_t=\mu(t) S_t dt +\sigma(t) S_t dB_t} are assumed
to be functions of time.}
}
