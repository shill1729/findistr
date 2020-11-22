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
    mus <- mcfit$parameters$mean*scale
    sigmas <- sqrt(mcfit$parameters$variance$sigmasq*scale)
  } else if(continuous)
  {
    sigmas <- sqrt(mcfit$parameters$variance$sigmasq*scale)
    mus <- mcfit$parameters$mean*scale+0.5*sigmas^2
  }

  parameters <- rbind(probs, mus, sigmas)
  return(parameters)
}
