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




