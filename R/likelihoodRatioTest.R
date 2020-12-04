#' Likelihood ratio hypothesis test
#'
#' @param f denominator density, true under H1
#' @param g numerator density, the alternative choice, true under H2
#' @param alpha1 error of choosing H2 when H1 is true
#' @param alpha2 error of choosing H1 when H2 is true
#'
#' @description {Hypothesis test for choosing among
#' two densities for fitting data.}
#'
#' @return list
#' @export likelihood_ratio
likelihood_ratio <- function(f, g, alpha1 = 10^-3, alpha2 = 10^-3)
{
  z <- prod(g/f)
  a <- alpha2/(1-alpha1)
  b <- (1-alpha2)/alpha1
  if(z <= a)
  {
    return(list(z = z, interval = c(a, b), result = "Choose H1: X~f"))
  } else if(z >= b)
  {
    return(list(z = z, interval = c(a, b), result = "Choose H2: X~g"))
  } else
  {
    return(list(z = NA, interval = NA, result = "Test inconclusive; need more samples"))
  }
}
