#' The PDF of the hitting time of for a GBM to reach a price level
#'
#' @param t time input
#' @param gbm the parameters defining the GBM, a drift and volatility
#' @param target the target price to reach in \code{t} time
#' @param spot the spot price
#'
#' @description {The well-known density for the hitting time of a given level
#' in a GBM starting from a known price.}
#' @return numeric or vector
#' @export dhit_gbm
dhit_gbm <- function(t, gbm, target, spot)
{

  mu <- gbm[1]
  volat <- gbm[2]
  a <- log(target/spot)
  b <- (mu-volat^2/2)
  C <- volat
  f <- (a/(sqrt(2*pi*(C^2)*(t^3))))*exp(-(a-b*t)^2/(2*(C^2)*t))
  return(f)
}


#' The CDF of the hitting time of for a GBM to reach a price level
#'
#' @param t time input
#' @param gbm the parameters defining the GBM, a drift and volatility
#' @param target the target price to reach in \code{t} time
#' @param spot the spot price
#'
#' @description {The well-known CDF for the hitting time of a given level
#' in a GBM starting from a known price.}
#' @return numeric or vector
#' @export phit_gbm
phit_gbm <- function(t, gbm, target, spot)
{
  mu <- gbm[1]
  volat <- gbm[2]
  a <- log(target/spot)
  b <- (mu-volat^2/2)
  C <- volat
  z <- (a-b*t)/(C*sqrt(t))
  z2 <- (-a-b*t)/(C*sqrt(t))
  f <- 1-stats::pnorm(z) + stats::pnorm(z2)*exp(2*a*b/(C^2))
  return(f)
}

#' The mode of the hitting time to a price level for a GBM
#'
#' @param gbm the parameters defining the GBM, drift and volatility
#' @param target the target price
#' @param spot the initial price
#'
#' @description {The mode of the hitting time density for a GBM to
#' reach a given level.}
#' @return numeric
#' @export hitmode_gbm
hitmode_gbm <- function(gbm, target, spot)
{
  mu <- gbm[1]
  volat <- gbm[2]
  a <- log(target/spot)
  b <- (mu-0.5*volat^2)
  C <- volat
  f <- (sqrt(4*(a*b)^2+9*C^4)-3*C^2)/(2*b^2)
  return(f)
}
