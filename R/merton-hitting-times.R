
#' The exponent used in the method of images for hitting time densities
#'
#' @param a the boundary to hit above
#' @param t the current time
#' @param param the dynamics of the jump-diffusion, a vector of 5, the drift, volatility
#' mean-jump rate, mean jump-size, sd of jump size
#'
#' @description {The exponent used in the method of images.}
#' @return numeric or vector
eta <- function(a, t, param)
{
  p <- param
  p[1] <- p[1]+2*a/t
  w <- findistr::dmerton(a, t, p)
  w <- w[is.finite(w)]
  ww <- findistr::dmerton(a, t, param)
  ww <- ww[is.finite(ww)]
  ee <- log(w)-log(ww)
  return(ee)
}

#' The transition density used and reflect in the method of images
#'
#' @param x the space variable
#' @param t the current time
#' @param a the boundary to hit above
#' @param param the dynamics of the jump-diffusion, a vector of 5, the drift, volatility
#' mean-jump rate, mean jump-size, sd of jump size
#'
#' @description {The transition density figuring in the method of images.}
#' @return numeric or vector
trans_den_merton <- function(x, t, a, param)
{
  p <- param
  p[1] <- p[1]+2*a/t
  findistr::dmerton(x, t, param)-exp(-eta(a, t, param))*findistr::dmerton(x, t, p)
}

#' The survival probability of a boundary by a given time under Merton's jump diffusion
#'
#' @param t the current time
#' @param a the boundary to hit above
#' @param param the dynamics of the jump-diffusion, a vector of 5, the drift, volatility
#' mean-jump rate, mean jump-size, sd of jump size
#'
#' @description {The survival probability of hitting \code{a} by time \code{t}. By survival,
#' we mean the process is not been stopped yet, i.e. has not hit above the boundary.}
#' @return numeric or vector
survival_merton <- function(t, a, param)
{
  y <- matrix(0, nrow = length(t))
  for(i in 1:length(t))
  {
    y[i] <- stats::integrate(trans_den_merton, lower = -Inf, upper = a, t = t[i], a = a, param = param)$value
  }
  return(y)
}

#' The hitting time CDF of a Merton jump-diffusion
#'
#' @param t the current time
#' @param a the boundary to hit above
#' @param param the dynamics of the jump-diffusion, a vector of 5, the drift, volatility
#' mean-jump rate, mean jump-size, sd of jump size
#'
#' @description {The CDF of the hitting time of \code{a} under Merton's jump-diffusion, i.e.
#' a geometric Levy process with the log being an Ito diffusion plus a compound Poisson
#' process of Gaussians.}
#' @return numeric or vector
#' @export phit_merton
phit_merton <- function(t, a, param)
{
  return(1-survival_merton(t, a, param))
}

#' The hitting time PDF of a Merton jump-diffusion
#'
#' @param t the current time
#' @param a the boundary to hit above
#' @param param the dynamics of the jump-diffusion, a vector of 5, the drift, volatility
#' mean-jump rate, mean jump-size, sd of jump size
#' @param n the number used in the discretization for differentiating the CDF with
#' forward differences
#'
#' @description {The PDF of the hitting time of \code{a} under Merton's jump-diffusion, i.e.
#' a geometric Levy process with the log being an Ito diffusion plus a compound Poisson
#' process of Gaussians. Most likely requires smoothing.}
#' @return numeric or vector
#' @export dhit_merton
dhit_merton <- function(t, a, param, n = 200)
{
  tt <- seq(0.0001, t, length.out = n+1)
  ff <- survival_merton(tt, a, param)
  rawDiff <- -diff(ff)/diff(tt)[1]

  return(data.frame(t = tt[-1], f = rawDiff))
}



