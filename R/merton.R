#' The probability density function under Merton's jump-diffusion model for stock-dynamics
#'
#' @param x Input value for the CDF
#' @param t The time input, for daily use 1/252
#' @param drift The drift rate
#' @param volat the diffusion coefficient
#' @param lambda Mean rate of jumps per year
#' @param a The mean jump size
#' @param b the mean jump volatility
#' @description {Auxillary function for PDF for a Jump Diffusion process. Not necessarily under the martingale measure.}
#' @return Numerical
dmerton1 <- function(x, t, drift, volat, lambda, a, b)
{
  # 200 should be enough for negligible probabilities in the Poisson sum
  n <- 0:200
  p <- stats::dpois(n, lambda = lambda*t)
  phi <- stats::dnorm(x, drift*t+n*a, sqrt(t*volat^2+n*b^2))
  sum(p*phi)
}

#' The probability density function under Merton's jump-diffusion model for stock-dynamics
#'
#' @param x Input value for the CDF
#' @param t The time input, for daily use 1/252
#' @param drift The drift rate
#' @param volat the diffusion coefficient
#' @param lambda Mean rate of jumps per year
#' @param a The mean jump size
#' @param b the mean jump volatility
#' @description {PDF for a Jump Diffusion process. Not necessarily under the martingale measure.}
#' @return Numerical
#' @export dmerton
dmerton <- function(x, t, drift, volat, lambda, a, b)
{
  if(length(x)==1)
  {
    dmerton1(x, t, drift, volat, lambda, a, b)
  } else
  {
    mapply(function(X){
      dmerton1(X, t, drift, volat, lambda, a, b)
    }, X = x)
  }
}


#' Estimate parameters for Merton jump-diffusion model via MLE
#'
#' @param log_returns log returns of data
#' @param thresh_hold threshhold for jump size
#' @param tn the length in years of period the data was observed in
#' @param time_step time step of data, for real stock data use 1/252, for simulations use \code{tn/n}.
#' @param compensate whether to treat the mean as compensated by the jump-drift or not
#' @param direction "down", "up" or "both" for direction of jumps to estimate
#'
#' @description {Estimate parameters for Merton's jump-diffusion model given a set of variates, using MLE and R's basic \code{optim} function.}
#' @return list
#' @export fitMerton
fitMerton <- function(log_returns, thresh_hold = 0.02, tn = 1, time_step = 1/252, compensate = FALSE, direction = "both")
{

  if(direction == "both")
  {
    jumps <- abs(log_returns) > thresh_hold
  } else if(direction == "down")
  {
    jumps <- log_returns < -thresh_hold
  } else if(direction == "up")
  {
    jumps <- log_returns > thresh_hold
  }

  lambda1 <- sum(jumps)/tn
  if(!compensate)
  {
    volat1 <- stats::sd(log_returns[!jumps])/(sqrt(time_step))
    mu1 <- (2*mean(log_returns[!jumps]+(time_step)*volat1^2))/(2*time_step)
    volatj <- sqrt(stats::var(log_returns[jumps])-(time_step)*volat1^2)
    muj <- mean(log_returns[jumps])-(mu1-volat1^2/2)*(time_step)
  } else{
    volat1 <- stats::sd(log_returns[!jumps])/(sqrt(time_step))
    volatj <- stats::sd(log_returns[jumps])
    muj <- mean(log_returns[jumps])
    eta1 <- exp(muj+0.5*volatj^2)-1
    mu1 <- mean(log_returns[!jumps])/time_step+0.5*volat1^2+lambda1*eta1
  }

  initial_v <- c(mu1, volat1, lambda1, muj, volatj)
  print("Initial estimate")
  print(initial_v)

  Sys.sleep(2)

  log_lik <- function(v)
  {
    eta <- exp(v[4]+0.5*v[5]^2)-1
    -sum(log(dmerton(log_returns, t = time_step, drift = v[1]-v[3]*eta-0.5*v[2]^2, volat = v[2], lambda = v[3], a = v[4], b = v[5])))
  }
  mle <- stats::optim(par = initial_v, fn = log_lik, method = "L-BFGS-B", control = list(trace = 5), lower = c(-1, 0.001, 0, -1, 0.001))
  return(mle)
}
