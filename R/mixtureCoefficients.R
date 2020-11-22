#' Time-homotopy between a constant and a given function of (s,t)
#'
#' @param s the current price
#' @param t the time
#' @param f function to homotopy towards
#' @param initial initial constant
#' @param epsilon time-cut off
#' @param ... additional parameters for the function \code{f(s, t, ...)}
#'
#' @return numeric
#' @export time_homotopy
time_homotopy <- function(s, t, f, initial, epsilon = 0.5/252,...)
{
  if(t >= 0 && t < epsilon)
  {
    return(initial)
  } else if(t >= epsilon && t < 2*epsilon)
  {
    homotopy <- initial*(1-(t-epsilon)/epsilon)+f(s, t,...)*(t-epsilon)/epsilon
  } else if(t >= 2*epsilon)
  {
    return(f(s, t, ...))
  }
}


#' local drift mixture coefficient function
#'
#' @param s the price
#' @param t the time, in trading years
#' @param probs the probability of each mixing component
#' @param mus the drift of each mixing component
#' @param sigmas the volatility of each mixing component
#' @param spot the spot value
#' @description {A local volatility mixture of constant volatilities, with linear behavior near time zero.}
#' @details {See .pdf for derivation}
#' @return numeric
drift_lvm1 <- function(s, t, probs, mus, sigmas, spot)
{
  # Component PDF
  p_i <- stats::dnorm(log(s/spot), mean = (mus-0.5*sigmas^2)*t, sd = sqrt(t)*sigmas)/s
  # Mixture PDFs
  p <- sum(probs*p_i)
  # Mixture weight for coef
  mix <- probs*p_i/p
  # drift coefficient
  mu <- sum(mix*mus)
  return(mu)
}

#' local drift mixture coefficient function
#'
#' @param s the price
#' @param t the time, in trading years
#' @param probs the probability of each mixing component
#' @param mus the drift of each mixing component
#' @param sigmas the volatility of each mixing component
#' @param spot the spot value
#' @param spot_cut the truncation for behavior near \code{s=0}
#'
#' @description {A local volatility mixture of constant volatilities, with linear behavior near time zero.}
#' @details {See .pdf for derivation}
#' @return numeric
#' @export drift_lvm
drift_lvm <- function(s, t, probs, mus, sigmas, spot, spot_cut = 0.5)
{
  if(s < spot_cut)
  {
    return(min(mus))
  } else{
    return(time_homotopy(s, t, f = drift_lvm1, initial = min(mus), epsilon = 0.5/252, probs = probs, mus = mus, sigmas = sigmas, spot = spot))
  }

}




#' local volatility mixture coefficient function
#'
#' @param s the price
#' @param t the time, in trading years
#' @param probs the probability of each mixing component
#' @param mus the drift or risk-free rate (depending on measure) can also be a mixture
#' @param sigmas the volatility of each mixing component
#' @param spot the spot value
#' @description {A local volatility mixture of constant volatilities, with linear behavior near time zero.}
#' @details {See .pdf for derivation}
#' @return numeric
volat_lvm1 <- function(s, t, probs, mus, sigmas, spot)
{
  # Component PDF
  p_i <- stats::dnorm(log(s/spot), mean = (mus-0.5*sigmas^2)*t, sd = sqrt(t)*sigmas)/s
  # Mixture PDF
  p <- sum(probs*p_i)
  # Mixture weight for coef
  mix <- probs*p_i/p
  # Volatility coefficient
  volat <- sqrt(sum(mix*(sigmas^2)))
  return(volat)
}

#' local volatility mixture coefficient function
#'
#' @param s the price
#' @param t the time, in trading years
#' @param probs the probability of each mixing component
#' @param mus the drift or risk-free rate (depending on measure) can also be a mixture
#' @param sigmas the volatility of each mixing component
#' @param spot the spot value
#' @param spot_cut the truncation of the price near \code{s=0}
#'
#' @description {A local volatility mixture of constant volatilities, with linear behavior near time zero.}
#' @details {See .pdf for derivation}
#' @return numeric
#' @export volat_lvm
volat_lvm <- function(s, t, probs, mus, sigmas, spot, spot_cut = 0.5)
{
  if(s < spot_cut)
  {
    return(max(sigmas))
  } else
  {
    return(time_homotopy(s, t, f = volat_lvm1, initial = max(sigmas), epsilon = 0.5/252, probs = probs, mus = mus, sigmas = sigmas, spot = spot))
  }

}
