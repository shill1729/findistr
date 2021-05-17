#' Conditional mean and volatility functions of state and observation
#'
#' @param s the volatility level
#' @param param the vector of parameters, see details
#' @param h the time step
#'
#' @description {The conditional mean and volatility functions for the Heston state and
#' observation.}
#' @details {The argument \code{param} must be a vector whose entries represent in order,
#' \eqn{(\kappa, \theta, \xi, \mu)}, the mean-reversion speed, mean-reversion level, the vol-of-vol,
#' and mean-drift of the stock price.}
#'
#' @return numeric
state_mean <- function(s, param, h = 1/252)
{
  # Volatility parameters: reversion speed, reversion level, vol-of-vol
  kappa <- param[1]
  theta <- param[2]
  return(kappa*theta*h+(1-kappa*h)*s)
}

#' Conditional mean and volatility functions of state and observation
#'
#' @param s the volatility level
#' @param param the vector of parameters, see details
#' @param h the time step
#'
#' @description {The conditional mean and volatility functions for the Heston state and
#' observation.}
#' @details {The argument \code{param} must be a vector whose entries represent in order,
#' \eqn{(\kappa, \theta, \xi, \mu)}, the mean-reversion speed, mean-reversion level, the vol-of-vol,
#' and mean-drift of the stock price.}
#'
#' @return numeric
state_vol <- function(s, param, h = 1/252)
{
  # Vol of vol
  xi <- param[3]
  return(xi*sqrt(h)*sqrt(s))
}


#' Conditional mean and volatility functions of state and observation
#'
#' @param s the volatility level
#' @param param the vector of parameters, see details
#' @param h the time step
#'
#' @description {The conditional mean and volatility functions for the Heston state and
#' observation.}
#' @details {The argument \code{param} must be a vector whose entries represent in order,
#' \eqn{(\kappa, \theta, \xi, \mu)}, the mean-reversion speed, mean-reversion level, the vol-of-vol,
#' and mean-drift of the stock price.}
#'
#' @return numeric
obs_mean <- function(s, param, h = 1/252)
{
  # mean drift of stock price
  mu <- param[4]
  return((mu-0.5*s)*h)
}


#' Conditional mean and volatility functions of state and observation
#'
#' @param s the volatility level
#' @param param the vector of parameters, see details
#' @param h the time step
#'
#' @description {The conditional mean and volatility functions for the Heston state and
#' observation.}
#' @details {The argument \code{param} must be a vector whose entries represent in order,
#' \eqn{(\kappa, \theta, \xi, \mu)}, the mean-reversion speed, mean-reversion level, the vol-of-vol,
#' and mean-drift of the stock price.}
#'
#' @return numeric
obs_vol <- function(s, param, h = 1/252)
{
  # Volatility of stock price
  return(sqrt(h)*sqrt(s))
}


#' Particle filter for filtering volatility in the Heston model
#'
#' @param y time-series of log-price increments
#' @param param vector of parameters defining Heston dynamics, see details
#' @param N number of particles to use
#' @param nthresh threshold of particles to use
#' @param h the time-step to use
#'
#' @description {Sequential importance (re)-sampling implementation of the particle method
#' to estimate the hidden volatility of a Heston model.}
#' @details {The argument \code{param} must be a vector whose entries represent in order,
#' \eqn{(\kappa, \theta, \xi, \mu)}, the mean-reversion speed, mean-reversion level, the vol-of-vol,
#' and mean-drift of the stock price. The algorithm's details are on Wikipedia.}
#'
#' @return vector/numeric
#' @export hestonParticleFilter
hestonParticleFilter <- function(y, param, N, nthresh, h = 1/252)
{
  # Volatility parameters: reversion speed, reversion level, vol-of-vol
  kappa <- param[1]
  theta <- param[2]
  # Feller condition for positivity
  xi <- param[3]
  # mean drift of stock price
  mu <- param[4]

  v0 <- stats::var(y)*252
  n <- length(y)
  v <- matrix(data = 0, nrow = n, ncol = N)
  w <- matrix(data = 0, nrow = n, ncol = N)
  v[1, ] <- stats::runif(N, min = v0/2, max = 2*v0)
  w[1, ] <- 1/N
  # Over every observation data-point
  for(k in 2:n)
  {
    # Generate N variates representing the state, conditional on
    # previous state: assumed Gaussian with specified m/v
    v[k, ] <- stats::rnorm(N, mean = state_mean(v[k-1, ], param, h), sd = state_vol(v[k-1, ], param, h))
    v[k, ] <- abs(v[k, ]) # ! positivity hack for the square-root

    # Weights are recursively defined with conditional distr of observation given
    # the previous states
    w[k, ] <- w[k-1, ]*stats::dnorm(y[k], mean = obs_mean(v[k-1, ], param, h), sd = obs_vol(v[k-1, ], param, h))
    # If all the weights are zero, then following two computations present issues and
    # will error out, typically with NA or NaN values.
    if(identical(sum(w[k, ]), 0)) # So let's check for that and enforce non-zero weights
    {
      w[k, ] <- rep(1/N, N)
    }
    # Normalize the weights
    w[k, ] <- w[k, ]/sum(w[k, ])
    # Resampling procedure
    neff <- 1/sum(w[k, ]^2)
    if(neff < nthresh)
    {
      v[k, ] <- sample(v[k, ], N, TRUE, prob = w[k, ])
      w[k, ] <- rep(1/N, N)
    }
  }
  # Now we compute our filtering estimate by
  # MC integration of conditional expectation
  vhat <- matrix(0, nrow = n)
  for(k in 1:(n))
  {
    vhat[k] <- (v[k, ])%*%w[k, ]
  }
  # Return volatility instead of variance
  vhat <- sqrt(vhat)
  return(vhat)
}

#' Quasi log-likelihood function for Heston model
#'
#' @param p vector of parameters, see details
#' @param y time-series of log-price increments
#' @param v time-series of volatility
#'
#' @description {The log-likelihood function of the Heston model, computed using
#' Bayesian recursion.}
#' @details {The argument \code{p} must be a vector whose entries represent in order,
#' \eqn{(\kappa, \theta, \xi, \mu)}, the mean-reversion speed, mean-reversion level, the vol-of-vol,
#' and mean-drift of the stock price.}
#' @return numeric
hestonLogLikelihood <- function(p, y, v)
{
  n <- length(v)
  x1 <- stats::dnorm(v[-1], state_mean(v[-n], p), state_vol(v[-n], p))
  x2 <- stats::dnorm(y[-1], obs_mean(v[-n], p), obs_vol(v[-n], p))
  t1 <- sum(log(x1[x1 > 0])) # avoiding undefined log-likelihood
  t2 <- sum(log(x2[x2 > 0 ]))
  z <- t1+t2
  return(-z)
}

#' Quasi MLE for Heston dynamics
#'
#' @param y time-series of log-price increments
#' @param iterations number of iterations to step through for filtering and
#' @param N number of particles to use in the particle-filter
#' @param nthresh the threshold of particles in the particle-filter
#' performing the MLE estimates
#'
#' @description {An ad-hoc step-wise QMLE routine for estimating parameters
#' for the Heston dynamics given an observational time-series of log-price increments.}
#' @details {At each step, hidden volatility is filtered into an estimate based on log-price observations,
#' and conditional on these observations, the likelihood is maximized over the parameter space. Then
#' on the next step, the volatility is re-filtered under the updated parameters until reaching the last observation.
#' The filtered volatility is tracked step-wise, so that each iteration produces a new filtered-state
#' appended to the previous filtered-states produced under the previous parameters.}
#' @return vector/numeric
#' @export hestonMLE
hestonMLE <- function(y, iterations = 1, N = 100, nthresh = 100)
{
  n <- length(y)
  k <- n-iterations
  # Initial guess for Heston parameters
  kappa <- 1
  theta <- stats::var(y[1:(k-1)])*252
  xi <- sqrt(2*kappa*theta)/2
  mu <- mean(y[1:(k-1)])*252+theta/2
  param <- c(kappa, theta, xi, mu)
  # Iteratively perform MLE
  for(i in k:(k+iterations))
  {
    # Check Feller condition
    if(2*param[1]*param[2] < param[3]^2)
    {
      warning("Feller condition violated; enforcing arbitrarily")
      param[3] <- sqrt(2*param[1]*param[2])/exp(1)
    }

    # Re-filter variance and only keep latest state
    if(i == k)
    {
      # Remember our PF returns volatility, so transform back to variance
      vhat <- hestonParticleFilter(y[1:i], param, N, nthresh)^2
    } else if(i > k)
    {
      # Remember our PF returns volatility, so transform back to variance
      w <- hestonParticleFilter(y[1:i], param, N, nthresh)^2
      # appending the latest state with previous time-series
      vhat <- c(vhat, w[i])
    }

    # Optimize the conditional likelihood
    # parameters are in order (kappa, theta, xi, mu)
    # mu can be any real number, the rest must be positive,
    # and xi ought to be bounded so as to ensure the Feller condition
    result <- stats::optim(par = param,
                    fn = hestonLogLikelihood,
                    y = y[1:i],
                    v = vhat,
                    method = "L-BFGS-B",
                    lower = c(0.01, 0.01, 0.01,-Inf),
                    upper = c(Inf, Inf, 2, Inf)
                    )
    # Update parameters
    param <- result$par
  }
  # TODO consider option for returning the filtered state history?
  return(param)
}
