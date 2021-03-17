
# findistr

<!-- badges: start -->
<!-- badges: end -->

The package findistr provides a collection of distributions and fitting routines used in financial mathematics. Probability density functions, cumulative distribution functions, and random-generators are available for most of the distributions implemented here. Also provided are various fitting rountines.

## Installation

You can install the github version via devtools:

``` r
devtools::install_github("shill1729/findistr")
```


## Discrete Time
The following distributions are included for discrete-time modeling:
Uniform, Gaussian, Gaussian-mixture, and the stable distribution.
## Continuous time
For continuous-time stochastic processes we can fit the log-normal distribution as parameterized by Black-Scholes, the Geometric Poisson model, jump-diffusions with Gaussian, Kou, or uniform jumps, a continuous-time Gaussian mixture, and multi-variate Black-Scholes.
### Hitting times
A test for hitting time densities to go above a log-price level under Merton's jump-diffusion and comparison to GBM
hitting times:
```r
mu <- 0.7027
volat <- 0.5483
lambda <- 21.74
alpha <- 0.0
beta <- 0.1253
a <- log(1+0.10)
maturity <- 30/252
n <- 200
param <- c(mu, volat, lambda, alpha, beta)
g <- dhit_merton(maturity, a, param, n)
# Compute smooth and GBM denisty
# g <- head(g, 100)
smoothed <- smooth(g$f)
gbmHit <- findistr::dhit_gbm(g$t, c(mu, volat), exp(a), 1)
# Plotting PDF and CDF
par(mfrow = c(1, 2))

plot(g, type = "l", ylim = c(0, max(g$f, gbmHit, smoothed)), main = "Hitting time density")
lines(g$t, smoothed, col = "red", lty = "dashed")
lines(g$t, gbmHit, col = "blue", lty = "dashed")
legend("topright", legend = c("Raw Merton", "Smoothed Merton", "GBM"),
       col = c("black", "red", "blue"), lty = 1, cex = 0.5)
plot(g$t, phit_merton(g$t, a, param), type = "l", main = "CDF")
lines(g$t, findistr::phit_gbm(g$t, c(mu, volat), exp(a), 1), col = "blue", lty = "dashed")
legend("topleft", legend = c("Merton", "GBM"),
       col = c("black", "blue"), lty = 1, cex = 0.5)
```
As one can see from running this, the raw density is periodically unstable but the CDF is fine.
## Fitting routines
Fitting routines are only available for the discrete time models, Black-Scholes, Gaussian-mixture, Merton's jump diffusion (ad-hoc), and the multivariate Black-Scholes.
