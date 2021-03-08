
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
## Fitting routines
Fitting routines are only available for the discrete time models, Black-Scholes, Gaussian-mixture, Merton's jump diffusion (ad-hoc), and the multivariate Black-Scholes.
