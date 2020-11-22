#' Simulate double-exponential jumps
#'
#' @param n number of variates to simulate
#' @param p mixing probability
#' @param alpha mean size of positive jumps
#' @param beta mean size of negative jumps
#'
#' @description {Simulates mixture of two exponentials one negative.}
#' @return vector
#' @export rkou
#' @references The double-exponential mixture distribution or the Kou distribution was introduced (in mathematical finance) in the paper \href{http://www.columbia.edu/~sk75/MagSci02.pdf}{http://www.columbia.edu/~sk75/MagSci02.pdf} by S.G. Kou.
rkou <- function(n, p, alpha, beta)
{
  if(p <= 0 || p > 1)
  {
    stop("argument 'p' must be in the unit-interval")
  }
  if(alpha <= 0 || beta <= 0)
  {
    stop("means must be positive")
  }

  X <- stats::rexp(n, rate = 1/alpha)
  Y <- stats::rexp(n, rate = 1/beta)
  B <- stats::rbinom(n, size = 1, prob = p)
  Z <- B*X-(1-B)*Y
  return(Z)
}

#' PDF of double-exponential jumps
#'
#' @param x real-input of the PDF
#' @param p mixing probability
#' @param alpha mean size of positive jumps
#' @param beta mean size of negative jumps
#'
#' @description {PDF of mixture of two exponentials one negative.}
#' @return vector
#' @export dkou
#' @references The double-exponential mixture distribution or the Kou distribution was introduced (in mathematical finance) in the paper \href{http://www.columbia.edu/~sk75/MagSci02.pdf}{http://www.columbia.edu/~sk75/MagSci02.pdf} by S.G. Kou.
dkou <- function(x, p, alpha, beta)
{
  if(p <= 0 || p > 1)
  {
    stop("argument 'p' must be in the unit-interval")
  }
  if(alpha <= 0 || beta <= 0)
  {
    stop("means must be positive")
  }
  p*stats::dexp(x, rate = 1/alpha)*ifelse(x >= 0, 1, 0)+(1-p)*stats::dexp(-x, rate = 1/beta)*ifelse(x < 0, 1, 0)
}

#' CDF of double-exponential jumps
#'
#' @param x real-input of the PDF
#' @param p mixing probability
#' @param alpha mean size of positive jumps
#' @param beta mean size of negative jumps
#'
#' @description {CDF of mixture of two exponentials one negative.}
#' @return vector
#' @export pkou
#' @references The double-exponential mixture distribution or the Kou distribution was introduced (in mathematical finance) in the paper \href{http://www.columbia.edu/~sk75/MagSci02.pdf}{http://www.columbia.edu/~sk75/MagSci02.pdf} by S.G. Kou.
pkou <- function(x, p, alpha, beta)
{
  if(p <= 0 || p > 1)
  {
    stop("argument 'p' must be in the unit-interval")
  }
  if(alpha <= 0 || beta <= 0)
  {
    stop("means must be positive")
  }
  p*stats::pexp(x, rate = 1/alpha)*ifelse(x >= 0, 1, 0)+(1-p)*stats::pexp(-x, rate = 1/beta)*ifelse(x < 0, 1, 0)
}

#' Simulate displaced double-exponential jumps
#'
#' @param n number of variates to simulate
#' @param p mixing probability
#' @param alpha mean size of positive jumps
#' @param beta mean size of negative jumps
#' @param ku the displacement of upward jumps away from the origin
#' @param kd the displacement of downward jumps away from the origin
#'
#' @description {Simulates mixture of two exponentials one negative.}
#' @return vector
#' @export rdkou
#' @references The double-exponential mixture distribution or the Kou distribution was introduced (in mathematical finance) in the paper \href{http://www.columbia.edu/~sk75/MagSci02.pdf}{http://www.columbia.edu/~sk75/MagSci02.pdf} by S.G. Kou.
rdkou <- function(n, p, alpha, beta, ku, kd)
{
  if(p <= 0 || p > 1)
  {
    stop("argument 'p' must be in the unit-interval")
  }
  if(alpha <= 0 || beta <= 0)
  {
    stop("means must be positive")
  }

  X <- stats::rexp(n, rate = 1/alpha)+ku
  Y <- stats::rexp(n, rate = 1/beta)-kd
  B <- stats::rbinom(n, size = 1, prob = p)
  Z <- B*X-(1-B)*Y
  return(Z)
}

#' PDF of displaced double-exponential jumps
#'
#' @param x real-input of the PDF
#' @param p mixing probability
#' @param alpha mean size of positive jumps
#' @param beta mean size of negative jumps
#' @param ku the displacement of upward jumps away from the origin
#' @param kd the displacement of downward jumps away from the origin
#'
#'
#' @description {PDF of mixture of two exponentials one negative.}
#' @return vector
#' @export ddkou
#' @references The double-exponential mixture distribution or the Kou distribution was introduced (in mathematical finance) in the paper \href{http://www.columbia.edu/~sk75/MagSci02.pdf}{http://www.columbia.edu/~sk75/MagSci02.pdf} by S.G. Kou.
ddkou <- function(x, p, alpha, beta, ku, kd)
{
  if(p <= 0 || p > 1)
  {
    stop("argument 'p' must be in the unit-interval")
  }
  if(alpha <= 0 || beta <= 0)
  {
    stop("means must be positive")
  }
  if(ku < 0)
  {
    stop("Upward displacement 'ku' must be non-negative")
  }
  if(kd > 0)
  {
    stop("Downward displacement 'kd' must be non-positive")
  }
  p*stats::dexp(x-ku, rate = 1/alpha)*ifelse(x >= ku, 1, 0)+(1-p)*stats::dexp(-(x-kd), rate = 1/beta)*ifelse(x < kd, 1, 0)
}

#' CDF of displaced double-exponential jumps
#'
#' @param x real-input of the PDF
#' @param p mixing probability
#' @param alpha mean size of positive jumps
#' @param beta mean size of negative jumps
#' @param ku the displacement of upward jumps away from the origin
#' @param kd the displacement of downward jumps away from the origin
#'
#' @description {CDF of mixture of two exponentials one negative.}
#' @return vector
#' @export pdkou
#' @references The double-exponential mixture distribution or the Kou distribution was introduced (in mathematical finance) in the paper \href{http://www.columbia.edu/~sk75/MagSci02.pdf}{http://www.columbia.edu/~sk75/MagSci02.pdf} by S.G. Kou.
pdkou <- function(x, p, alpha, beta, ku, kd)
{
  if(p <= 0 || p > 1)
  {
    stop("argument 'p' must be in the unit-interval")
  }
  if(alpha <= 0 || beta <= 0)
  {
    stop("means must be positive")
  }
  if(ku < 0)
  {
    stop("Upward displacement 'ku' must be non-negative")
  }
  if(kd > 0)
  {
    stop("Downward displacement 'kd' must be non-positive")
  }
  p*stats::pexp(x-ku, rate = 1/alpha)*ifelse(x >= ku, 1, 0)+(1-p)*stats::pexp(-(x-kd), rate = 1/beta)*ifelse(x < kd, 1, 0)
}
