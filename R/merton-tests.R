#' Verify the empirical density of Merton simulated variates
#'
#' @param n number of variates to simulate
#' @param t the terminal time-horizon
#' @param param the parameter set defining the jump diffusion dynamics
#'
#' @description {Compares the empirical density of simulated variates against
#' the exact density per the infinite mixture formula, for Merton log dynamics.}
#'
#' @return vector of errors
#' @export testMertonSimDensity
testMertonSimDensity <- function(n, t, param)
{
  x <- rmerton(n, t, param)
  epdf <- stats::density(x)
  exact <- dmerton(epdf$x, t = t, param)
  ae <- abs(exact-epdf$y)
  se <- ae^2
  mae <- mean(ae)
  mse <- mean(se)
  mre <- mean(ae[ae>0]/exact[exact>0])
  errors <- data.frame(mae, mre, mse)
  # Plot estimated density from simulation vs model density
  plot(epdf$x, epdf$y, type = "l", ylim = c(0, max(epdf$y, exact)))
  graphics::lines(epdf$x, exact, col = "blue", lty = "dashed")
  return(errors)
}
