# # Instability issues near  (0,0)?
# # The MLE is not bad for long time periods,
# # but it is hard to get good initial guesses for (a,b)
# # while lambda is not too difficult to produce an ad-hoc
# # estimate.
# library(findistr)
# n <- 10000
# tt <- 1
# spot <- 100
# a <- 1
# b <- 0.2
# lambda <- 30
# # Simulate terminal variates: S_t = S_0 exp(X_t)
# s <- rgcpp(n, tt, spot, a, b, lambda)
# # X_t = log(S_t/S_0)=aN_t-bt has range for t>0
# # -bt, a-bt, 2a-bt, 3a-bt, ...
# ab <- range(log(s/spot))
#
# x <- seq(ab[1], ab[2], by = a)
# # The empirical and model PMF for X_t
# epmf <- unlist(lapply(x, function(y) mean(log(s/spot)==y)))
# mpmf <- dgcpp(exp(x), tt, 1, a, b, lambda)
#
# # Log-likelihood function:
# logLikeGCPP <- function(v)
# {
#   aa <- v[1]
#   bb <- v[2]
#   ll <- v[3]
#   f <- dgcpp(s, tt, spot, aa, bb, ll)
#   l <- sum(log(f[f>0]))
#   return(l)
# }
# # Base R numerical optimization of LL for MLE:
# # Initial guesses are tough!?
# w <- optim(par = c(1, 1, mean(log(s/spot)>x[1])/tt), fn = logLikeGCPP,
#       method = "L-BFGS-B", lower = c(0.01, 0.01, 0.01),
#       control = list(fnscale = -1))
# # Compute PMF with estimated parameters
# estpmf <- dgcpp(exp(x), tt, 1, w$par[1], w$par[2], w$par[3])
# # Plotting
# par(mfrow = c(1, 1))
# plot(x, epmf, type = "h", ylim = c(0, max(epmf, mpmf, estpmf)))
# lines(x, mpmf, type = "h", col = "blue", lty = "dashed")
# lines(x, estpmf, type = "h", col = "red", lty = "dashed")
# # Print fit:
# parameters <- rbind(c(a, b, lambda), w$par)
# rownames(parameters) <- c("exact", "fit")
# print(parameters)
#
#
#
#
