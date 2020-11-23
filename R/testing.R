# drift <- 0.1
# volat <- 0.33
# lambda <- 5
# a <- 0
# b <- 0.02
# x <- 0.03
# t <- 0.03
# n <- 150
#
# param <- c(drift, volat, lambda, a, b)
# # Normal
# findistr:::pdfNormal(x, drift, volat)
# dnorm(x, drift, volat)
# # poisson
# findistr:::pdfPoisson(n, lambda*t)
# dpois(n, lambda*t)
# # Merton
# findistr:::pdfMerton(x, t, param)
# dmerton(x, t, drift, volat, lambda, a, b)
#
# # Microbenchmark
# y <- findistr:::pdfMerton(x, t, param)
# w <- dmerton(x, t, drift, volat, lambda, a, b)
# microbenchmark::microbenchmark(cpp = y <- findistr:::pdfMerton(x, t, param),
#                                r = w <- dmerton(x, t, drift, volat, lambda, a, b),
#                                times = 30)
