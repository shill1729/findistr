% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kou-family.R
\name{rdkou}
\alias{rdkou}
\title{Simulate displaced double-exponential jumps}
\usage{
rdkou(n, p, alpha, beta, ku, kd)
}
\arguments{
\item{n}{number of variates to simulate}

\item{p}{mixing probability}

\item{alpha}{mean size of positive jumps}

\item{beta}{mean size of negative jumps}

\item{ku}{the displacement of upward jumps away from the origin}

\item{kd}{the displacement of downward jumps away from the origin}
}
\value{
vector
}
\description{
{Simulates mixture of two exponentials one negative.}
}
\references{
The double-exponential mixture distribution or the Kou distribution was introduced (in mathematical finance) in the paper \href{http://www.columbia.edu/~sk75/MagSci02.pdf}{http://www.columbia.edu/~sk75/MagSci02.pdf} by S.G. Kou.
}
