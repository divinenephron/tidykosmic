% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kosmic_math.R
\name{boxcox}
\alias{boxcox}
\alias{boxcox_inverse}
\title{Box-Cox Transform}
\usage{
boxcox(y, lambda, tolerance = 1e-06)

boxcox_inverse(ty, lambda, tolerance = 1e-06)
}
\arguments{
\item{y}{A numeric vector of untransformed values.}

\item{lambda}{The lambda parameter for the Box-Cox transform.}

\item{tolerance}{A numeric. If the absolute value of lambda is less than the
tolerance, it is treated as if it were zero.}

\item{ty}{A numeric vector of transformed values.}
}
\value{
A vector the same length as \code{y} or \code{ty}.
}
\description{
Box-Cox Transform
}
\keyword{internal}
