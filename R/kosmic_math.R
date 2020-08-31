#' Box-Cox Transform
#'
#' @param y A numeric vector of untransformed values.
#' @param ty A numeric vector of transformed values.
#' @param lambda The lambda parameter for the Box-Cox transform.
#' @param tolerance A numeric. If the absolute value of lambda is less than the
#'   tolerance, it is treated as if it were zero.
#'
#' @return A vector the same length as `y` or `ty`.
#'
#' @examples
#' y <- seq(1:10)
#' ty <- boxcox(y, 2.5)
#' y2 <- boxcox_inverse(ty, 2.5)
boxcox <- function(y, lambda, tolerance = 1e-6) {
  if (abs(lambda) < tolerance) {
    log(y)
  } else {
    ((y ^ lambda) - 1) / lambda
  }
}

#' @rdname boxcox
boxcox_inverse <- function(ty, lambda, tolerance = 1e-6) {
  if (abs(lambda) < tolerance) {
    exp(ty)
  } else {
    (lambda * ty + 1) ^ (1 / lambda)
  }
}


#' @rdname boxcox-distribution
#' @export
dboxcox <- function(x, mean, sd, lambda) {
  dnorm(boxcox(x, lambda), mean, sd)
}

#' Box-Cox Transformed Normal Distribution
#' 
#' @rdname boxcox-distribution
#' @export
pboxcox <- function(q, mean, sd, lambda, lower.tail = TRUE) {
  pnorm(boxcox(q, lambda), mean, sd, lower.tail)
}

#' @rdname boxcox-distribution
#' @export
qboxcox <- function(p, mean, sd, lambda, lower.tail = TRUE) {
  boxcox_inverse(qnorm(p, mean, sd, lower.tail), lambda)
}

#' @rdname boxcox-distribution
#' @export
rboxcox <- function(n, mean, sd, lambda) {
  boxcox_inverse(rnorm(n, mean, sd),)
}