#' Create a new kosmic result
#'
#' @description
#'
#' This creates an S3 object to hold the results of the kosmic library being run
#' on a set of data. It should be called by user-facing methods such as
#' `kosmic()`.
#'
#' @param n A positive number. The number of results used to estimate the
#'   distribution.
#' @param l A number. The lambda parameter of the Box-Cox transformation for the
#'   estimated distribution. It which describes the skewness.
#' @param mean A number. The mean parameter shows the central point of the
#'   estimated distribution before it undergoes Box-Cox transformation.
#' @param sd A number. The spread of the estimated distribution before it
#'   undergoes Box-Cox transformation.
#' @param T1 A number.
#' @param T2 A number.
#' @param decimals A number. The number of digits of precision the kosmic
#'   algorithm was configured for.
#' @param T1min A number.
#' @param T1max A number.
#' @param T2min A number.
#' @param T2max A number.
#'
#' @return
#'
#' A list with the class `kosmic`, containing the runtime parameters and results
#' of the kosmic algorithm.
new_kosmic <- function(n,
                       l,
                       mean,
                       sd,
                       T1,
                       T2,
                       decimals,
                       T1min,
                       T1max,
                       T2min,
                       T2max) {
  if(!is_bare_numeric(n, n=1) | n <= 0) {
    abort("`n` must be a single positive integer.")
  }
  if(!is_bare_numeric(decimals, n=1) | decimals <= 0) {
    abort("`decimals` must be a single positive integer.")
  }
  for (arg in exprs(mean, sd, T1, T2, decimals,
                    T1min, T1max, T2min, T2max)) {
    if(!is_bare_numeric(eval(arg), n=1)) {
      abort(glue("`{arg}` must be a single numeric value."))
    }
  }
  
  elems <- as.list(match.call())
  
  structure(elems, class = c("kosmic"))
}
