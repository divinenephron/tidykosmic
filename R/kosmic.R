#' @rdname kosmic
#' @export
kosmic <- function(data, ...) {
  UseMethod("kosmic")
}

#' @export
kosmic.default <- function(data, ...) {
  abort(glue("No `kosmic` method is defined for the class `{class(data)[1]}'."))
}

#' Estimate a Distribution of Physiological Results Using Kosmic
#'
#' @description Estimates the distribution of physiological results from a mixed
#'   distribution of physiological and abnormal results, such as those found in
#'   laboratory databases.
#'
#' @details The kosmic algorithm works by minimizes the difference between an
#'   estimated parametric distribution and a truncated part of the observed
#'   distribution, specifically, the Kolmogorov-Smirnov-distance between a
#'   hypothetical Gaussian distribution and the observed distribution of test
#'   results after Box-Cox-transformation.
#'
#' @param data A numeric vector. A mixed distribution of physiological and
#'   abnormal results.
#' @param decimals A positive integer. The number of digits of precision to
#'   calculate the results to. Increasing this makes the algorithm take longer.
#' @param t1min A quantile. Start of the search range for T1.
#' @param t1max A quantile. End of the search range for T1.
#' @param t2min A quantile. Start of the search range for T2.
#' @param t2max A quantile. End of the search range for T1.
#' @param sd_guess A quantile. The quantile used for the initial guess of the
#'   standard deviation.
#' @param tol The absolute convergence tolerance for the optimizer. The
#'   algorithm stops if it is unable to reduce the cost by more than this
#'   amount.
#' @param threads A positive integer >= 1. The number of threads to use. Default
#'   value is 1. Setting this to a higher number may speed up computation, but
#'   not always. If a value higher than 1 is used the results may not be
#'   reproduceable from run-to-run.
#' @param na.rm Logical. If true, any NA and NaNs are removed from `data` before
#'   calling kosmic.
#'
#' @return `kosmic` returns an object of class "kosmic". This contains the
#'   following components:
#'
#'   \code{n} The number of data points used to estimate the distribution.
#'
#'   \code{estimates} A named vector of estimates for the ditribution:
#'   \code{lambda}, \code{mean}, \code{sd}.
#'
#'   \code{truncation} A named vector of estimate for truncating the original
#'   data: \code{t1} and \code{t2}.
#'
#'   \code{settings} A named vector of the settings passed to \code{kosmic}.
#'
#' @examples
#' set.seed(1)
#' k <- kosmic(rnorm(10000, 16, 1), 1)
#' quantile(k)
#'
#' @rdname kosmic
#' @export
kosmic.numeric <- function(data,
                           decimals,
                           t1min = 0.05,
                           t1max = 0.30,
                           t2min = 0.70,
                           t2max = 0.95,
                           sd_guess = 0.80,
                           abstol = 1e-7,
                           threads=1L,
                           na.rm = FALSE,
                           ...) {
  if (na.rm) 
    data <- data[!is.na(data)]
  else if (anyNA(data)) 
    stop("missing values and NaN's not allowed if 'na.rm' is FALSE")
  
  kosmic_bridge(data, decimals, t1min, t1max, t2min, t2max,
                sd_guess, abstol, threads)
}

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
#' @param lambda A number. The lambda parameter of the Box-Cox transformation
#'   for the estimated distribution. It which describes the skewness.
#' @param mean A number. The mean parameter shows the central point of the
#'   estimated distribution before it undergoes Box-Cox transformation.
#' @param sd A number. The spread of the estimated distribution before it
#'   undergoes Box-Cox transformation.
#' @param t1 A number.
#' @param t2 A number.
#' @param decimals A number. The number of digits of precision the kosmic
#'   algorithm was configured for.
#' @param t1min A number.
#' @param t1max A number.
#' @param t2min A number.
#' @param t2max A number.
#' @param sd_guess A number.
#' @param abstol A number.
#'
#' @return
#'
#' A list with the class `kosmic`, containing a named vector of parameters for
#' the estimated ditribution, a named vector of truncation limits, a named
#' vector of kosmic settings, the number of original data points, and a
#' frequency table of the original data.
new_kosmic <- function(data,
                       n,
                       lambda,
                       mean,
                       sd,
                       t1,
                       t2,
                       decimals,
                       t1min,
                       t1max,
                       t2min,
                       t2max,
                       sd_guess,
                       abstol) {
  if(!is.data.frame(data)) {
    abort("`data` must be a data frame.")
  }
  if(!is_bare_numeric(n, n=1) | n <= 0) {
    abort("`n` must be a single positive integer.")
  }
  if(!is_bare_numeric(decimals, n=1) | decimals <= 0) {
    abort("`decimals` must be a single positive integer.")
  }
  for (arg in exprs(mean, sd, t1, t2, decimals,
                    t1min, t1max, t2min, t2max,
                    sd_guess, abstol)) {
    if(!is_bare_numeric(eval(arg), n=1)) {
      abort(glue("`{arg}` must be a single numeric value."))
    }
  }
  
  estimates <- c(lambda = lambda,
                mean = mean,
                sd = sd)
  truncation <- c(t1 = t1,
                  t2 = t2)
  settings <- c(decimals = decimals,
                t1min = t1min,
                t1max = t1max,
                t2min = t2min,
                t2max = t2max,
                sd_guess = sd_guess,
                abstol = abstol)
  elems <- list(data = data,
                n = n,
                estimates = estimates,
                truncation = truncation,
                settings = settings)
  structure(elems, class = c("kosmic"))
}

#' Run Kosmic and Create an Object to Hold the Results
#'
#' @return
#' A `kosmic` object.
kosmic_bridge <- function(data,
                          decimals,
                          t1min,
                          t1max,
                          t2min,
                          t2max,
                          sd_guess,
                          abstol,
                          threads) {
  if(!is.numeric(data)) {
    abort("`data` must be a numeric vector.")
  }
  if(!is_bare_numeric(decimals, n = 1) | decimals < 0) {
    abort("`decimals` must be a single number >= 0.")
  }
  if(!is_bare_numeric(abstol, n = 1) | abstol <= 0) {
    abort("`abstol` must be a single number > 0.")
  }
  if(!is_bare_numeric(threads, n = 1) | threads < 1) {
    abort("`threads` must be a single integer >= 1.")
  }
  # Check quantile are quantiles
  for (arg in exprs(t1min, t1max, t2min, t2max, sd_guess)) {
    if(!is_bare_numeric(eval(arg), n=1)) {
      abort(glue("`{arg}` must be a single number."))
    }
    if(eval(arg) > 1 | eval(arg) < 0) {
      abort(glue("`{arg}` must be between 0 and 1 (inclusive)."))
    }
  }

  # Run the Kosmic algorithm, written in C++
  bootstrap <- 0
  bootstrap_seed <- get_kosmic_seed()
  impl_result <- kosmic_impl(data,
                             trunc(decimals),
                             trunc(bootstrap),
                             bootstrap_seed,
                             trunc(threads),
                             t1min, t1max,
                             t2min, t2max,
                             sd_guess, abstol)
  res <- impl_result$result
  
  # Include a frequency table of the original data to allow plotting
  # later on
  freqs <- data.frame(result = round(data, trunc(decimals))) %>%
    mutate(result = round(result, trunc(decimals))) %>%
    count(result)
  
  # Create a kosmic object to hold the results
  new_kosmic(freqs,
             n = length(data),
             lambda = res[1],
             mean = res[2],
             sd = res[3],
             t1 = res[5],
             t2 = res[6],
             decimals = decimals,
             t1min = t1min,
             t1max = t1max,
             t2min = t2min,
             t2max = t2max,
             sd_guess = sd_guess,
             abstol = abstol)
}

#' @export
print.kosmic <- function(x, ...) {
  if (!inherits(x, "kosmic")) {
    abort("Use only with `kosmic` objects")
  }
  
  cat("Distribution of physiological results estimated using kosmic \n\n")
  cat(glue("Number of data points: {x$n}"), "\n\n")
  cat("Distribution estimates:\n")
  print(data.frame("estimate" = x$estimates), ...)
  cat("\n")
  cat("Truncation estimates:\n")
  print(data.frame("estimate" = x$truncation), ...)
  cat("\n")
  cat("Settings:\n")
  print(data.frame("setting" = x$settings), ...)
  cat("\n")
  
  invisible(x)
}


#' Calculate Quantile for an Estimated Distribution of Physiological Results
#'
#' @param x A kosmic result.
#' @param probs Numeric vector. Probabilities with a value between 0 and 1
#'   (exclusive).
#' @param names Logical. If true, the result has a names attribute. Set to false
#'   to speed up the calculation of many probabilities.
#' @param ...
#'
#' @return
#' A vector of 
#' 
#' @export
quantile.kosmic <- function(x,
                            probs = c(0.025, 0.500, 0.975),
                            names = TRUE,
                            ...) {
  if (!inherits(x, "kosmic")) {
    abort("Use only with `kosmic` objects")
  }
  if(!is.numeric(probs) | any(probs >= 1) | any(probs <= 0)) {
    abort(glue("`probs` must all be numbers between 0 and 1 (exclusive)."))
  }
  
  res <- qboxcox(probs, x$estimates["mean"], x$estimates["sd"], x$estimates["lambda"])
  
  if (names && length(probs) > 0L) {
    names(res) <- format_perc(probs)
  }
  
  res
}

# Change probabilities like 0.1 into percentage strings like "10%"
format_perc <- function (x,
                         digits = max(2L,getOption("digits")),
                         probability = TRUE, 
                         ...) {
  if (length(x)) {
    x <- 100 * x
    paste0(format(x, trim = TRUE, digits = digits, ...), "%")
  }
  else character(0)
}


