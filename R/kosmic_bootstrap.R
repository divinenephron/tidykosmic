#' @rdname kosmic_bootstrap
#' @export
kosmic_bootstrap <- function(data, ...) {
  UseMethod("kosmic_bootstrap")
}

#' @rdname kosmic_bootstrap
#' @export
kosmic_bootstrap.default <- function(data, ...) {
  abort(glue("No `kosmic_bootstrap` method is defined for the class `{class(data)[1]}'."))
}

#' @rdname kosmic_bootstrap
#' @export
kosmic_bootstrap.numeric <- function(data,
                                   decimals,
                                   replicates,
                                   t1min = 0.05,
                                   t1max = 0.30,
                                   t2min = 0.70,
                                   t2max = 0.95,
                                   sd_guess = 0.80,
                                   abstol = 1e-7,
                                   na.rm = FALSE,
                                   ...) {
  if (missing(decimals)) {
    abort("Argument `decimals` is required")
  }
  if (missing(replicates)) {
    abort("Argument `replicates` is required")
  }
  if (na.rm) 
    data <- data[!is.na(data)]
  else if (anyNA(data)) 
    abort("missing values and NaN's not allowed if 'na.rm' is FALSE")
  
  kosmic_bootstrap_bridge(data, decimals, replicates, t1min, t1max, t2min, t2max,
                sd_guess, abstol)
}

kosmic_bootstrap_bridge <- function(data, decimals, replicates,
                                    t1min, t1max, t2min, t2max,
                                    sd_guess, abstol) {
  if(!is.numeric(data)) {
    abort("`data` must be a numeric vector.")
  }
  if(!is_bare_numeric(decimals, n = 1) | decimals < 0) {
    abort("`decimals` must be a single number >= 0.")
  }
  if(!is_bare_numeric(abstol, n = 1) | abstol <= 0) {
    abort("`abstol` must be a single number > 0.")
  }
  if(!is_bare_numeric(replicates, n = 1) | replicates < 1) {
    abort("`replicates` must be a single number >= 1.")
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
  
  # Get the starting seed
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  
  # Run the Kosmic algorithm on the original data
  impl_result <- kosmic_impl(data,
                             trunc(decimals),
                             trunc(replicates),
                             t1min, t1max,
                             t2min, t2max,
                             sd_guess, abstol)
  res <- impl_result$result
  # Leave the cost column out of the bootstrap estimates.
  bootstrap_estimates <- impl_result$boot[,-4]
  colnames(bootstrap_estimates) <- c("lambda", "mean", "sd", "t1", "t2")
  
  # Include a frequency table of the original data to allow plotting
  # later on
  freqs <- data.frame(result = round(data, trunc(decimals))) %>%
    count(result)
  
  # Create a kosmic object to hold the results
  new_kosmic_bootstrap(freqs,
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
             abstol = abstol,
             replicates = replicates,
             seed = seed,
             bootstrap_estimates = bootstrap_estimates,)
}

new_kosmic_bootstrap <- function(data, n,
                                 lambda, mean, sd, t1, t2,
                                 decimals, t1min, t1max, t2min, t2max, sd_guess, abstol,
                                 replicates, seed, bootstrap_estimates,
                                 class = character()) {
  if(!is_bare_numeric(replicates, n=1) | replicates <= 0) {
    abort("`replicates` must be a single positive integer.")
  }
  if(!is.matrix(bootstrap_estimates)) {
    abort("`bootstrap_estimates` must be a matrix.")
  }
  if(ncol(bootstrap_estimates) != 5) {
    abort("`bootstrap_extimates` must have five columns")
  }
  res <- new_kosmic(data, n,
                    lambda, mean, sd, t1, t2,
                    decimals, t1min, t1max, t2min, t2max, sd_guess, abstol,
                    class = c(class, "kosmic_bootstrap"))
  res$replicates <- replicates
  res$seed <- seed
  res$bootstrap_estimates <- bootstrap_estimates
  res
}

#' Resamples of a Kosmic Bootstrap
#'
#' @description This function takes an object created by \code{kosmic_bootstrap}
#'   and returns the resampled data it uses.
#'
#' @param object
#'
#' @return
#' @export
kosmic_resamples <- function(object) {
  # Save the current seed and restore the object's seed
  if (exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE))
    temp <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
  else temp<- NULL
  assign(".Random.seed",  object$seed, envir=.GlobalEnv)
  
  # Perform the resampling
  resamples <- kosmic_resamples_impl(results = object$data$result,
                                     counts = object$data$n,
                                     replicates = object$replicates,
                                     settings = object$settings)
  # Rearrange resamples into a data.frame
  res <- data.frame(id = rep(1:object$replicates, each = resamples$classes),
                    result = rep(resamples$result,
                                 times = object$replicates),
                    n = resamples$frequencies)
  
  # Restore the previous seed
  if (!is.null(temp)) assign(".Random.seed", temp, envir=.GlobalEnv)
  else rm(.Random.seed, pos=1)
  res
}
