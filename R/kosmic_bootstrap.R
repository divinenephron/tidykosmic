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
                                   threads=1L,
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
                sd_guess, abstol, threads)
}

kosmic_bootstrap_bridge <- function(data, decimals, replicates,
                                    t1min, t1max, t2min, t2max,
                                    sd_guess, abstol, threads) {
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
  
  bootstrap_seed <- get_kosmic_seed()
  # Run the Kosmic algorithm on the original data
  impl_result <- kosmic_impl(data,
                             trunc(decimals),
                             trunc(replicates),
                             bootstrap_seed,
                             trunc(threads),
                             t1min, t1max,
                             t2min, t2max,
                             sd_guess, abstol)
  res <- impl_result$result
  # Leave the cost column out of the bootstrap estimates.
  # Also transpose it so that each estimates is a columns and each
  # replicate is a row.
  bootstrap_estimates <- t(impl_result$boot[,-4])
  colnames(bootstrap_estimates) <- c("lambda", "mean", "sd", "t1", "t2")
  
  # Include a frequency table of the original data to allow plotting
  # later on
  freqs <- data.frame(result = round(data, trunc(decimals))) %>%
    mutate(result = round(result, trunc(decimals))) %>%
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
             bootstrap_estimates = bootstrap_estimates,
             data_resamples = NULL)
}

new_kosmic_bootstrap <- function(data, n,
                                 lambda, mean, sd, t1, t2,
                                 decimals, t1min, t1max, t2min, t2max, sd_guess, abstol,
                                 replicates, bootstrap_estimates, data_resamples,
                                 ...,
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
  new_kosmic(data, n,
             lambda, mean, sd, t1, t2,
             decimals, t1min, t1max, t2min, t2max, sd_guess, abstol,
             replicates, bootstrap_estimates, data_resamples,
             ...,
             class = c(class, "kosmic_bootstrap"))
}