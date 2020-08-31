# R doesn't allow us to modify package variables, so we'll create an environment
# to store the seed in.
seed_env <- new_environment()
seed_env$seed <- as.integer(Sys.time())

#' Random Number Generation for Kosmic
#'
#' @description To make your kosmic results reproduceable from run-to-run, set
#' the random number generator seed.
#'
#' @param seed An single positive integer.
#'
#' @export
#'
#' @examples
#' set_kosmic_seed(42)
set_kosmic_seed <- function(seed) {
  if (!is_bare_numeric(seed) | trunc(seed) < 0) {
    stop("`seed` must be a single positive integer.")
  }
  invisible(seed_env$seed <- trunc(seed))
}

get_kosmic_seed <- function() {
  seed_env$seed
}
