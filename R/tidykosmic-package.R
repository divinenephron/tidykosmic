#' Estimate Reference Intervals from Routinely Collected Laboratory Data
#' 
#' To learn more about Kosmic, start with [`vignette("tidykosmic")`].
#' 
#' @aliases tidykosmic-package
#' @keywords internal
#'
#' @useDynLib tidykosmic
#' @importFrom Rcpp sourceCpp
#' @importFrom glue glue
#' @importFrom stats quantile
#' @importFrom rlang abort
#' @import rlang
#' @import dplyr
"_PACKAGE"

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot