#' @aliases kosmic-package
#' @keywords internal
#'
#' @useDynLib kosmic
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