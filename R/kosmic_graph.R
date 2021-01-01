#' Data to Easily Plot a Distribution Estimated by Kosmic
#'
#' @param k 
#'
#' @return
#' @export
kosmic_plot_data <- function(k, probs = c(0.025, 0.975)) {
  observed <- k$data
  m <- k$estimates[["mean"]]
  s <- k$estimates[["sd"]]
  l <- k$estimates[["lambda"]]
  t1 <- k$estimates[["t1"]]
  t2 <- k$estimates[["t2"]]
  decimals <- k$settings[["decimals"]]
  binwidth <- 10^-decimals
  
  # Generate the range of `result` values
  # They need to be rounded to avoid floating point errors
  lower_limit <- min(observed$result)
  upper_limit <- max(observed$result)
  data <- data.frame(
    result = round(seq(from = min(observed$result),
                       to = max(observed$result),
                       by = binwidth),
                   decimals)
  )
  
  # Add observed frequency and cumulative frequency
  data <- data %>%
    left_join(observed, by = "result") %>%
    rename(observed.freq = n) %>%
    mutate(
      observed.freq = if_else(is.na(observed.freq), 0L, observed.freq),
      observed.cum = cumsum(observed.freq)
    )
  
  # Add density (d) and cumulative distribution (p) of the estimated
  # distribution.
  # The Box-Cox transform means the area under the density function is not
  # one, so we bodge it by calculating the area then dividing through by it.
  f <- function(x) { dnorm(boxcox(x, l), m, s)}
  d_area <- integrate(f, 0, Inf)$value
  data <- data %>%
    mutate(
      estimated.d = dnorm(boxcox(result, l), m, s) / d_area,
      estimated.p = pnorm(boxcox(result, l), m, s)
    )
  
  # Add estimated frequencies
  # The estimated frequency is fitted to the truncated region (>= t1 and < t2),
  # so scale the area under the estimated frequencies so that they are equal to
  # the area under the observed frequencies.
  trunc_area_observed <- sum(data[data$result >= t1 & data$result < t2, "observed.freq"]) * binwidth
  trunc_area_estimated <- sum(data[data$result >= t1 & data$result < t2, "estimated.d"]) * binwidth
  scale_freq <- trunc_area_observed / trunc_area_estimated
  data <- data %>%
    mutate(
      estimated.freq = estimated.d * scale_freq
    )
  
  # Add cumulative frequencies for the estimate.
  # The estimate is fitted to the truncated region (>= t1 and < t2), so the height
  # of this part of the curve should be equal to the number of observed results
  # in the truncated region. The observed data also contains extra abnormal results
  # below the estimate, so raise the estimate up so that its cumulative frequency
  # is equal to the observed at t1.
  trunc_count_observed <- sum(data[data$result >= t1 & data$result < t2, "observed.freq"])
  scale_cum <- trunc_count_observed / trunc_area_estimated
  below_t1_observed <- sum(data[data$result < t1, "observed.freq"])
  below_t1_estimated <- sum(data[data$result < t1, "estimated.freq"])
  data <- data %>%
    mutate(
      estimated.cum =
        estimated.p * scale_cum + below_t1_observed - below_t1_estimated
    )
  
  quantiles <- data.frame(
    prob = probs,
    quantile = quantile(k, probs = probs, names = FALSE)
  )
  
  list(data = data,
       quantiles = quantiles)
}

#' Plot a Distribution Estimated by Kosmic
#' 
#' @param x 
#'
#' @param type 
#'
#' @rdname plot.kosmic
#' @importFrom graphics plot
#' @export
plot.kosmic <- function(x, type = c("frequency", "cumulative"),
                        probs = c(0.025, 0.975), ...) {
  print(autoplot(x, type, ...))
  invisible(x)
}

#' @rdname plot.kosmic
#' @importFrom ggplot2 ggplot geom_bar geom_line aes .data
#' @export
autoplot.kosmic <- function(object, type = c("frequency", "cumulative"),
                            probs = c(0.025, 0.975), ...) {
  type <- arg_match(type)
  if (!inherits(object, "kosmic")) {
    abort("Use only with `kosmic` objects")
  }
  
  plot_data <- kosmic_plot_data(object)
  binwidth <- 10^-object$settings[["decimals"]]
  
  g <- ggplot2::ggplot(plot_data$data) +
    ggplot2::theme_classic()
  if (type == "frequency") {
    g +
      ggplot2::geom_bar(aes(x = .data$result, y = .data$observed.freq),
                        stat = "identity",
                        width = binwidth,
                        fill = "#FC5A45") +
      ggplot2::geom_line(aes(x = .data$result, y = .data$estimated.freq),
                         color = "black") +
      ggplot2::geom_vline(aes(xintercept = quantile),
                          data = plot_data$quantiles,
                          linetype = "dashed")
  } else if (type == "cumulative") {
    g + 
      ggplot2::geom_bar(aes(x = .data$result, y = .data$observed.cum),
                        stat = "identity",
                        width = binwidth,
                        fill = "#FC5A45") +
      ggplot2::geom_line(aes(x = .data$result, y = .data$estimated.cum),
                         color = "black")
  }
}