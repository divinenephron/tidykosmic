#' Data to Easily Plot a Distribution Estimated by Kosmic
#'
#' @param k 
#'
#' @return
#' @export
kosmic_plot_data <- function(k) {
  observed <- k$data
  m <- k$estimates["mean"]
  s <- k$estimates["sd"]
  l <- k$estimates["lambda"]
  t1 <- k$estimates["t1"]
  t2 <- k$estimates["t2"]
  decimals <- k$settings["decimals"]
  binwidth <- 10^-decimals
  
  # Generate the range of `result` values
  # They need to be rounded to avoid floatin poigt errors
  lower_limit <- min(observed$result)
  upper_limit <- max(observed$result)
  data <- tibble(
    result = round(seq(from = min(observed$result),
                       to = max(observed$result),
                       by = binwidth),
                   decimals)
  )
  
  # Add observed frequencies for each `result` value
  data <- data %>%
    left_join(observed, by = "result") %>%
    rename(observed.freq = n) %>%
    mutate(
      observed.freq = if_else(is.na(observed.freq), 0L, observed.freq)
    )
  
  # Add estimated frequencies for each `result` value.
  # The estimated frequency is fitted to the truncated region between t1 and t2,
  # so scale the area under the estimated frequencies so that they are equal to
  # the area under the observed frequencies.
  data <- data %>%
      mutate(
        result.transformed = boxcox(result, l),
        estimated.dens = dnorm(result.transformed, m, s)
      )
  truncated_area_observed <- sum(data[data$result >= t1 & data$result <= t2, "observed.freq"])
  truncated_area_estimated <- sum(data[data$result >= t1 & data$result <= t2, "estimated.dens"])
  data <- data %>%
    mutate(
      estimated.freq = estimated.dens * truncated_area_observed / truncated_area_estimated
    ) %>%
    select(-result.transformed, -estimated.dens)
  
  # Add cumulative frequencies for the estimate and observed frequencies.
  # The cumulative frequencies are fitted to the truncated area betwen t1 and t2,
  # so raise them up on the graph so they align at t1.
  below_t1_observed <- sum(data[data$result < t1, "observed.freq"])
  below_t1_estimated <- sum(data[data$result < t1, "estimated.freq"])
  data <- data %>%
    mutate(
      observed.cum = cumsum(observed.freq),
      estimated.cum = cumsum(estimated.freq) + below_t1_observed - below_t1_estimated
    )
  
  list(data = data)
}

# library(ggplot2)
# observed <- hemoglobin
# k <- kosmic(observed$result, 1)
# d <- kosmic_plot_data(k)
# ggplot(d$data) +
#   geom_bar(aes(x = result, y = observed.freq),
#            stat = "identity",
#            width = 0.1) +
#   geom_line(aes(x = result, y = estimated.freq),
#             color = "red")
# ggplot(d$data) +
#   geom_bar(aes(x = result, y = observed.cum),
#            stat = "identity",
#            width = 0.1) +
#   geom_line(aes(x = result, y = estimated.cum),
#             color = "red")
