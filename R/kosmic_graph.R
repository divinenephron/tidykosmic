#' Data to Easily Plot a Distributed Estimated by Kosmic
#'
#' @param k 
#'
#' @return
#' @export
kosmic_plot_data <- function(k) {
  data <- k$data
  m <- k$estimates["mean"]
  s <- k$estimates["sd"]
  l <- k$estimates["lambda"]
  t1 <- k$estimates["t1"]
  t2 <- k$estimates["t2"]
  decimals <- k$settings["decimals"]
  binwidth <- 10^-decimals
  
  lower_limit <- min(data$result, qboxcox(0.01, m, s, l))
  upper_limit <- max(data$result, qboxcox(0.99, m, s, l))
  
  list(lower_limit = lower_limit, upper_limit = upper_limit)
}



# observed <- hemoglobin
# k <- kosmic(observed$result, 1)
# 
# quantile(k)
# 

# 
# # Calculate the expected frequency of normal results.
# # This is:
# # probability density * total number of result in truncated region / total area under density curve in transformed region
# estimated <- data.frame(result = seq(from = min(observed$result), to = max(observed$result), by = binwidth)) %>%
#   mutate(result_transformed = kosmic:::boxcox(result, l),
#          density = dnorm(result_transformed, m, s))
# truncated_region_samples <- sum(observed$result > t1 & observed$result < t2)
# truncated_region_density <- sum(filter(estimated, result > t1, result < t2)$density)
# estimated <- estimated %>%
#   mutate(frequency = density * truncated_region_samples / truncated_region_density)
# 
# ggplot(observed, aes(result)) +
#   geom_histogram(binwidth = binwidth) +
#   geom_line(aes(x = result,
#                 y = frequency),
#             data = estimated,
#             color = "red")
