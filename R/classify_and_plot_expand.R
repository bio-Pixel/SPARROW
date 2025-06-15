
#' Classify continuous values into binary labels using Gaussian mixture model
#'
#' @param score_vector A numeric vector (named or unnamed)
#' @return A binary vector (0 = low, 1 = high)
#' @export
classify_continuous_vector <- function(score_vector) {
  if (!is.numeric(score_vector)) {
    stop("Input must be a numeric vector")
  }

  model <- mclust::Mclust(score_vector, G = 2)
  labels <- model$classification
  group_means <- tapply(score_vector, labels, mean)
  positive_label <- which.max(group_means)
  binary_label <- as.integer(labels == positive_label)
  names(binary_label) <- names(score_vector)
  return(binary_label)
}
