find_main_valley_threshold <- function(score_vector) {
  dens <- density(score_vector)
  y_vals <- dens$y
  x_vals <- dens$x

  peaks_idx <- which(diff(sign(diff(y_vals))) == -2) + 1
  valleys_idx <- which(diff(sign(diff(y_vals))) == 2) + 1

  if (length(peaks_idx) < 2) {
    warning("Less than two peaks detected, returning median as threshold")
    return(median(score_vector))
  }

  peak_heights <- y_vals[peaks_idx]
  top_peaks_idx <- sort(peaks_idx[order(peak_heights, decreasing = TRUE)[1:2]])

  valleys_between <- valleys_idx[valleys_idx > top_peaks_idx[1] & valleys_idx < top_peaks_idx[2]]

  if (length(valleys_between) == 0) {
    warning("No valley between main peaks, returning median as threshold")
    return(median(score_vector))
  }

  valley_heights <- y_vals[valleys_between]
  main_valley_idx <- valleys_between[which.min(valley_heights)]

  threshold <- x_vals[main_valley_idx]
  return(threshold)
}

classify_by_main_valley <- function(score_vector, plot = FALSE) {
  threshold <- find_main_valley_threshold(score_vector)

  labels <- as.integer(score_vector > threshold)
  names(labels) <- names(score_vector)

  if (plot) {
    dens <- density(score_vector)
    plot(dens, main = "Density with main valley threshold", xlab = "Score")
    abline(v = threshold, col = "red", lwd = 2, lty = 2)
  }

  return(labels)
}
