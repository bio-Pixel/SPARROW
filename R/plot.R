
#' @import ggplot2

concave_dat <- function(coords, concavity = 1) {
  if (!requireNamespace("concaveman", quietly = TRUE)) install.packages("concaveman")
  if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
  df <- coords
  colnames(df) <- c("x", "y")
  hull_df <- as.data.frame(concaveman::concaveman(as.matrix(df), concavity = concavity))
  colnames(hull_df) <- c("V1", "V2")
  return(hull_df)
}

#' Plot cell type probability (continuous color)
#'
#' @param object CellProxobj
#' @param celltype Character, cell type to plot
#' @param outline Logical, whether to draw concave outline
#' @return ggplot object
#' @export
plotCellTypeProb <- function(object, celltype, outline = FALSE, pt.size = 1) {
  if (!celltype %in% object@celltypes) {
    stop(paste("Cell type", celltype, "not found in object@celltypes"))
  }
  df <- as.data.frame(object@meta.data)
  df$prob <- object@cell_prob[, celltype]

  p <- ggplot(df, aes(x = row, y = col, color = prob)) +
    geom_point(size = pt.size) +
    scale_color_gradientn(colours = c("gray90","#D9AAD7","#A765B1","#A765B1")) +
    theme_void() +
    ggtitle(paste("Cell type probability:", celltype)) +
    theme(legend.position = "right")

  if (outline) {
    hull_df <- concave_dat(coords = object@coords)
    p <- p + geom_polygon(data = hull_df, aes(x = V1, y = V2), color = "black", alpha = 0)
  }
  return(p)
}

#' Plot binary presence of a cell type
#'
#' @param object CellProxobj
#' @param celltype Character, cell type to plot
#' @param outline Logical, whether to draw concave outline
#' @return ggplot object
#' @export
plotCellType <- function(object, celltype, outline = FALSE, pt.size = 1, color = "red") {
  bin_col <- paste0("bin_", celltype)
  if (!bin_col %in% colnames(object@meta.data)) {
    stop(paste("Binarized column", bin_col, "not found in meta.data"))
  }
  df <- as.data.frame(object@meta.data)
  df$binary <- factor(df[[bin_col]], levels = c(0, 1), labels = c("absent", "present"))

  p <- ggplot(df, aes(x = row, y = col, color = binary)) +
    geom_point(size = pt.size) +
    scale_color_manual(values = c("grey80", color)) +
    theme_void() +
    ggtitle(paste("Cell type presence:", celltype)) +
    theme(legend.position = "right")

  if (outline) {
    hull_df <- concave_dat(coords = object@coords)
    p <- p + geom_polygon(data = hull_df, aes(x = V1, y = V2), color = "black", alpha = 0)
  }
  return(p)
}
