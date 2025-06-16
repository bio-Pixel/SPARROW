
#' Plot multi-celltype probability with independent color gradients
#'
#' @param object CellProxobj
#' @param celltype Character vector of cell types to plot (default = all)
#' @param pt.size Numeric, point size
#' @param outline Logical, whether to draw concave outline
#' @param color Named vector to override celltype colors (default NULL)
#' @param coord.fixed Logical, whether to apply fixed aspect ratio (default TRUE)
#' @return ggplot object
#' @export
plotMultiCellTypeProb <- function(object, celltype = NULL, pt.size = 1, outline = TRUE, color = NULL, coord.fixed = TRUE) {
  prob <- object@cell_prob
  meta <- as.data.frame(object@meta.data)
  coords <- meta[, c("row", "col")]

  if (is.null(celltype)) {
    celltype <- colnames(prob)
  }
  if (!all(celltype %in% colnames(prob))) {
    stop("Some celltypes not found in cell_prob matrix.")
  }

  color_map <- if (!is.null(color)) color else object@params$celltype_colors
  base <- ggplot() +
    geom_point(data = coords, aes(x = row, y = col), color = "gray96")

  for (ct in celltype) {
    prob_vec <- prob[, ct]

    rng <- range(prob_vec, na.rm = TRUE)
    if (diff(rng) == 0) {
      prob_vec <- rep(0, length(prob_vec))
    } else {
      prob_vec <- (prob_vec - rng[1]) / diff(rng)
    }

    threshold <- find_main_valley_threshold(prob_vec)
    valid_idx <- which(prob_vec > threshold)

    if (length(valid_idx) > 0) {
      prob_sel <- prob_vec[valid_idx]
      rows <- coords[valid_idx, , drop = FALSE]
      df <- data.frame(row = rows$row, col = rows$col, value = prob_sel)

      base <- base +
        geom_point(data = df, aes(x = row, y = col, alpha = value), size = pt.size, color = color_map[ct]) +
        ggnewscale::new_scale_color()
    }
  }

  base <- base +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())

  if (coord.fixed) {
    base <- base + coord_fixed()
  }

  if (outline) {
    hull_df <- concave_dat(coords)
    base <- base + geom_polygon(data = hull_df, aes(x = V1, y = V2), color = "black", alpha = 0)
  }

  # Construct legend once using dummy data
  legend_df <- data.frame(
    row = seq_along(celltype),
    col = 1,
    celltype = factor(celltype, levels = celltype)
  )

  legend_layer <- ggplot(legend_df, aes(x = row, y = col, color = celltype)) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    geom_point(size = 3) +
    scale_color_manual(values = color_map[celltype]) +
    theme_void() +
    theme(legend.position = "right")

  legend_plot <- cowplot::get_legend(legend_layer)
  base <- cowplot::ggdraw() +
    cowplot::draw_plot(base, 0, 0, 0.85, 1) +
    cowplot::draw_plot(legend_plot, 0.85, 0, 0.15, 1)

  return(base)
}
