setMethod("show", "CellProxobj", function(object) {
  cat("An object of class 'CellProxobj'\n\n")
  cat("Number of spots/cells: ", nrow(object@coords), "\n")
  cat("Number of genes:       ", nrow(object@expr), "\n")
  cat("Number of cell types:  ", length(object@celltypes), "\n")

  meta_cols <- colnames(object@meta.data)
  coord_cols <- colnames(object@coords)
  bin_cols <- grep("^bin_", meta_cols, value = TRUE)
  other_cols <- setdiff(setdiff(meta_cols, coord_cols), bin_cols)

  cat("Meta data columns:\n")
  cat("  • coords:    ", paste(coord_cols, collapse = ", "), "\n")
  if (length(bin_cols)) {
    cat("  • binarized: ", paste(bin_cols, collapse = ", "), "\n")
  }
  if (length(other_cols)) {
    cat("  • others:    ", paste(other_cols, collapse = ", "), "\n")
  }

  cat("\nUse @meta.data, @expr, @cell_prob, or accessor methods to explore.\n")
})
