#' Convert Seurat spatial object to CellProxobj
#'
#' @param seurat_obj A Seurat object with spatial coordinates
#' @param cell_prob A matrix of cell type probabilities (rows = spots, cols = cell types)
#' @param expr_slot Expression slot to extract (default = "data")
#' @param assay Assay name (default = NULL for default assay)
#' @param params List of optional parameters to include
#' @return A CellProxobj object
#' @export
convertSeuratToCellProx <- function(seurat_obj,
                                    cell_prob,
                                    expr_slot = "data",
                                    assay = NULL,
                                    params = list()) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required.")
  }

  # Expression matrix
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }
  expr <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = expr_slot)

  # Coordinates
  coords <- Seurat::GetTissueCoordinates(seurat_obj)
  coords <- as.matrix(coords)
  colnames(coords) <- c("row", "col")

  # Cell probabilities check
  if (!all(rownames(cell_prob) %in% rownames(coords))) {
    stop("Row names of cell_prob must match those of Seurat spatial coordinates.")
  }

  cell_prob <- cell_prob[rownames(coords), , drop = FALSE]

  # Construct object
  object <- createCellProxobj(coords = coords,
                              expr = expr,
                              cell_prob = cell_prob,
                              meta.data = seurat_obj@meta.data,
                              params = params)
  return(object)
}
