
setClass("SPARobj",
         slots = list(
           coords = "matrix",
           expr = "dgCMatrix",
           cell_prob = "matrix",
           meta.data = "DataFrame",
           celltypes = "character",
           params = "list"
         )
)

createSPARobj <- function(coords, expr, cell_prob, meta.data = NULL, params = list()) { 
  coords <- as.matrix(coords)
  colnames(coords) <- c("row", "col")
  expr <- methods::as(expr, "dgCMatrix")
  cell_prob <- as.matrix(cell_prob)

  # Check dimension consistency
  if (!all(rownames(coords) == colnames(expr))) {
    stop("Row names of coords and column names of expr must be identical and in the same order.")
  }
  if (!all(rownames(coords) == rownames(cell_prob))) {
    stop("Row names of coords and cell_prob must be identical and in the same order.")
  }

  # Initialize meta.data
  if (is.null(meta.data)) {
    meta.data <- S4Vectors::DataFrame(row.names = rownames(coords))
  } else if (!methods::is(meta.data, "DataFrame")) {
    meta.data <- S4Vectors::DataFrame(meta.data)
  }

  # Add binarized matrix
  binary_matrix <- apply(cell_prob, 2, classify_by_main_valley)
  colnames(binary_matrix) <- paste0("bin_", colnames(binary_matrix))

  # Combine into meta.data
  meta.data <- S4Vectors::cbind(
    S4Vectors::DataFrame(coords),
    meta.data,
    S4Vectors::DataFrame(binary_matrix)
  )

  # Assign default color palette
  color_50 <- c(
    "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#4a6fe3",
    "#b5bbe3", "#bec1d4", "#d6bcc0", "#bb7784", "#8dd593", "#f0b98d", "#f6c4e1",
    "#023fa5", "#7D58B9", "#11c638", "#F0E442", "#ef9708", "#ead3c6", "#FEC260",
    "#8e063b", "#d33f6a", "#e6afb9", "#0d6c0d", "#0fcfc0", "#7382BC", "#e07b91",
    "#F7C394", "#ade87c", "#2D81FF", "#FF6A00", "#00B37F", "#e07b91", "#A259FF",
    "#0099C6", "#FF33CC", "#66CC66", "#FFB347", "#9933FF", "#CC3366", "#33CCCC",
    "#FF6666", "#66B2FF", "#228B22", "#999933", "#CC99FF", "#FF99CC", "#669999",
    "#FFCC99", "#0066CC"
  )

  n_ct <- length(colnames(cell_prob))
  if (n_ct > length(color_50)) {
    extra_colors <- grDevices::rainbow(n_ct - length(color_50))
    all_colors <- c(color_50, extra_colors)
  } else {
    all_colors <- color_50
  }
  celltype_colors <- setNames(all_colors[seq_len(n_ct)], colnames(cell_prob))

  # Add to params
  params <- c(params, list(celltype_colors = celltype_colors))

  new("SPARobj", 
      coords = coords,
      expr = expr,
      cell_prob = cell_prob,
      meta.data = meta.data,
      celltypes = colnames(cell_prob),
      params = params
  )
}

# Accessor methods
setGeneric("getCoords", function(object) standardGeneric("getCoords"))
setMethod("getCoords", "SPARobj", function(object) object@coords)

setGeneric("getExpr", function(object) standardGeneric("getExpr"))
setMethod("getExpr", "SPARobj", function(object) object@expr)

setGeneric("getCellProb", function(object) standardGeneric("getCellProb"))
setMethod("getCellProb", "SPARobj", function(object) object@cell_prob)

setGeneric("getMeta", function(object) standardGeneric("getMeta"))
setMethod("getMeta", "SPARobj", function(object) object@meta.data)

setGeneric("getCelltypes", function(object) standardGeneric("getCelltypes"))
setMethod("getCelltypes", "SPARobj", function(object) object@celltypes)

setGeneric("getParams", function(object) standardGeneric("getParams"))
setMethod("getParams", "SPARobj", function(object) object@params)

# Generate a fixed palette of 100 distinct colors using viridis and hue
.generate_color_palette <- function(n = 100) {
  v <- viridisLite::viridis(n, option = "D")
  return(setNames(v, NULL))
}
