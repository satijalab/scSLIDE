#' Find Multi-Modal Nearest Neighbors
#'
#' A wrapper function that finds Multi-Modal weights and generate a Multi-Modal NN object.
#'
#' This function is essentially a wrapper function that perform the following 2 procedures:
#' \itemize{
#'   \item Seurat:::FindModalityWeights
#'   \item Seurat:::MultiModalNN
#' }
#'
#' @param object A Seurat object
#' @param sketch.assay The name of the sketched assay of the landmark cells.
#' @param reduction.list A list of dimensional reductions, one for each modality
#' @param dims.list A list containing the dimensions for each reduction to use
#' @param k.nn The number of nearest neighbors to compute for each modality
#' @param knn.range Range parameter for nearest neighbor search
#' @param l2.norm Perform L2 normalization on the cell embeddings
#' @param fix.wnn.weights Pre-specified modality weights. If provided, skips the calculation
#'   and uses these weights directly. Should be a list with the same length as reduction.list.
#' @param weighted.nn.name Multimodal neighbor object name
#' @param verbose Print progress bars and output
#' @param ...	Arguments passed to other methods
#' 
#' @return return a Seurat object that contains a weighted.nn Neighbor object between the landmark cells and all the other cells.
#' @export
#' @concept multimodal
#'
FindmmNN <- function(
    object,
    sketch.assay = NULL,
    reduction.list,
    dims.list,
    k.nn = 20,
    knn.range = 200,
    l2.norm = TRUE,
    fix.wnn.weights = NULL,
    weighted.nn.name = "weighted.nn",
    verbose = TRUE,
    ...
){
  # Get internal functions from Seurat
  FindModalityWeights <- getFromNamespace("FindModalityWeights", "Seurat")
  MultiModalNN <- getFromNamespace("MultiModalNN", "Seurat")
  
  if(k.nn <= 2) {
    stop("k.nn must be a integer larger than 1.")
  }
  # calculate WNN weights
  wnn.weight <- FindModalityWeights(object = object,
                                    reduction.list = reduction.list,
                                    dims.list = dims.list,
                                    fix.wnn.weights = fix.wnn.weights,
                                    snn.far.nn = FALSE,
                                    l2.norm = l2.norm,
                                    k.nn = ifelse(20 >= k.nn, 20, k.nn),
                                    ...)

  # prepare for the objects used to search for WNN
  if (is.null(sketch.assay) | !sketch.assay %in% names(object@assays)){
    stop("sketch.assay is not correctly defined. please check.")
  }
  sketch_cells <- Cells(object[[sketch.assay]])
  sketch_obj <- subset(object, cells = sketch_cells)
  # query_obj <- subset(object, cells = sketch_cells, invert = T)
  # run WNN
  weighted.nn <- MultiModalNN(object = sketch_obj,
                              query = object,
                              k.nn = ifelse(20 >= k.nn, 20, k.nn),
                              modality.weight = wnn.weight,
                              knn.range = knn.range,
                              verbose = T,
                              l2.norm = l2.norm,
                              ...)
  #
  object[[weighted.nn.name]] <- weighted.nn
  return (object)
}