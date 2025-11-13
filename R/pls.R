#' Run Partial Least Squares (PLS) on Seurat Objects
#'
#' Performs Partial Least Squares regression analysis on single-cell data.
#' This function provides methods for plsr, spls, and cppls from the pls and spls packages.
#'
#' @param object An object to run PLS on
#' @param assay Name of Assay PLS is being run on
#' @param ncomp Number of components to compute
#' @param Y a vector or matrix of responses, i.e., the dependent variable that PLS regresses X on. The length / number of rows should be the same as the number of cells.
#' Y can have multiple columns.
#' @param Y.add a vector or matrix of additional responses containing relevant information about the observations. Only used for cppls.
#' @param pls.function PLS function from pls package to run (options: plsr, spls, cppls)
#' @param verbose Print the top genes associated with high/low loadings for
#' the components
#' @param ndims.print components to print genes for
#' @param nfeatures.print Number of genes to print for each component
#' @param reduction.name the name of the DimReduc object
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names.
#' @param seed.use Set a random seed.  Setting NULL will not set a seed.
#' @param eta Thresholding parameter that controls the sparsity of the spls method (larger --> sparser). eta should be between 0 and 1.
#' @param features Features to compute PLS on
#' @param layer The layer in `assay` to use when running PLS analysis.
#' @param ... Additional arguments to be passed to the PLS function
#'
#' @return Returns a DimReduc object with PLS results
#' @export
#' @concept dimensional_reduction
#'
#' @examples
#' \dontrun{
#' # Run PLS on a Seurat object
#' seurat_obj <- RunPLS(seurat_obj, Y = "condition", ncomp = 10)
#' }
#'
#' @import pls
#' @importFrom spls spls
#' @importFrom stats model.matrix
#' @importFrom SeuratObject CreateDimReducObject DefaultAssay Assays Cells
#' @importFrom Seurat LogSeuratCommand
#' @importFrom utils getFromNamespace
#' @rdname RunPLS
#'
RunPLS <- function(object, ...) {
  UseMethod(generic = 'RunPLS', object = object)
}

#' @rdname RunPLS
#' @method RunPLS default
#' @export
RunPLS.default <- function(
    object,
    assay = NULL,
    ncomp = 20,
    Y = NULL,
    Y.add = NULL,
    pls.function = c("plsr", "spls", "cppls"),
    verbose = TRUE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.name = "pls",
    reduction.key = "PLS_",
    seed.use = 42,
    eta = 0.5,
    ...
) {
  # Get internal function from Seurat
  CheckDots <- getFromNamespace("CheckDots", "Seurat")
  
  CheckDots(..., fxns = pls.function)
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  ncomp <- min(ncomp, ncol(x = object))
  pls.fxn <- eval(expr = parse(text = pls.function))

  Y_mat <- model.matrix(~. + 0, data = Y)
  # run PLS using the selected method
  if (pls.function == "cppls" & !is.null(Y.add)){
    requireNamespace("pls", quietly = FALSE)
    message("Fitting Y.add in CPPLS...")
    Y.add <- model.matrix(~. + 0, data = Y.add)
    pls.results <- pls.fxn(Y_mat ~ t(object), ncomp = ncomp, validation = "none", Y.add = Y.add, ...)
  } else if (pls.function == "spls") {
    requireNamespace("spls", quietly = FALSE)
    pls.results <- pls.fxn(x = t(object), y = Y_mat, K = ncomp, fit = "kernelpls", eta = eta, ...)
  } else {
    requireNamespace("pls", quietly = FALSE)
    pls.results <- pls.fxn(Y_mat ~ t(object), ncomp = ncomp, validation = "none", ...)
  }
  #
  if (pls.function == "spls"){
    feature.loadings <- pls.results$projection
    colnames(feature.loadings) <- paste0(reduction.key, 1:ncol(feature.loadings))
    rownames(feature.loadings) <- rownames(object)[pls.results$A]
    cell.embeddings <- t(as.matrix(object[pls.results$A, ])) %*% feature.loadings
    colnames(cell.embeddings) <- paste0(reduction.key, 1:ncol(cell.embeddings))
    stdev <- numeric()
    misc <- list()
  } else {
    feature.loadings <- unclass(pls.results$projection)
    colnames(feature.loadings) <- paste0(reduction.key, 1:ncol(feature.loadings))
    rownames(feature.loadings) <- rownames(object)
    cell.embeddings <- unclass(pls.results$scores)
    colnames(cell.embeddings) <- paste0(reduction.key, 1:ncol(cell.embeddings))
    # if we are using plsr or cppls, we are able to retrieve the var_explained for each comp
    stdev <- pls.results$Xvar
    # other miscellaneous data might be useful for model selection
    misc <- list(R2 = pls::R2(pls.results),
                 RMSEP = pls::RMSEP(pls.results))
    # a diagnosis component
    R2 <- matrix(pls::R2(pls.results)$val, byrow = T, ncol = dim(pls::R2(pls.results)$val)[2])
    mean_final_R2 <- mean(R2[nrow(R2), ])
    message("The average R2 of the PLS model is ", mean_final_R2)
  }
  #
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    stdev = stdev,
    assay = assay,
    key = reduction.key,
    misc = misc
  )
  if (verbose) {
    print(x = reduction.data, dims = ndims.print, nfeatures = nfeatures.print)
  }
  return(reduction.data)
}

#' @rdname RunPLS
#' @method RunPLS Assay
#' @export
RunPLS.Assay <- function(
    object,
    assay = NULL,
    features = NULL,
    ncomp = 20,
    Y = NULL,
    Y.add = NULL,
    pls.function = c("plsr", "spls", "cppls"),
    verbose = TRUE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.name = "pls",
    reduction.key = "PLS_",
    seed.use = 42,
    eta = 0.5,
    ...
) {
  # Get internal function from Seurat
  PrepDR <- getFromNamespace("PrepDR", "Seurat")
  
  data.use <- PrepDR(
    object = object,
    features = features,
    verbose = verbose
  )
  reduction.data <- RunPLS(
    object = data.use,
    assay = assay,
    ncomp = ncomp,
    Y = Y,
    Y.add = Y.add,
    pls.function = pls.function,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    eta = eta,
    ...

  )
  return(reduction.data)
}

#' @rdname RunPLS
#' @method RunPLS StdAssay
#' @export
RunPLS.StdAssay <- function(
    object,
    assay = NULL,
    features = NULL,
    layer = 'scale.data',
    ncomp = 20,
    Y = NULL,
    Y.add = NULL,
    pls.function = c("plsr", "spls", "cppls"),
    verbose = TRUE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.name = "pls",
    reduction.key = "PLS_",
    seed.use = 42,
    eta = 0.5,
    ...
) {
  # Get internal function from Seurat
  PrepDR5 <- getFromNamespace("PrepDR5", "Seurat")
  
  data.use <- PrepDR5(
    object = object,
    features = features,
    layer = layer,
    verbose = verbose
  )
  reduction.data <- RunPLS(
    object = data.use,
    assay = assay,
    ncomp = ncomp,
    Y = Y,
    Y.add = Y.add,
    pls.function = pls.function,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    eta = eta,
    ...

  )
  return(reduction.data)
}

#' @rdname RunPLS
#' @method RunPLS Seurat
#' @export
RunPLS.Seurat <- function(
    object,
    assay = NULL,
    features = NULL,
    ncomp = 20,
    Y = NULL,
    Y.add = NULL,
    pls.function = c("plsr", "spls", "cppls"),
    verbose = TRUE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.name = "pls",
    reduction.key = "PLS_",
    seed.use = 42,
    eta = 0.5,
    ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  pls.function <- match.arg(arg = pls.function)

  #
  Y <- object[[Y]][Cells(object[[assay]]), , drop = F]
  if(!is.null(Y.add)) Y.add <- object[[Y.add]][Cells(object[[assay]]), , drop = F]

  reduction.data <- RunPLS(
    object = object[[assay]],
    assay = assay,
    features = features,
    ncomp = ncomp,
    Y = Y,
    Y.add = Y.add,
    pls.function = pls.function,
    verbose = verbose,
    ndims.print = ndims.print,
    nfeatures.print = nfeatures.print,
    reduction.key = reduction.key,
    seed.use = seed.use,
    eta = eta,
    ...
  )
  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}