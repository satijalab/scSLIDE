#' Predict Responses from a PLS Model
#'
#' Performs response prediction on new data using a PLS model previously
#' fitted with \code{\link{RunPLS}(..., save.model = TRUE)}. Mirrors the
#' core behaviour of \code{predict.mvr()} from the \pkg{pls} package.
#'
#' @param object A \code{DimReduc} object or \code{Seurat} object containing
#'   a PLS model fitted with \code{save.model = TRUE}.
#' @param newdata New observations: a cells x features numeric matrix (for
#'   the DimReduc method) or a Seurat object (for the Seurat method).
#' @param reduction (Seurat method only) Name of the PLS reduction slot.
#'   Default \code{"pls"}.
#' @param ncomp Integer vector of component counts to predict with. Each
#'   value produces a separate set of predictions using the cumulative
#'   coefficients for that many components. Default: maximum available.
#' @param assay (Seurat method only) Assay to use from \code{newdata}.
#'   Default: \code{DefaultAssay(newdata)}.
#' @param layer (Seurat method only) Layer to use from the assay. Default
#'   \code{"data"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A matrix (n x q) if a single \code{ncomp} value, or a 3D array
#'   (n x q x length(ncomp)) if multiple values are given.
#'
#' @export
#' @rdname PredictPLS
#'
#' @importFrom SeuratObject Loadings DefaultAssay LayerData
#'
#' @examples
#' \dontrun{
#' # Fit PLS model with save.model = TRUE
#' obj <- RunPLS(obj, Y = "response", save.model = TRUE)
#'
#' # Predict on new Seurat object
#' predictions <- PredictPLS(obj, newdata = new_obj)
#' }
PredictPLS <- function(object, newdata, ...) {
  UseMethod("PredictPLS")
}

#' @rdname PredictPLS
#' @method PredictPLS DimReduc
#' @export
PredictPLS.DimReduc <- function(
    object,
    newdata,
    ncomp = NULL,
    ...
) {
  # Validate model exists

if (is.null(object@misc$model)) {
    stop("No saved model found in the DimReduc object. ",
         "Re-run RunPLS() with save.model = TRUE.")
  }

  # Extract model components
  coefficients <- object@misc$model$coefficients
  Xmeans <- object@misc$model$Xmeans
  Ymeans <- object@misc$model$Ymeans

  # Get feature names from loadings
  model.features <- rownames(Loadings(object))
  p <- dim(coefficients)[1]

  # Coerce newdata to matrix if needed
  if (!is.matrix(newdata)) {
    newdata <- as.matrix(newdata)
  }

  # Feature alignment
  if (length(model.features) > 0 && !is.null(colnames(newdata))) {
    missing.features <- setdiff(model.features, colnames(newdata))
    if (length(missing.features) > 0) {
      stop("newdata is missing ", length(missing.features),
           " features required by the model: ",
           paste(head(missing.features, 5), collapse = ", "),
           if (length(missing.features) > 5) ", ..." else "")
    }
    newdata <- newdata[, model.features, drop = FALSE]
  } else {
    if (ncol(newdata) != p) {
      stop("newdata has ", ncol(newdata),
           " columns but model expects ", p, " features")
    }
  }

  # Validate ncomp
  ncomp_max <- dim(coefficients)[3]
  if (is.null(ncomp)) {
    ncomp <- ncomp_max
  }
  if (any(ncomp < 1) || any(ncomp > ncomp_max)) {
    stop("ncomp values must be between 1 and ", ncomp_max)
  }

  nobs <- nrow(newdata)
  q <- length(Ymeans)

  if (length(ncomp) == 1) {
    B <- coefficients[, , ncomp, drop = FALSE]
    dim(B) <- dim(B)[1:2]
    B0 <- Ymeans - crossprod(Xmeans, B)
    pred <- newdata %*% B + rep(B0, each = nobs)
    # Set dimnames
    rownames(pred) <- rownames(newdata)
    colnames(pred) <- names(Ymeans)
  } else {
    pred <- array(NA_real_, dim = c(nobs, q, length(ncomp)))
    for (i in seq_along(ncomp)) {
      k <- ncomp[i]
      B <- coefficients[, , k, drop = FALSE]
      dim(B) <- dim(B)[1:2]
      B0 <- Ymeans - crossprod(Xmeans, B)
      pred[, , i] <- newdata %*% B + rep(B0, each = nobs)
    }
    # Set dimnames
    comp.names <- dimnames(coefficients)[[3]]
    if (!is.null(comp.names)) {
      comp.labels <- comp.names[ncomp]
    } else {
      comp.labels <- paste(ncomp, "comps")
    }
    dimnames(pred) <- list(
      rownames(newdata),
      names(Ymeans),
      comp.labels
    )
  }

  return(pred)
}

#' @rdname PredictPLS
#' @method PredictPLS Seurat
#' @export
PredictPLS.Seurat <- function(
    object,
    newdata,
    reduction = "pls",
    ncomp = NULL,
    assay = NULL,
    layer = "data",
    ...
) {
  # Extract the DimReduc object
  reduction.obj <- object[[reduction]]
  if (!inherits(reduction.obj, "DimReduc")) {
    stop("reduction '", reduction, "' not found or is not a DimReduc object")
  }

  # Validate model exists
  if (is.null(reduction.obj@misc$model)) {
    stop("No saved model found in reduction '", reduction, "'. ",
         "Re-run RunPLS() with save.model = TRUE.")
  }

  # Get model features
  features <- rownames(Loadings(reduction.obj))

  # Resolve assay/layer in newdata
  assay <- assay %||% DefaultAssay(newdata)

  # Extract data matrix for matching features
  available.features <- rownames(LayerData(newdata[[assay]], layer = layer))
  missing.features <- setdiff(features, available.features)
  if (length(missing.features) > 0) {
    stop("newdata is missing ", length(missing.features),
         " features required by the model: ",
         paste(head(missing.features, 5), collapse = ", "),
         if (length(missing.features) > 5) ", ..." else "")
  }

  new.mat <- t(as.matrix(
    LayerData(newdata[[assay]], layer = layer)[features, ]
  ))

  # Delegate to DimReduc method
  PredictPLS.DimReduc(
    object = reduction.obj,
    newdata = new.mat,
    ncomp = ncomp,
    ...
  )
}
