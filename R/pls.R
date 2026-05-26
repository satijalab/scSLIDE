# Resolve Y (or Y.add) into a numeric matrix.
#
# Accepts a numeric matrix (returned as-is), a bare numeric vector
# (promoted to single-column matrix), or a data.frame (converted via
# model.matrix(~. + 0, data = ...)).  Character vectors are NOT
# handled here; the Seurat method resolves those before calling this helper.
#
# @param Y Input response object.
# @param expected_nrow Expected number of rows (for validation); NULL skips.
# @param label Label used in error messages (e.g. "Y", "Y.add").
# @return A numeric matrix, or NULL if Y is NULL.
# @keywords internal
.resolve_Y_matrix <- function(Y, expected_nrow = NULL, label = "Y") {
  if (is.null(Y)) return(NULL)
  if (is.matrix(Y) && is.numeric(Y)) {
    Y_mat <- Y
  } else if (is.numeric(Y) && is.null(dim(Y))) {
    # bare numeric vector -> single-column matrix
    Y_mat <- matrix(Y, ncol = 1)
  } else if (is.data.frame(Y)) {
    Y_mat <- model.matrix(~. + 0, data = Y)
  } else {
    stop(sprintf(
      "'%s' must be a numeric matrix, a data.frame, or a character vector of metadata column names. Got: %s",
      label, paste(class(Y), collapse = ", ")
    ))
  }
  if (!is.null(expected_nrow) && nrow(Y_mat) != expected_nrow) {
    stop(sprintf(
      "nrow(%s) is %d but expected %d (number of cells/samples).",
      label, nrow(Y_mat), expected_nrow
    ))
  }
  Y_mat
}

# Compute R2 and RMSEP from raw Y and fitted values (no pls model object needed).
#
# @param Y n x q response matrix.
# @param fitted_values n x q x ncomp array of cumulative fitted values.
# @return A list with \code{$R2} and \code{$RMSEP}, each an \code{mvrVal}-class
#   object matching the structure of \code{pls::R2()} / \code{pls::RMSEP()}.
# @keywords internal
.compute_r2_rmsep <- function(Y, fitted_values) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  q <- ncol(Y)
  ncomp <- dim(fitted_values)[3]

  # SST per response (intercept-only model = column means)
  Y_mean <- matrix(colMeans(Y), nrow = n, ncol = q, byrow = TRUE)
  SST <- colSums((Y - Y_mean)^2)  # length q

  # R2$val layout from pls: 1 x (ncomp+1) x q
  # [1,1,] = intercept (0 comps), [1,a+1,] = a comps
  R2_val <- array(NA_real_, dim = c(1, ncomp + 1, q))
  RMSEP_val <- array(NA_real_, dim = c(1, ncomp + 1, q))

  # Intercept-only model
  R2_val[1, 1, ] <- 0
  RMSEP_val[1, 1, ] <- sqrt(SST / n)

  for (a in seq_len(ncomp)) {
    fitted_a <- fitted_values[, , a, drop = FALSE]
    dim(fitted_a) <- c(n, q)
    SSE <- colSums((Y - fitted_a)^2)
    R2_val[1, a + 1, ] <- 1 - SSE / SST
    RMSEP_val[1, a + 1, ] <- sqrt(SSE / n)
  }

  # Set dimnames matching pls convention
  comp_names <- c("(Intercept)", paste(seq_len(ncomp), "comps"))
  resp_names <- colnames(Y)
  if (is.null(resp_names)) resp_names <- paste0("Y", seq_len(q))

  dimnames(R2_val) <- list("train", comp_names, resp_names)
  dimnames(RMSEP_val) <- list("train", comp_names, resp_names)

  R2_out <- list(val = R2_val, type = "train")
  class(R2_out) <- "mvrVal"

  RMSEP_out <- list(val = RMSEP_val, type = "train")
  class(RMSEP_out) <- "mvrVal"

  list(R2 = R2_out, RMSEP = RMSEP_out)
}

#' Run Partial Least Squares (PLS) on Seurat Objects
#'
#' Performs Partial Least Squares regression analysis on single-cell data.
#' This function provides methods for plsr, spls, and cppls from the pls and spls packages.
#'
#' @param object An object to run PLS on
#' @param assay Name of Assay PLS is being run on
#' @param ncomp Number of components to compute
#' @param Y The response variable(s) for PLS regression.
#'   Accepted formats:
#'   \itemize{
#'     \item A character vector of column names in the Seurat object metadata
#'           (Seurat method only).
#'     \item A \code{data.frame}, which is converted via
#'           \code{model.matrix(~. + 0, data = Y)}.
#'     \item A numeric matrix (used as-is).
#'     \item A numeric vector (promoted to a single-column matrix).
#'   }
#'   The number of rows must match the number of cells/samples.
#' @param Y.add Additional response variable(s) containing relevant
#'   information about the observations. Only used for cppls. Accepts the same
#'   formats as \code{Y}.
#' @param pls.function PLS function from pls package to run (options: plsr, spls, cppls)
#' @param verbose Print progress messages, diagnostics (e.g. R2), and the top
#' genes associated with high/low loadings for the components
#' @param ndims.print components to print genes for
#' @param nfeatures.print Number of genes to print for each component
#' @param reduction.name the name of the DimReduc object
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names.
#' @param seed.use Set a random seed.  Setting NULL will not set a seed.
#' @param eta Thresholding parameter that controls the sparsity of the spls method (larger --> sparser). eta should be between 0 and 1.
#' @param features Features to compute PLS on
#' @param save.model Logical; if TRUE, save model components (coefficients, Xmeans, Ymeans, feature.mean, feature.sd) into the misc slot for later prediction. The per-gene \code{feature.mean} and \code{feature.sd} are computed from the training \code{data} layer so that \code{PredictPLS} can scale new data into the same coordinate system. Default FALSE.
#' @param threads Integer; number of threads for BPCells parallel computation
#'   (IterableMatrix path only). \code{1} (default) = single-threaded;
#'   \code{0} = auto-detect via \code{parallel::detectCores(logical = FALSE)};
#'   values \eqn{\ge 2} use that many threads. Ignored for in-memory data.
#' @param layer The layer in `assay` to use when running PLS analysis.
#' @param ... Additional arguments to be passed to the PLS function
#'
#' @return Returns a DimReduc object with PLS results
#' @export
#' @concept dimensional_reduction
#'
#' @examples
#' \dontrun{
#' # Run PLS on a Seurat object using metadata column names
#' seurat_obj <- RunPLS(seurat_obj, Y = "condition", ncomp = 10)
#'
#' # Run PLS with a pre-built numeric design matrix
#' design_mat <- model.matrix(~ age + sex, data = seurat_obj[[]])
#' seurat_obj <- RunPLS(seurat_obj, Y = design_mat, ncomp = 10)
#' }
#'
#' @import pls Matrix
#' @importFrom spls spls
#' @importFrom stats coef model.matrix median sd uniroot
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
    save.model = FALSE,
    threads = 1L,
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

  Y_mat <- .resolve_Y_matrix(Y, expected_nrow = ncol(object), label = "Y")
  # run PLS using the selected method
  if (pls.function == "cppls" & !is.null(Y.add)){
    requireNamespace("pls", quietly = FALSE)
    message("Fitting Y.add in CPPLS...")
    Y.add <- .resolve_Y_matrix(Y.add, expected_nrow = ncol(object), label = "Y.add")
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
    if (save.model) {
      warning("save.model is not supported for spls; model components will not be saved.")
    }
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
    if (save.model) {
      misc$model <- list(
        coefficients = pls.results$coefficients,
        Xmeans = pls.results$Xmeans,
        Ymeans = pls.results$Ymeans
      )
    }
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
#' @method RunPLS IterableMatrix
#' @export
RunPLS.IterableMatrix <- function(
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
    save.model = FALSE,
    threads = 1L,
    ...
) {
  pls.function <- match.arg(arg = pls.function)

  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }

  n.samples <- ncol(object)
  n.features <- nrow(object)
  max.ncomp <- min(n.samples, n.features) - 1L
  if (max.ncomp < 1L) {
    rlang::abort("IterableMatrix must have at least 2 samples/features to run PLS.")
  }
  ncomp <- min(ncomp, max.ncomp)

  Y_mat <- .resolve_Y_matrix(Y, expected_nrow = n.samples, label = "Y")

  if (pls.function == "cppls") {
    pls.results <- cppls_ondisk(X = object, Y = Y_mat, Y.add = Y.add,
                                ncomp = ncomp, threads = threads, ...)
  } else if (pls.function == "plsr") {
    if (!is.null(Y.add)) {
      warning("Y.add is not supported for plsr/kernelpls and will be ignored")
    }
    pls.results <- kernelpls(X = object, Y = Y_mat, ncomp = ncomp,
                             threads = threads, ...)
  } else if (pls.function == "spls") {
    pls.results <- spls_ondisk(X = object, Y = Y_mat, ncomp = ncomp,
                                eta = eta, threads = threads, ...)
  } else {
    rlang::abort("Unsupported pls.function for IterableMatrix")
  }

  if (pls.function == "spls") {
    # SPLS result extraction (mirrors RunPLS.default spls branch)
    feature.loadings <- pls.results$projection
    colnames(feature.loadings) <- paste0(reduction.key, seq_len(ncol(feature.loadings)))
    feature.names <- rownames(object)
    if (!is.null(feature.names)) {
      rownames(feature.loadings) <- feature.names[pls.results$A]
    }

    # Cell embeddings: use scores from the final PLS sub-fit (already computed
    # on-disk by kernelpls).  These are t(X[A,]_centered) %*% projection.
    cell.embeddings <- unclass(pls.results$scores)
    colnames(cell.embeddings) <- paste0(reduction.key, seq_len(ncol(cell.embeddings)))
    sample.names <- colnames(object)
    if (!is.null(sample.names) && length(sample.names) == nrow(cell.embeddings)) {
      rownames(cell.embeddings) <- sample.names
    }

    stdev <- numeric()
    r2_rmsep <- .compute_r2_rmsep(Y_mat, pls.results$fitted.values)
    misc <- list(R2 = r2_rmsep$R2,
                 RMSEP = r2_rmsep$RMSEP)
    R2 <- matrix(r2_rmsep$R2$val, byrow = TRUE, ncol = dim(r2_rmsep$R2$val)[2])
    mean_final_R2 <- mean(R2[nrow(R2), ])
    message("The average R2 of the PLS model is ", mean_final_R2)
    if (save.model) {
      warning("save.model is not supported for on-disk spls; model components will not be saved.")
    }
  } else {
    feature.loadings <- unclass(pls.results$projection)
    colnames(feature.loadings) <- paste0(reduction.key, seq_len(ncol(feature.loadings)))
    feature.names <- rownames(object)
    if (!is.null(feature.names) && length(feature.names) == nrow(feature.loadings)) {
      rownames(feature.loadings) <- feature.names
    }

    cell.embeddings <- unclass(pls.results$scores)
    colnames(cell.embeddings) <- paste0(reduction.key, seq_len(ncol(cell.embeddings)))
    sample.names <- colnames(object)
    if (!is.null(sample.names) && length(sample.names) == nrow(cell.embeddings)) {
      rownames(cell.embeddings) <- sample.names
    }

    stdev <- pls.results$Xvar
    r2_rmsep <- .compute_r2_rmsep(Y_mat, pls.results$fitted.values)
    misc <- list(R2 = r2_rmsep$R2,
                 RMSEP = r2_rmsep$RMSEP)
    # diagnostic message matching RunPLS.default()
    R2 <- matrix(r2_rmsep$R2$val, byrow = TRUE, ncol = dim(r2_rmsep$R2$val)[2])
    mean_final_R2 <- mean(R2[nrow(R2), ])
    message("The average R2 of the PLS model is ", mean_final_R2)
    if (save.model) {
      misc$model <- list(
        coefficients = pls.results$coefficients,
        Xmeans = pls.results$Xmeans,
        Ymeans = pls.results$Ymeans
      )
    }
  }

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
    save.model = FALSE,
    threads = 1L,
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
    save.model = save.model,
    threads = threads,
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
    save.model = FALSE,
    threads = 1L,
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
    save.model = save.model,
    threads = threads,
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
    save.model = FALSE,
    threads = 1L,
    ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  pls.function <- match.arg(arg = pls.function)

  # Resolve Y: character -> metadata lookup; matrix/data.frame -> use directly
  cells <- Cells(object[[assay]])
  if (is.character(Y)) {
    Y <- object[[Y]][cells, , drop = FALSE]
  } else if (is.matrix(Y) || is.data.frame(Y) || (is.numeric(Y) && is.null(dim(Y)))) {
    Y <- .resolve_Y_matrix(Y, expected_nrow = length(cells), label = "Y")
    # If rownames exist, reorder to match assay cell order
    if (!is.null(rownames(Y))) {
      missing <- setdiff(cells, rownames(Y))
      if (length(missing) > 0) {
        stop(sprintf("Y is missing %d cells present in the assay.", length(missing)))
      }
      Y <- Y[cells, , drop = FALSE]
    }
  } else {
    stop("'Y' must be a character vector of metadata column names, a numeric matrix/vector, or a data.frame.")
  }

  # Resolve Y.add with the same logic
  if (!is.null(Y.add)) {
    if (is.character(Y.add)) {
      Y.add <- object[[Y.add]][cells, , drop = FALSE]
    } else if (is.matrix(Y.add) || is.data.frame(Y.add) || (is.numeric(Y.add) && is.null(dim(Y.add)))) {
      Y.add <- .resolve_Y_matrix(Y.add, expected_nrow = length(cells), label = "Y.add")
      if (!is.null(rownames(Y.add))) {
        missing <- setdiff(cells, rownames(Y.add))
        if (length(missing) > 0) {
          stop(sprintf("Y.add is missing %d cells present in the assay.", length(missing)))
        }
        Y.add <- Y.add[cells, , drop = FALSE]
      }
    } else {
      stop("'Y.add' must be a character vector of metadata column names, a numeric matrix/vector, or a data.frame.")
    }
  }

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
    save.model = save.model,
    threads = threads,
    ...
  )

  # Store per-gene scaling parameters from the training data layer so that
  # PredictPLS.Seurat can transform new (unscaled) data into training-scaled
  # space before applying PLS coefficients.
  if (save.model && !is.null(reduction.data@misc$model)) {
    pls.features <- rownames(Loadings(reduction.data))
    train.data <- LayerData(object[[assay]], layer = "data")
    train.data <- train.data[pls.features, , drop = FALSE]

    if (inherits(train.data, "IterableMatrix")) {
      stats <- BPCells::matrix_stats(train.data, row_stats = "variance")
      feat.mean <- stats$row_stats["mean", ]
      feat.sd <- sqrt(stats$row_stats["variance", ])
    } else {
      feat.mean <- Matrix::rowMeans(train.data)
      feat.sd <- apply(train.data, 1, sd)
    }

    # Edge cases: NA means -> 0, zero/NA SDs -> 1
    feat.mean[is.na(feat.mean)] <- 0
    feat.sd[is.na(feat.sd) | feat.sd == 0] <- 1

    names(feat.mean) <- pls.features
    names(feat.sd) <- pls.features
    reduction.data@misc$model$feature.mean <- feat.mean
    reduction.data@misc$model$feature.sd <- feat.sd
  }

  object[[reduction.name]] <- reduction.data
  object <- LogSeuratCommand(object = object)
  return(object)
}


# Univariate soft thresholding (port of spls::ust).
#
# Sets elements of b to zero when their absolute value is below
# eta * max(|b|), and shrinks survivors toward zero.
#
# @param b Numeric vector or single-column matrix.
# @param eta Thresholding parameter in [0, 1).
# @return A single-column matrix of thresholded values.
# @keywords internal
.ust <- function(b, eta) {
  b_ust <- matrix(0, length(b), 1)
  if (eta < 1) {
    valb <- abs(b) - eta * max(abs(b))
    b_ust[valb >= 0] <- valb[valb >= 0] * (sign(b))[valb >= 0]
  }
  b_ust
}

# Sparse PLS direction vector (port of spls:::spls.dv).
#
# Finds a sparse direction vector via soft thresholding and, for multivariate
# Y, iterative SVD/ridge updates.
#
# @param Z Cross-product t(X) %*% Y (p x q).
# @param eta Sparsity parameter in [0, 1).
# @param kappa Ridge parameter in (0, 0.5].
# @param eps Convergence tolerance.
# @param maxstep Maximum number of iterations.
# @return A single-column matrix of direction weights (p x 1).
# @keywords internal
.spls_dv <- function(Z, eta, kappa, eps, maxstep) {
  p <- nrow(Z)
  q <- ncol(Z)
  Znorm1 <- stats::median(abs(Z))
  if (Znorm1 == 0) return(matrix(0, p, 1))
  Z <- Z / Znorm1

  if (q == 1) {
    return(.ust(Z, eta))
  }

  # q > 1: multivariate response

  M <- Z %*% t(Z)
  dis <- 10
  i <- 1

  if (kappa == 0.5) {
    cc <- matrix(10, p, 1)
    cc_old <- cc
    while (dis > eps && i <= maxstep) {
      mcsvd <- svd(M %*% cc)
      a <- mcsvd$u %*% t(mcsvd$v)
      cc <- .ust(M %*% a, eta)
      dis <- max(abs(cc - cc_old))
      cc_old <- cc
      i <- i + 1
    }
  } else if (kappa > 0 && kappa < 0.5) {
    kappa2 <- (1 - kappa) / (1 - 2 * kappa)
    cc <- matrix(10, p, 1)
    cc_old <- cc
    h <- function(lambda) {
      alpha <- solve(M + lambda * diag(p)) %*% M %*% cc
      obj <- t(alpha) %*% alpha - 1 / kappa2^2
      return(obj)
    }
    if (h(eps) * h(1e+30) > 0) {
      while (h(eps) <= 1e+05) {
        M <- 2 * M
        cc <- 2 * cc
      }
    }
    while (dis > eps && i <= maxstep) {
      if (h(eps) * h(1e+30) > 0) {
        while (h(eps) <= 1e+05) {
          M <- 2 * M
          cc <- 2 * cc
        }
      }
      lambdas <- stats::uniroot(h, c(eps, 1e+30))$root
      a <- kappa2 * solve(M + lambdas * diag(p)) %*% M %*% cc
      cc <- .ust(M %*% a, eta)
      dis <- max(abs(cc - cc_old))
      cc_old <- cc
      i <- i + 1
    }
  }
  cc
}

# Compute t(X_centered) %*% M for an on-disk IterableMatrix X and in-memory M.
#
# Uses BPCells' %*% operator and adjusts for centering:
#   t(X - 1 %*% t(Xmeans)) %*% M  =  t(X) %*% M - Xmeans %*% colSums(M)
#
# @param X IterableMatrix (features x samples).
# @param M Dense matrix (samples x q).
# @param Xmeans Numeric vector of feature means (length p).
# @return Dense matrix (p x q).
# @keywords internal
.ondisk_crossprod <- function(X, M, Xmeans) {
  # t(X) %*% M is (features x samples) transposed to (samples x features),
  # then multiplied by M (samples x q) => features x q ... but we want
  # t(X_centered) %*% M  where X_centered is (features x samples) centered
  # so each row has its mean subtracted.
  # t(X_centered) %*% M = t(X) %*% M - Xmeans %*% t(colSums(M))
  # But BPCells: X is features x samples, so t(X) is samples x features.
  # We want (samples x features)' %*% M(samples x q) = features x q
  # That is:  X %*% M  (since X is features x samples and M is samples x q)
  raw <- as.matrix(X %*% M)
  # centering correction: subtract outer(Xmeans, colSums(M))
  raw - Xmeans %*% t(colSums(M))
}

#' Sparse Partial Least Squares for On-Disk Matrices
#'
#' Generic function for sparse PLS. The \code{IterableMatrix} method computes
#' sparse PLS for BPCells on-disk matrices using a hybrid R/C++ approach.
#'
#' @param X Predictor matrix (or an IterableMatrix).
#' @param ... Arguments passed to methods.
#' @export
spls_ondisk <- function(X, ...) UseMethod("spls_ondisk")

#' Sparse Partial Least Squares for IterableMatrix
#'
#' @description
#' Compute sparse PLS decomposition for large on-disk BPCells matrices.
#' Uses BPCells for the expensive cross-product operations on disk, while
#' the sparsity step (direction vector + soft thresholding) and the PLS
#' sub-regression on active columns run in pure R on small in-memory subsets.
#'
#' @param X An IterableMatrix (features x samples, the BPCells convention).
#' @param Y A matrix or vector of responses (samples x q).
#' @param ncomp Integer; number of SPLS components to compute.
#' @param eta Numeric; sparsity parameter in [0, 1). Larger values produce
#'   sparser solutions. Default 0.5.
#' @param kappa Numeric; ridge parameter in (0, 0.5]. Default 0.5.
#' @param fit Character; reserved for compatibility. The on-disk sub-regression
#'   always uses kernel PLS via \code{\link{kernelpls}}. Default
#'   \code{"kernelpls"}.
#' @param center Logical; whether to center X and Y (default \code{TRUE}).
#'   Note: scaling is not applied for on-disk SPLS.
#' @param eps Numeric; convergence tolerance for direction vector fitting.
#'   Default 1e-4.
#' @param maxstep Integer; maximum iterations for direction vector fitting.
#'   Default 100.
#' @param threads Integer; number of threads for BPCells parallel computation.
#'   Currently unused (reserved for future optimisation). Default 1L.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with components:
#' \describe{
#'   \item{betahat}{Regression coefficients (p x q matrix).}
#'   \item{A}{Active set: indices of selected features.}
#'   \item{projection}{Projection matrix from the final PLS sub-fit
#'     (|A| x ncomp_actual matrix).}
#'   \item{Xmeans}{Feature means (length p).}
#'   \item{Ymeans}{Response means (length q).}
#'   \item{fitted.values}{Fitted Y values (n x q x ncomp array).}
#'   \item{scores}{X-scores from the final PLS sub-fit.}
#' }
#'
#' @details
#' Implements the SPLS algorithm (Chun & Keles, 2010) adapted for on-disk
#' data. Only \code{"pls2"} deflation (Y-deflation) is supported.
#' The predictor matrix is not scaled (only centered) when on disk; a message
#' is printed to inform users.
#'
#' @references
#' Chun, H. and Keles, S. (2010) Sparse partial least squares regression for
#' simultaneous dimension reduction and variable selection.
#' \emph{Journal of the Royal Statistical Society: Series B}, \bold{72}, 3--25.
#'
#' @export
#' @method spls_ondisk IterableMatrix
spls_ondisk.IterableMatrix <- function(X, Y, ncomp, eta = 0.5, kappa = 0.5,
                                        fit = "kernelpls", center = TRUE,
                                        eps = 1e-4, maxstep = 100L,
                                        threads = 1L, ...) {
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    rlang::abort("Package 'BPCells' is required for spls_ondisk.IterableMatrix but is not installed.")
  }
  if (!requireNamespace("pls", quietly = TRUE)) {
    rlang::abort("Package 'pls' is required for spls_ondisk.IterableMatrix but is not installed.")
  }

  # --- input validation ---
  ncomp <- as.integer(ncomp)
  if (length(ncomp) != 1 || ncomp <= 0) {
    rlang::abort("ncomp must be a single positive integer")
  }
  if (!is.numeric(eta) || length(eta) != 1 || eta < 0 || eta >= 1) {
    rlang::abort("eta must be a numeric value in [0, 1)")
  }
  if (!is.numeric(kappa) || length(kappa) != 1 || kappa <= 0 || kappa > 0.5) {
    rlang::abort("kappa must be a numeric value in (0, 0.5]")
  }

  Y <- as.matrix(Y)
  if (!is.numeric(Y)) {
    rlang::abort("Y must be numeric")
  }

  n <- ncol(X)    # samples
  p <- nrow(X)    # features

  if (n != nrow(Y)) {
    rlang::abort(sprintf("Number of samples in X (%d) must match nrow(Y) (%d)", n, nrow(Y)))
  }
  q <- ncol(Y)

  # Cap ncomp
  max_ncomp <- min(n - 1L, p)
  if (ncomp > max_ncomp) {
    ncomp <- max_ncomp
    message(sprintf("ncomp capped to %d (min(n-1, p))", ncomp))
  }

  # On-disk: center only, no scaling
  message("On-disk SPLS: centering is applied but scaling is not. ",
          "Results may differ from in-memory spls::spls() which scales by default.")

  # --- compute means and center Y ---
  if (center) {
    stats <- BPCells::matrix_stats(X, row_stats = "mean")
    Xmeans <- stats$row_stats["mean", ]
    Ymeans <- as.numeric(colMeans(Y))
    Y_c <- scale(Y, center = Ymeans, scale = FALSE)
  } else {
    Xmeans <- rep(0, p)
    Ymeans <- rep(0, q)
    Y_c <- Y
  }

  # --- SPLS iteration ---
  ip <- seq_len(p)
  betahat <- matrix(0, p, q)
  Y1 <- Y_c   # deflated Y (centered)
  fitted_array <- array(0, dim = c(n, q, ncomp))
  projection <- NULL
  A <- integer(0)
  scores <- NULL

  for (k in seq_len(ncomp)) {
    # Step 1: compute Z = t(X_centered) %*% Y1  (on-disk pass)
    Z <- .ondisk_crossprod(X, Y1, Xmeans)

    # Step 2: sparse direction vector (pure R, in-memory)
    what <- .spls_dv(Z, eta, kappa, eps, maxstep)

    # Step 3: build active set
    A <- unique(ip[what != 0 | betahat[, 1] != 0])

    if (length(A) == 0) {
      warning(sprintf("No active features at component %d; stopping early.", k))
      ncomp <- k - 1L
      if (ncomp == 0) {
        rlang::abort("SPLS selected zero features at the first component. Try decreasing eta.")
      }
      fitted_array <- fitted_array[, , seq_len(ncomp), drop = FALSE]
      break
    }

    # Step 4: subset on-disk X to active features (stays on-disk)
    X_A <- X[A, , drop = FALSE]  # IterableMatrix, |A| x n

    # Step 5: run kernel-PLS on active subset
    # kernelpls() requires ncomp < min(n, |A|).  pls::plsr silently
    # caps at min(n-1, |A|), so we mirror that behaviour here.
    n_sub_comp <- min(k, length(A), min(n, length(A)) - 1L)

    if (n_sub_comp >= 1L) {
      # Normal path: use on-disk kernelpls
      kpls_fit <- kernelpls(X_A, Y_c, ncomp = n_sub_comp,
                             center = center, threads = threads)

      betahat <- matrix(0, p, q)
      betahat[A, ] <- kpls_fit$coefficients[, , n_sub_comp, drop = TRUE]
      projection <- kpls_fit$projection
      fitted_k <- kpls_fit$fitted.values[, , n_sub_comp, drop = TRUE]
      dim(fitted_k) <- c(n, q)
      scores <- kpls_fit$scores
    } else {
      # Edge case: |A| == 1 so on-disk kernelpls can't fit even 1 component.
      # Materialize the single row (negligible memory) and use pls::plsr.
      X_A_dense <- as.matrix(X_A)
      if (center) X_A_dense <- X_A_dense - Xmeans[A]
      plsfit <- pls::plsr(Y_c ~ t(X_A_dense), ncomp = 1L,
                            method = "kernelpls", scale = FALSE)
      betahat <- matrix(0, p, q)
      betahat[A, ] <- matrix(coef(plsfit), length(A), q)
      projection <- plsfit$projection
      fitted_k <- t(X_A_dense) %*% betahat[A, , drop = FALSE]
      dim(fitted_k) <- c(n, q)
      scores <- unclass(plsfit$scores)
    }

    # Step 7: deflate Y (pls2 deflation)
    Y1 <- Y_c - fitted_k

    # Store fitted values (cumulative)
    fitted_array[, , k] <- fitted_k
  }

  if (ncomp == 0) {
    rlang::abort("SPLS produced no components.")
  }

  list(
    betahat = betahat,
    A = A,
    projection = projection,
    Xmeans = Xmeans,
    Ymeans = Ymeans,
    fitted.values = fitted_array,
    scores = scores
  )
}


#' Kernel Partial Least Squares
#'
#' Generic function for kernel-PLS. The \code{IterableMatrix} method computes
#' kernel-PLS for BPCells on-disk matrices.
#'
#' @param X Predictor matrix (or an IterableMatrix).
#' @param ... Arguments passed to methods.
#' @export
kernelpls <- function(X, ...) UseMethod("kernelpls")

#' Kernel Partial Least Squares for IterableMatrix
#'
#' @description
#' Compute kernel-PLS decomposition for large on-disk BPCells matrices using
#' the algorithm of Dayal and MacGregor (1997).
#'
#' @param X An IterableMatrix (features x samples, the BPCells convention).
#'   The transpose is handled internally.
#' @param Y A matrix or vector of responses (samples x q).
#' @param ncomp Integer; number of PLS components to compute.
#' @param center Logical; whether to center X and Y (default \code{TRUE}).
#' @param stripped Logical; if \code{TRUE} return only coefficients, Xmeans, and
#'   Ymeans for speed (default \code{FALSE}).
#' @param threads Integer; number of threads for BPCells parallel computation.
#'   \code{1} (default) = single-threaded;
#'   \code{0} = auto-detect via \code{parallel::detectCores(logical = FALSE)};
#'   values \eqn{\ge 2} use that many threads.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with components:
#' \describe{
#'   \item{coefficients}{Regression coefficients (p x q x ncomp array).}
#'   \item{scores}{X-scores (n x ncomp matrix).}
#'   \item{loadings}{X-loadings (p x ncomp matrix).}
#'   \item{loading.weights}{Loading weights (p x ncomp matrix).}
#'   \item{Yscores}{Y-scores (n x ncomp matrix).}
#'   \item{Yloadings}{Y-loadings (q x ncomp matrix).}
#'   \item{projection}{Projection matrix (p x ncomp matrix).}
#'   \item{Xmeans}{Column means of X (length p).}
#'   \item{Ymeans}{Column means of Y (length q).}
#'   \item{fitted.values}{Fitted Y values (n x q x ncomp array).}
#'   \item{residuals}{Residuals (n x q x ncomp array).}
#'   \item{Xvar}{Variance explained in X per component (length ncomp).}
#'   \item{Xtotvar}{Total variance in X (scalar).}
#' }
#' If \code{stripped = TRUE}, only \code{coefficients}, \code{Xmeans}, and
#' \code{Ymeans} are returned.
#'
#' @details
#' Implements the kernel-PLS algorithm (Dayal & MacGregor, 1997) for BPCells
#' on-disk matrices. The large predictor matrix X is accessed from disk
#' column-by-column; the response matrix Y is kept in memory.
#'
#' @references
#' Dayal, B. S. and MacGregor, J. F. (1997) Improved PLS algorithms.
#' \emph{Journal of Chemometrics}, \bold{11}, 73--85.
#'
#' @export
#' @method kernelpls IterableMatrix
kernelpls.IterableMatrix <- function(X, Y, ncomp, center = TRUE,
                                      stripped = FALSE, threads = 1L, ...) {
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    rlang::abort("Package 'BPCells' is required for kernelpls.IterableMatrix but is not installed.")
  }

  if (!is.numeric(ncomp) || length(ncomp) != 1 || ncomp != as.integer(ncomp)) {
    rlang::abort("ncomp must be a single whole number")
  }
  ncomp <- as.integer(ncomp)
  if (ncomp <= 0) {
    rlang::abort("ncomp must be a positive integer")
  }

  threads <- as.integer(threads)
  if (threads < 0L) rlang::abort("threads must be a non-negative integer")

  # Convert factor/character Y to a dummy indicator matrix
  if (is.factor(Y) || is.character(Y)) {
    Y <- stats::model.matrix(~ 0 + factor(Y))
    colnames(Y) <- sub("^factor\\(Y\\)", "", colnames(Y))
  }

  Y <- as.matrix(Y)

  if (!is.numeric(Y)) {
    rlang::abort("Y must be numeric (or a factor/character vector for classification)")
  }

  # BPCells' nrow()/ncol() already account for @transpose, so they always
  # return the user-facing dimensions: nrow = features, ncol = samples.
  n <- ncol(X)
  p <- nrow(X)
  q <- ncol(Y)

  if (n != nrow(Y)) {
    rlang::abort(sprintf("Number of samples in X (%d) must match nrow(Y) (%d)", n, nrow(Y)))
  }
  if (ncomp >= min(n, p)) {
    rlang::abort(sprintf("ncomp must be less than min(n, p) = %d", min(n, p)))
  }

  # The transpose flag tells C++ whether to swap row/col semantics.
  # BPCells default (features x samples) needs transpose=TRUE so that C++
  # treats rows as features and columns as samples.
  is_transposed <- !X@transpose

  if (threads == 0L) {
    threads <- max(parallel::detectCores(logical = FALSE), 1L, na.rm = TRUE)
  }

  n_splits <- min(threads * 4L, n)

  parallel_split <- getFromNamespace("parallel_split", "BPCells")
  iterate_matrix <- getFromNamespace("iterate_matrix", "BPCells")
  it <- X %>%
    BPCells::convert_matrix_type("double") %>%
    parallel_split(threads, n_splits) %>%
    iterate_matrix()

  result <- kernelpls_cpp(it, Y, as.integer(ncomp), center, stripped, is_transposed)

  # Reshape 2D coefficient matrix to 3D array (p x q x ncomp)
  comp_names <- paste(seq_len(ncomp), "comps")
  result$coefficients <- array(result$coefficients, dim = c(p, q, ncomp),
                               dimnames = list(NULL, NULL, comp_names))

  if (!stripped) {
    # Reshape fitted.values and residuals to 3D arrays (n x q x ncomp)
    result$fitted.values <- array(result$fitted.values, dim = c(n, q, ncomp),
                                  dimnames = list(NULL, NULL, comp_names))
    result$residuals     <- array(result$residuals, dim = c(n, q, ncomp),
                                  dimnames = list(NULL, NULL, comp_names))

    # Add class attributes for compatibility with pls package conventions
    class(result$scores)  <- class(result$Yscores) <- "scores"
    class(result$loadings) <- class(result$loading.weights) <-
      class(result$Yloadings) <- "loadings"
  }

  result
}


#' CPLS (Canonical PLS) for On-Disk Matrices
#'
#' Generic function for CPLS. The \code{IterableMatrix} method computes
#' CPLS for BPCells on-disk matrices.
#'
#' @param X Predictor matrix (or an IterableMatrix).
#' @param ... Arguments passed to methods.
#' @export
cppls_ondisk <- function(X, ...) UseMethod("cppls_ondisk")

#' CPLS (Canonical PLS) for IterableMatrix
#'
#' @description
#' Compute CPLS decomposition for large on-disk BPCells matrices.
#' CPLS uses canonical correlation analysis (CCA) to determine loading weights,
#' which can improve performance with multi-response or mixed-type response data.
#'
#' @param X An IterableMatrix (features x samples, the BPCells convention).
#'   The transpose is handled internally.
#' @param Y A matrix or vector of primary responses (samples x q).
#' @param Y.add A matrix of additional responses (samples x q_add), or NULL.
#'   Y.add provides information that guides loading weights but does not
#'   enter the regression coefficients.
#' @param ncomp Integer; number of PLS components to compute.
#' @param center Logical; whether to center X and Y (default \code{TRUE}).
#' @param stripped Logical; if \code{TRUE} return only coefficients, Xmeans, and
#'   Ymeans for speed (default \code{FALSE}).
#' @param w.tol Numeric; threshold for zeroing small loading weights (default 0).
#' @param X.tol Numeric; threshold for small-norm variable detection (default 1e-12).
#' @param threads Integer; number of threads for BPCells parallel computation.
#'   \code{1} (default) = single-threaded;
#'   \code{0} = auto-detect via \code{parallel::detectCores(logical = FALSE)};
#'   values \eqn{\ge 2} use that many threads.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with components:
#' \describe{
#'   \item{coefficients}{Regression coefficients (p x q x ncomp array).}
#'   \item{scores}{X-scores (n x ncomp matrix).}
#'   \item{loadings}{X-loadings (p x ncomp matrix).}
#'   \item{loading.weights}{Loading weights (p x ncomp matrix).}
#'   \item{Yscores}{Y-scores (n x ncomp matrix).}
#'   \item{Yloadings}{Y-loadings (q x ncomp matrix).}
#'   \item{projection}{Projection matrix (p x ncomp matrix).}
#'   \item{Xmeans}{Column means of X (length p).}
#'   \item{Ymeans}{Column means of Y (length q).}
#'   \item{fitted.values}{Fitted Y values (n x q x ncomp array).}
#'   \item{residuals}{Residuals (n x q x ncomp array).}
#'   \item{Xvar}{Variance explained in X per component (length ncomp).}
#'   \item{Xtotvar}{Total variance in X (scalar).}
#'   \item{gammas}{Power values per component (0.5 for all in this implementation).}
#'   \item{canonical.correlations}{Squared canonical correlations per component.}
#'   \item{A}{CCA weight vectors (q_full x ncomp matrix).}
#'   \item{smallNorm}{Indices of near-zero-norm variables.}
#' }
#' If \code{stripped = TRUE}, only \code{coefficients}, \code{Xmeans},
#' \code{Ymeans}, \code{gammas}, \code{canonical.correlations}, \code{A},
#' and \code{smallNorm} are returned.
#'
#' @export
#' @method cppls_ondisk IterableMatrix
cppls_ondisk.IterableMatrix <- function(X, Y, Y.add = NULL, ncomp,
                                         center = TRUE, stripped = FALSE,
                                         w.tol = 0, X.tol = 1e-12,
                                         threads = 1L, ...) {
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    rlang::abort("Package 'BPCells' is required for cppls_ondisk.IterableMatrix but is not installed.")
  }

  if (!is.numeric(ncomp) || length(ncomp) != 1 || ncomp != as.integer(ncomp)) {
    rlang::abort("ncomp must be a single whole number")
  }
  ncomp <- as.integer(ncomp)
  if (ncomp <= 0) {
    rlang::abort("ncomp must be a positive integer")
  }

  threads <- as.integer(threads)
  if (threads < 0L) rlang::abort("threads must be a non-negative integer")

  # Convert factor/character Y to a dummy indicator matrix
  if (is.factor(Y) || is.character(Y)) {
    Y <- stats::model.matrix(~ 0 + factor(Y))
    colnames(Y) <- sub("^factor\\(Y\\)", "", colnames(Y))
  }

  Y <- as.matrix(Y)

  if (!is.numeric(Y)) {
    rlang::abort("Y must be numeric (or a factor/character vector for classification)")
  }

  # Handle Y.add
  if (is.null(Y.add)) {
    Y_add_mat <- matrix(0.0, nrow = nrow(Y), ncol = 0)
  } else {
    if (is.factor(Y.add) || is.character(Y.add)) {
      Y.add <- stats::model.matrix(~ 0 + factor(Y.add))
      colnames(Y.add) <- sub("^factor\\(Y\\.add\\)", "", colnames(Y.add))
    }
    Y_add_mat <- as.matrix(Y.add)
    if (!is.numeric(Y_add_mat)) {
      rlang::abort("Y.add must be numeric (or a factor/character vector)")
    }
    if (nrow(Y_add_mat) != nrow(Y)) {
      rlang::abort(sprintf("nrow(Y.add) (%d) must match nrow(Y) (%d)",
                           nrow(Y_add_mat), nrow(Y)))
    }
  }

  # BPCells' nrow()/ncol() already account for @transpose, so they always
  # return the user-facing dimensions: nrow = features, ncol = samples.
  n <- ncol(X)
  p <- nrow(X)
  q <- ncol(Y)

  if (n != nrow(Y)) {
    rlang::abort(sprintf("Number of samples in X (%d) must match nrow(Y) (%d)", n, nrow(Y)))
  }
  if (ncomp >= min(n, p)) {
    rlang::abort(sprintf("ncomp must be less than min(n, p) = %d", min(n, p)))
  }

  is_transposed <- !X@transpose

  if (threads == 0L) {
    threads <- max(parallel::detectCores(logical = FALSE), 1L, na.rm = TRUE)
  }

  n_splits <- min(threads * 4L, n)

  parallel_split <- getFromNamespace("parallel_split", "BPCells")
  iterate_matrix <- getFromNamespace("iterate_matrix", "BPCells")
  it <- X %>%
    BPCells::convert_matrix_type("double") %>%
    parallel_split(threads, n_splits) %>%
    iterate_matrix()

  result <- cppls_cpp(it, Y, Y_add_mat, as.integer(ncomp), center, stripped,
                       is_transposed, as.double(w.tol), as.double(X.tol))

  # Reshape 2D coefficient matrix to 3D array (p x q x ncomp)
  comp_names <- paste(seq_len(ncomp), "comps")
  result$coefficients <- array(result$coefficients, dim = c(p, q, ncomp),
                               dimnames = list(NULL, NULL, comp_names))

  if (!stripped) {
    # Reshape fitted.values and residuals to 3D arrays (n x q x ncomp)
    result$fitted.values <- array(result$fitted.values, dim = c(n, q, ncomp),
                                  dimnames = list(NULL, NULL, comp_names))
    result$residuals     <- array(result$residuals, dim = c(n, q, ncomp),
                                  dimnames = list(NULL, NULL, comp_names))

    # Add class attributes for compatibility with pls package conventions
    class(result$scores)  <- class(result$Yscores) <- "scores"
    class(result$loadings) <- class(result$loading.weights) <-
      class(result$Yloadings) <- "loadings"
  }

  result
}
