#' Batch-specific Mean Correction
#'
#' Performs a simple batch-specific mean correction on a sample-level Seurat object.
#' For each batch, the mean vector of the control samples in that batch is computed
#' and subtracted from all samples (including controls) belonging to that batch.
#'
#' When a batch contains no control samples, the function falls back to the global
#' control mean (across all batches) and emits a warning.
#'
#' @param object A Seurat object (sample-level)
#' @param assay Name of the assay to correct. Default is `DefaultAssay(object)`.
#' @param layer Layer within the assay that contains the data to correct. Default is `"data"`.
#' @param batch_key Column name in `object@meta.data` identifying the batch (e.g., plate, pool).
#' @param condition_key Column name in `object@meta.data` identifying the experimental condition.
#' @param control_label Value(s) in the `condition_key` column that mark control samples.
#' @param verbose Print progress messages. Default `TRUE`.
#'
#' @return The Seurat object with the corrected data written back to the specified assay layer.
#' @export
#' @concept batch_correction
#'
#' @examples
#' \dontrun{
#' obj <- CorrectBatchMean(obj, assay = "LMC",
#'                         batch_key = "plate",
#'                         condition_key = "Condition",
#'                         control_label = "ctrl-inj")
#' }
#'
#' @importFrom SeuratObject LayerData LayerData<- DefaultAssay Cells
#'
CorrectBatchMean <- function(object,
                             assay = NULL,
                             layer = "data",
                             batch_key = NULL,
                             condition_key = NULL,
                             control_label = NULL,
                             verbose = TRUE) {

  # --- input validation ---
  if (!inherits(object, "Seurat")) {
    stop("'object' must be a Seurat object")
  }

  assay <- assay %||% DefaultAssay(object = object)
  if (!assay %in% names(object@assays)) {
    stop("Assay '", assay, "' not found in object")
  }

  if (is.null(batch_key)) {
    stop("Please specify 'batch_key': the metadata column identifying batches")
  }
  if (!batch_key %in% colnames(object@meta.data)) {
    stop("Column '", batch_key, "' not found in object metadata")
  }

  if (is.null(condition_key)) {
    stop("Please specify 'condition_key': the metadata column identifying conditions")
  }
  if (!condition_key %in% colnames(object@meta.data)) {
    stop("Column '", condition_key, "' not found in object metadata")
  }

  if (is.null(control_label)) {
    stop("Please specify 'control_label': the value(s) marking control samples in '", condition_key, "'")
  }

  # --- retrieve the data matrix ---
  data_mat <- LayerData(object = object[[assay]], layer = layer)

  all_cells <- Cells(object[[assay]])
  conditions <- object@meta.data[all_cells, condition_key]
  batches    <- object@meta.data[all_cells, batch_key]
  is_control <- conditions %in% control_label

  # pre-compute global control mean as fallback
  global_ctrl_idx <- all_cells[is_control]
  if (length(global_ctrl_idx) == 0) {
    stop("No control samples found with '", condition_key, "' in ",
         paste(control_label, collapse = ", "))
  }
  if (length(global_ctrl_idx) >= 2) {
    global_mean_vec <- rowMeans(data_mat[, global_ctrl_idx])
  } else {
    global_mean_vec <- data_mat[, global_ctrl_idx]
  }

  # --- batch-specific mean correction ---
  unique_batches <- unique(batches)

  for (b in unique_batches) {
    batch_ctrl_idx <- all_cells[batches == b & is_control]

    if (length(batch_ctrl_idx) >= 2) {
      mean_vec <- rowMeans(data_mat[, batch_ctrl_idx])
    } else if (length(batch_ctrl_idx) == 1) {
      mean_vec <- data_mat[, batch_ctrl_idx]
    } else {
      warning("Batch '", b, "' has no control samples; using global control mean as fallback")
      mean_vec <- global_mean_vec
    }

    batch_all_idx <- all_cells[batches == b]
    data_mat[, batch_all_idx] <- sweep(data_mat[, batch_all_idx], 1, mean_vec)

    if (isTRUE(verbose)) {
      message("Done with batch: ", b)
    }
  }

  # --- write corrected data back ---
  LayerData(object[[assay]], layer = layer) <- data_mat
  if (isTRUE(verbose)) {
    message("Batch mean correction completed for assay '", assay, "'")
  }

  return(object)
}


#' Shrinkage + Regression Batch-Mean Correction
#'
#' Performs a shrinkage-based batch-mean correction on a sample-level Seurat
#' object.
#'
#' For every gene the function estimates the between-batch variance
#' (\eqn{\tau^2}) and within-batch sampling variance (\eqn{\sigma^2}) of the
#' control means using a one-way random-effects model.
#' The batch-specific control mean is then shrunk towards the global control
#' mean:
#' \deqn{Z^*_{gb} = \lambda_{gb}\,Z_{gb} + (1-\lambda_{gb})\,\bar{Z}_g}
#' where
#' \eqn{\lambda_{gb}=\tau^2/(\tau^2 + \sigma^2/n_b)}.
#'
#' A per-sample regression of the expression profile on the shrunk batch mean
#' is used to estimate a sample-specific scaling factor \eqn{\beta_i}, and
#' the correction applied is
#' \deqn{x^{\text{corrected}}_{gi} = x_{gi} - \hat\beta_i\,Z^*_{gb_i}.}
#'
#' Optionally an attenuation-correction is applied to the OLS estimate of
#' \eqn{\beta} to account for measurement error in \eqn{Z^*}.
#'
#' @param object A Seurat object (sample-level).
#' @param assay Name of the assay to correct.
#'   Default is \code{DefaultAssay(object)}.
#' @param layer Layer within the assay that contains the data to correct.
#'   Default is \code{"data"}.
#' @param batch_key Column name in \code{object@@meta.data} identifying the
#'   batch (e.g., plate, pool).
#' @param condition_key Column name in \code{object@@meta.data} identifying the
#'   experimental condition.
#' @param control_label Value(s) in the \code{condition_key} column that mark
#'   control samples.
#' @param shrink.floor Non-negative floor for the \eqn{\tau^2} estimates.
#'   Default \code{0}.
#' @param correct.attenuation Logical; if \code{TRUE} the OLS regression slope
#'   is corrected for attenuation caused by measurement error in the shrunk
#'   batch means.  Default \code{FALSE} (the shrinkage already reduces most of
#'   the attenuation).
#' @param floor.beta Floor applied to per-sample \eqn{\beta} estimates.
#'   Default \code{0} (prevent negative betas, i.e., never correct in the wrong
#'   direction).  Set to \code{NULL} to disable.
#' @param verbose Print progress messages. Default \code{TRUE}.
#'
#' @return The Seurat object with the corrected data written back to the
#'   specified assay layer.
#'
#' @export
#' @concept batch_correction
#'
#' @seealso \code{\link{CorrectBatchMean}}
#'
#' @examples
#' \dontrun{
#' obj <- CorrectBatchMeanShrink(obj, assay = "LMC",
#'                                batch_key = "plate",
#'                                condition_key = "Condition",
#'                                control_label = "ctrl-inj")
#' }
#'
#' @importFrom SeuratObject LayerData LayerData<- DefaultAssay Cells
#' @importFrom stats median setNames
#'
CorrectBatchMeanShrink <- function(object,
                                   assay = NULL,
                                   layer = "data",
                                   batch_key = NULL,
                                   condition_key = NULL,
                                   control_label = NULL,
                                   shrink.floor = 0,
                                   correct.attenuation = FALSE,
                                   floor.beta = 0,
                                   verbose = TRUE) {

  # ---- 1. Input validation + data extraction --------------------------------
  if (!inherits(object, "Seurat")) {
    stop("'object' must be a Seurat object")
  }

  assay <- assay %||% DefaultAssay(object = object)
  if (!assay %in% names(object@assays)) {
    stop("Assay '", assay, "' not found in object")
  }

  if (is.null(batch_key)) {
    stop("Please specify 'batch_key': the metadata column identifying batches")
  }
  if (!batch_key %in% colnames(object@meta.data)) {
    stop("Column '", batch_key, "' not found in object metadata")
  }

  if (is.null(condition_key)) {
    stop("Please specify 'condition_key': the metadata column identifying conditions")
  }
  if (!condition_key %in% colnames(object@meta.data)) {
    stop("Column '", condition_key, "' not found in object metadata")
  }

  if (is.null(control_label)) {
    stop("Please specify 'control_label': the value(s) marking control samples in '",
         condition_key, "'")
  }

  if (!is.numeric(shrink.floor) || length(shrink.floor) != 1 || shrink.floor < 0) {
    stop("'shrink.floor' must be a single non-negative number")
  }
  if (!is.null(floor.beta) && (!is.numeric(floor.beta) || length(floor.beta) != 1)) {
    stop("'floor.beta' must be NULL or a single numeric value")
  }

  # --- retrieve the data matrix ---
  data_mat <- LayerData(object = object[[assay]], layer = layer)

  all_cells  <- Cells(object[[assay]])
  conditions <- object@meta.data[all_cells, condition_key]
  batches    <- object@meta.data[all_cells, batch_key]
  is_control <- conditions %in% control_label

  # global control indices
  global_ctrl_idx <- all_cells[is_control]
  if (length(global_ctrl_idx) == 0) {
    stop("No control samples found with '", condition_key, "' in ",
         paste(control_label, collapse = ", "))
  }

  G <- nrow(data_mat)
  unique_batches <- unique(batches)
  B <- length(unique_batches)

  # ---- 2. Per-batch control means Z_b (G x B matrix) -----------------------
  ctrl_idx_list <- setNames(
    lapply(unique_batches, function(b) all_cells[batches == b & is_control]),
    unique_batches
  )
  n_ctrl <- vapply(ctrl_idx_list, length, integer(1))  # named integer vector

  # Z_batch: G x B matrix (NA for batches with 0 controls)
  Z_batch <- matrix(NA_real_, nrow = G, ncol = B,
                    dimnames = list(rownames(data_mat), unique_batches))
  for (b in unique_batches) {
    nb <- n_ctrl[b]
    if (nb >= 2) {
      Z_batch[, b] <- rowMeans(data_mat[, ctrl_idx_list[[b]]])
    } else if (nb == 1) {
      Z_batch[, b] <- data_mat[, ctrl_idx_list[[b]]]
    }
    # nb == 0 stays NA
  }

  # Global control mean (across all individual control cells)
  if (length(global_ctrl_idx) >= 2) {
    Z_global <- rowMeans(data_mat[, global_ctrl_idx])
  } else {
    Z_global <- data_mat[, global_ctrl_idx]
  }

  # ---- 3. Estimate sigma^2_g (pooled within-batch variance) -----------------
  # Only batches with n_ctrl >= 2 contribute
  batches_for_sigma <- unique_batches[n_ctrl >= 2]
  df_within <- sum(n_ctrl[batches_for_sigma] - 1)

  if (df_within >= 1) {
    SS_within <- rep(0, G)
    for (b in batches_for_sigma) {
      ctrl_cells <- ctrl_idx_list[[b]]
      # centre each column by the batch control mean, then sum squares
      resid <- data_mat[, ctrl_cells, drop = FALSE] - Z_batch[, b]
      SS_within <- SS_within + rowSums(resid^2)
    }
    sigma2 <- SS_within / df_within
  } else {
    warning("df_within < 1: cannot estimate within-batch variance; ",
            "setting sigma^2 = 0 (degrades to unshrunk subtraction)")
    sigma2 <- rep(0, G)
  }

  # ---- 4. Estimate tau^2_g (between-batch variance of control means) --------
  batches_obs <- unique_batches[n_ctrl >= 1]
  B_obs <- length(batches_obs)

  if (B_obs >= 2) {
    n_vec  <- n_ctrl[batches_obs]  # n_b for batches with controls
    N_total <- sum(n_vec)
    # Weighted global mean (weighted by n_b)
    Z_weighted <- Z_batch[, batches_obs, drop = FALSE] %*% n_vec / N_total

    # MSB
    MSB <- rep(0, G)
    for (b in batches_obs) {
      MSB <- MSB + n_ctrl[b] * (Z_batch[, b] - Z_weighted)^2
    }
    MSB <- MSB / (B_obs - 1)

    # n_0 (effective sample size for unbalanced ANOVA)
    n_0 <- (N_total - sum(n_vec^2) / N_total) / (B_obs - 1)

    tau2 <- pmax((MSB - sigma2) / n_0, shrink.floor)
  } else {
    warning("Fewer than 2 batches with controls; setting tau^2 = 0 ",
            "(use global mean everywhere)")
    tau2 <- rep(0, G)
  }

  # ---- 5. Shrinkage weights and shrunk means Z*_b ---------------------------
  # lambda[g, b] = tau2_g / (tau2_g + sigma2_g / n_b)
  # Z*[g, b]     = lambda * Z_b + (1 - lambda) * Z_global
  lambda_mat <- matrix(0, nrow = G, ncol = B,
                       dimnames = list(rownames(data_mat), unique_batches))
  Z_star <- matrix(NA_real_, nrow = G, ncol = B,
                   dimnames = list(rownames(data_mat), unique_batches))

  for (b in unique_batches) {
    nb <- n_ctrl[b]
    if (nb >= 1) {
      denom <- tau2 + sigma2 / nb
      # avoid 0/0: when both tau2 and sigma2/nb are 0, set lambda = 1
      lam <- ifelse(denom > 0, tau2 / denom, 1)
      lambda_mat[, b] <- lam
      Z_star[, b] <- lam * Z_batch[, b] + (1 - lam) * Z_global
    } else {
      # no controls: lambda = 0, use global mean
      lambda_mat[, b] <- 0
      Z_star[, b] <- Z_global
    }
  }

  # ---- 6. Per-sample regression of Z*_b with attenuation correction ----------
  beta_vec <- setNames(rep(NA_real_, length(all_cells)), all_cells)

  # Track per-batch attenuation diagnostics
  atten_raw_denom <- setNames(rep(NA_real_, B), unique_batches)
  atten_me_var    <- setNames(rep(NA_real_, B), unique_batches)
  atten_denom     <- setNames(rep(NA_real_, B), unique_batches)
  atten_clamped   <- setNames(rep(FALSE, B), unique_batches)

  for (b in unique_batches) {
    batch_cells <- all_cells[batches == b]
    if (length(batch_cells) == 0) next

    z_star_b <- Z_star[, b]       # G-length vector
    z_bar    <- mean(z_star_b)
    z_cent   <- z_star_b - z_bar  # centred Z*

    raw_denom <- sum(z_cent^2)
    atten_raw_denom[b] <- raw_denom

    # Attenuation correction: subtract measurement error variance from denom
    if (isTRUE(correct.attenuation) && any(sigma2 > 0)) {
      nb <- n_ctrl[b]
      lam_b <- lambda_mat[, b]
      me_var <- sum(lam_b^2 * sigma2 / max(nb, 1))
      corrected_denom <- raw_denom - me_var
      # Clamp: at least 50% of raw denom
      clamped <- corrected_denom < 0.5 * raw_denom
      denom <- max(corrected_denom, 0.5 * raw_denom)
      atten_me_var[b]  <- me_var
      atten_denom[b]   <- denom
      atten_clamped[b] <- clamped
    } else {
      atten_me_var[b] <- 0
      atten_denom[b]  <- raw_denom
      denom <- raw_denom
    }

    if (denom < .Machine$double.eps) {
      # Z* is essentially constant across genes: beta = 0
      beta_vec[batch_cells] <- 0
      next
    }

    for (cell in batch_cells) {
      x_i    <- data_mat[, cell]
      x_bar  <- mean(x_i)
      numer  <- sum(z_cent * (x_i - x_bar))
      beta_vec[cell] <- numer / denom
    }
  }

  # Apply floor.beta if requested
  n_floored <- setNames(integer(B), unique_batches)
  if (!is.null(floor.beta)) {
    for (b in unique_batches) {
      batch_cells <- all_cells[batches == b]
      below <- beta_vec[batch_cells] < floor.beta
      n_floored[b] <- sum(below, na.rm = TRUE)
      beta_vec[batch_cells[below]] <- floor.beta
    }
    if (isTRUE(verbose)) {
      total_floored <- sum(n_floored)
      message("floor.beta = ", floor.beta, ": floored ", total_floored,
              " / ", length(all_cells), " samples")
    }
  }

  # ---- 7. Apply correction --------------------------------------------------
  for (b in unique_batches) {
    batch_cells <- all_cells[batches == b]
    if (length(batch_cells) == 0) next
    z_star_b <- Z_star[, b]
    for (cell in batch_cells) {
      data_mat[, cell] <- data_mat[, cell] - beta_vec[cell] * z_star_b
    }
  }

  # ---- 8. Verbose output and write back -------------------------------------
  if (isTRUE(verbose)) {
    # Shrinkage summary
    shrink_df <- data.frame(
      batch  = unique_batches,
      n_ctrl = n_ctrl[unique_batches],
      median_lambda = vapply(unique_batches, function(b) median(lambda_mat[, b]), numeric(1)),
      mean_lambda   = vapply(unique_batches, function(b) mean(lambda_mat[, b]),   numeric(1)),
      stringsAsFactors = FALSE
    )
    message("--- Shrinkage summary ---")
    message(paste(utils::capture.output(print(shrink_df, row.names = FALSE)), collapse = "\n"))

    # Regression summary (betas shown AFTER floor.beta, if applied)
    reg_df <- data.frame(
      batch     = unique_batches,
      n_samples = vapply(unique_batches,
                         function(b) sum(batches == b), integer(1)),
      median_beta = vapply(unique_batches,
                           function(b) median(beta_vec[all_cells[batches == b]]), numeric(1)),
      mean_beta   = vapply(unique_batches,
                           function(b) mean(beta_vec[all_cells[batches == b]]),   numeric(1)),
      sd_beta     = vapply(unique_batches,
                           function(b) {
                             vals <- beta_vec[all_cells[batches == b]]
                             if (length(vals) < 2) return(NA_real_)
                             sd(vals)
                           }, numeric(1)),
      stringsAsFactors = FALSE
    )
    if (!is.null(floor.beta)) {
      reg_df$n_floored <- n_floored[unique_batches]
    }
    message("--- Regression summary ---")
    message(paste(utils::capture.output(print(reg_df, row.names = FALSE)), collapse = "\n"))

    # Attenuation correction summary
    if (isTRUE(correct.attenuation) && any(sigma2 > 0)) {
      atten_df <- data.frame(
        batch        = unique_batches,
        total_var    = round(atten_raw_denom[unique_batches], 4),
        noise_var    = round(atten_me_var[unique_batches], 4),
        signal_var   = round(atten_denom[unique_batches], 4),
        noise_frac   = round(atten_me_var[unique_batches] /
                               pmax(atten_raw_denom[unique_batches],
                                    .Machine$double.eps), 4),
        clamped      = atten_clamped[unique_batches],
        stringsAsFactors = FALSE
      )
      message("--- Attenuation correction summary ---")
      message("  total_var  = Var_genes(Z*_b)  [raw OLS denominator]")
      message("  noise_var  = sum(lambda^2 * sigma^2 / n_b)  [measurement error in Z*]")
      message("  signal_var = total_var - noise_var  [corrected denominator, clamped at 50%]")
      message(paste(utils::capture.output(print(atten_df, row.names = FALSE)),
                    collapse = "\n"))
    }
  }

  LayerData(object[[assay]], layer = layer) <- data_mat
  if (isTRUE(verbose)) {
    message("Shrinkage + regression batch correction completed for assay '",
            assay, "'")
  }

  return(object)
}


#' Perform CellAnova Batch correction
#'
#' CellAnova (Zhaojun Zhang et al, 2024 Nat Biotech) is a new method that can remove/mitigate batch effect in single-cell data
#' and return users with a "corrected" data matrix (instead of just a low-dim embedding). Here we implement the 2nd step of the CellAnova.
#' Specifically, this function is an R-implemention of the CellAnova python package's calc_BE() function
#' (https://github.com/Janezjz/cellanova).
#'
#' This function takes a Seurat object and its pre-computed integrated embedding from
#' methods like Harmony or Seurat-CCA, a batch-effect index, and a case-control index as input, to
#' estimate the batch effect from the control samples, and correct for it from the full original
#' expression data.
#' Most of the procedures are kept the same, with the following modifications:
#' \itemize{
#'   \item currently we only support one control group.
#'   \item we have additionally implemented a future_lapply() and a more efficient regression framework to enhance the efficiency.\cr
#'   \item the procedure can be done to a "sketched" data and later project to the whole data (for the purposes of efficiency and data balance).\cr
#' }
#'
#' @param object A Seurat object
#' @param assay the name of the assay to perform the CellAnova correction.
#' @param layer the name of the layer to be used to correct for the batch effect. Should be scale.data.
#' @param integrate_key A string indicating the smallest batch unit in the meta-data (e.g., library, donor, etc.),
#' which will be used for integration later.
#' @param features features to compute corrected expression for. Defaults to the variable features set in the assay specified.
#' @param control_dict A list indicating the control-group assignment of the controls. The name of each element in the
#' list should correspond to the batch name in the 'integrate_key' column.
#' @param reduction the name of the DimReduc object we use as the integrated embeddings. Should be from methods like
#' Harmony or Seurat-integration methods (e.g., CCA).
#' @param var_cutoff the fraction of explained variance to determine the optimal value of k in truncated SVD when calculating
#' the basis of the batch effect. Default is 0.9.
#' @param k_max the maximum of singular values and vectors to compute.
#' @param k_select the user-defined number of singular values and vectors to compute (override var_cutoff and k_max).
#' Default is NULL.
#' @param new.assay.name the name for the new assay to store the corrected expression matrix
#' @param max_core the number of cores to use for parallel processing. Default is 1 (sequential).
#' @param future.memory.per.core the maximum memory (in MB) allowed per core for parallel processing. Default is 2000.
#' @param verbose display progress + messages
#' @return Returns a Seurat object with a new assay added containing the batch-corrected expression matrix.
#' When the input layer is an on-disk BPCells IterableMatrix, the corrected assay is stored as a lazy
#' BPCells transform (no data materialised to disk).
#'
#' @importFrom irlba irlba
#' @importFrom RSpectra svds
#' @importFrom Seurat DefaultAssay VariableFeatures Idents CellsByIdentities CreateAssayObject
#' @importFrom SeuratObject LayerData Embeddings CreateAssay5Object
#' @importFrom future.apply future_lapply
#' @importFrom methods slot slot<- new
#' @importMethodsFrom BPCells t %*%
#' @noRd
cellanova_calc_BE <- function(object = NULL, assay = NULL, layer = "scale.data", integrate_key = NULL,
                              features = NULL, control_dict = NULL, reduction = NULL, var_cutoff = 0.9, k_max = 1500, k_select = NULL,
                              new.assay.name = "CORRECTED", max_core = 1, future.memory.per.core = 2000,
                              verbose = TRUE){

  # Input validation
  if (is.null(object) || !inherits(object, "Seurat")) {
    stop("'object' must be a Seurat object")
  }
  if (!is.null(assay) && !assay %in% names(object@assays)) {
    stop("Assay '", assay, "' not found in object")
  }
  if (is.null(x = reduction)) {
    stop('Please provide dimensionality reduction name')
  }
  if (!reduction %in% names(object@reductions)) {
    stop("Reduction '", reduction, "' not found in object")
  }
  if (is.null(x = integrate_key)) {
    stop("Please specify the colname of the integration unit")
  }
  if (!integrate_key %in% colnames(object@meta.data)) {
    stop("Column '", integrate_key, "' not found in object metadata")
  }
  if (is.null(x = control_dict)) {
    stop("Please provide a list or a vector indicating the control samples")
  }
  if (!is.numeric(var_cutoff) || var_cutoff <= 0 || var_cutoff > 1) {
    stop("'var_cutoff' must be a numeric value between 0 and 1")
  }
  if (!is.numeric(k_max) || k_max <= 0) {
    stop("'k_max' must be a positive integer")
  }

  assay <- assay %||% DefaultAssay(object = object )
  features <- features %||% VariableFeatures(object = object[[assay]])
  if (length(x = features) == 0) {
    features <- rownames(x = LayerData(object = object[[assay]], layer = layer))
  }
  if (length(features) == 0) {
    stop("No features found. Please specify features or ensure the assay has variable features or data.")
  }

  # Detect whether the layer is an on-disk BPCells IterableMatrix
  GEX_full <- LayerData(object = object[[assay]], layer = layer)
  ondisk <- inherits(GEX_full, "IterableMatrix")

  if (ondisk && !requireNamespace("BPCells", quietly = TRUE)) {
    stop("Package 'BPCells' is required for on-disk CellAnova correction. ",
         "Please install it with: BiocManager::install('BPCells')")
  }

  # Check that the resulting matrix will not exceed R's 2^31-1 element limit (dense path only)
  n_features <- nrow(GEX_full)
  n_cells <- ncol(GEX_full)
  n_elements <- as.double(n_features) * as.double(n_cells)
  if (!ondisk && n_elements > .Machine$integer.max) {
    stop("The expression matrix (", n_features, " features x ", n_cells, " cells = ",
         format(n_elements, big.mark = ",", scientific = FALSE), " elements) exceeds R's ",
         "maximum dense matrix size (2^31 - 1 = ", format(.Machine$integer.max, big.mark = ","),
         " elements). Many internal R and Seurat operations use 32-bit indexing and will fail ",
         "on matrices this large. Consider reducing the number of cells (e.g., sketching) or ",
         "features before running CellAnova.")
  }

  if (ondisk && isTRUE(verbose)) {
    message("On-disk mode: BPCells IterableMatrix detected")
  }

  if(!is.list(control_dict)){
    control_dict <- list(g1 = control_dict)
  }
  control_groups <- names(control_dict)

  # Validate control samples exist in the data
  all_batches <- unique(object@meta.data[[integrate_key]])
  missing_batches <- setdiff(unlist(control_dict), all_batches)
  if (length(missing_batches) > 0) {
    stop("Control batches not found in data: ", paste(missing_batches, collapse = ", "))
  }

  # get the integration embeddings and scale each one
  LL <- Embeddings(object = object, reduction = reduction)

  # Remove columns that contain all zeros (more efficient vectorized approach)
  zero_cols <- which(colSums(abs(LL) > .Machine$double.eps) == 0)
  if(length(zero_cols) > 0){
    LL <- LL[, -zero_cols, drop = FALSE]
    if(isTRUE(verbose)) {
      message("Removed ", length(zero_cols), " zero columns from embeddings")
    }
  }
  # if the Embeddings are from a different assay, the cells might be different from those in 'assay'
  # need to check for this
  if( nrow(LL) != ncol(object[[assay]]) ){
    warning("Some of the cells in the reduction do not match to those in the assay.")
    LL <- LL[colnames(object[[assay]]), ]
  }

  # get the cell idendity by donor id
  Idents(object) <- integrate_key
  cell_donor_list <- CellsByIdentities(object = object, cells = colnames(object[[assay]]))

  ### regression 1: calculate the overall mean effect
  batch_list <- unique(object[[integrate_key]][colnames(object[[assay]]),1])

  # GEX_full was already extracted above (during on-disk detection)
  # For the dense path, pre-transpose once: tGEX_full is (n_cells x n_genes)
  if (!ondisk) {
    tGEX_full <- t(GEX_full)
  }

  # Set up future parallelization
  original_plan <- future::plan()
  user_has_custom_plan <- !inherits(original_plan, "sequential")
  future_modified <- FALSE

  if (!user_has_custom_plan) {
    if (max_core > 1) {
      if (isTRUE(verbose)) {
        if (.Platform$OS.type == "windows") {
          message("Setting up future multisession with ", max_core, " workers (Windows)")
        } else {
          message("Setting up future multicore with ", max_core, " workers")
        }
        message("Setting future.globals.maxSize to ",
                max_core * future.memory.per.core, " MB")
      }
      if (.Platform$OS.type == "windows") {
        future::plan(future::multisession, workers = max_core)
      } else {
        future::plan(future::multicore, workers = max_core)
      }
      options(future.globals.maxSize = max_core * future.memory.per.core * 1024^2)
      future_modified <- TRUE
    } else {
      if (isTRUE(verbose)) {
        message("Using sequential future plan (no parallelization)")
      }
      future::plan("sequential")
      future_modified <- TRUE
    }
  } else {
    if (isTRUE(verbose)) {
      message("Using user-specified future plan")
    }
  }

  if (future_modified) {
    on.exit({
      future::plan(original_plan)
      if (isTRUE(verbose)) message("Restored original future plan")
    }, add = TRUE)
  }

  if (ondisk) {
    # On-disk path: materialise each batch subset individually (no full transpose)
    # Use sequential lapply — future::multicore cannot serialise IterableMatrix pointers
    overall_bloc_coef <- lapply(batch_list, function(d) {
      cells <- cell_donor_list[[d]]
      if (length(cells) < 1) return(NULL)
      LL_tmp <- LL[cells, , drop = FALSE]
      qr_LL <- qr(LL_tmp)
      # Materialise only this batch's subset from on-disk directly in transposed form
      tGEX_mat <- as.matrix(t(GEX_full[, cells]))
      base::t(qr.coef(qr_LL, tGEX_mat))
    })
  } else {
    # Dense path: use pre-transposed full matrix with future_lapply
    overall_bloc_coef <- future_lapply(batch_list, function(d) {
      if (length(cell_donor_list[[d]]) < 1) return(NULL)
      LL_tmp <- LL[cell_donor_list[[d]], , drop = FALSE]
      qr_LL <- qr(LL_tmp)
      tGEX_mat <- tGEX_full[cell_donor_list[[d]], , drop = FALSE]
      t(qr.coef(qr_LL, tGEX_mat))
    })
    # Free the transposed matrix — no longer needed
    rm(tGEX_full)
  }
  names(overall_bloc_coef) <- batch_list
  overall_bloc_coef <- Filter(Negate(is.null), overall_bloc_coef)

  if (isTRUE(verbose)) {
    message("Done with regression for ", length(overall_bloc_coef), " batches in 'integrate_key'")
  }
  M_overall <- Reduce(`+`, overall_bloc_coef) / length(overall_bloc_coef)

  ### regression 2: calculate the batch specific effect from controls
  # an empty object to store the regression outputs
  res <- list()
  #
  for(g in control_groups){
    control_batch <- control_dict[[g]]
    counter <- 0
    bloc_coef <- list()
    # loop through each batch to compute the regression effect
    for(d in control_batch){
      # LL_tmp <- LL[cell_donor_list[[d]], ]
      # qr_LL <- qr(LL_tmp)
      # # get the coef for all the m genes
      # coef <- qr.coef(qr_LL, t(LayerData(object = object[[assay]], layer = layer, cells = cell_donor_list[[d]])))
      bloc_coef[[d]] <- overall_bloc_coef[[d]]
      counter <- counter + 1
    }

    # calculate the mean of it
    M <- Reduce(`+`, bloc_coef) / length(control_batch)
    #
    res_temp <- do.call(cbind, bloc_coef) - do.call(cbind, replicate(length(control_batch), M, simplify = FALSE))
    res[[g]] <- res_temp
  }
  #
  res_combined <- t(do.call(cbind, res))

  # Clean up intermediate objects to free memory
  rm(overall_bloc_coef, bloc_coef, res, res_temp, M)
  gc()  # Force garbage collection

  ## Perform SVD decomposition for res_combined
  # Helper: try RSpectra::svds(), fall back to base::svd() when ARPACK fails
  # (e.g. "TridiagEigen" errors on numerically rank-deficient matrices)
  safe_svds <- function(A, k, verbose = FALSE) {
    tryCatch(
      RSpectra::svds(A, k = k),
      error = function(e) {
        if (isTRUE(verbose)) {
          message("RSpectra::svds() failed (", conditionMessage(e),
                  "), falling back to base::svd()")
        }
        res <- base::svd(A, nu = 0, nv = k)
        list(d = res$d[seq_len(k)], v = res$v[, seq_len(k), drop = FALSE])
      }
    )
  }

  # Compute true total variance (Frobenius norm squared = sum of ALL singular values squared)
  # This is needed because truncated SVD only returns the top-k singular values,
  # so sum(top_k_sv^2) / sum(top_k_sv^2) always reaches 1.0 at the last component.
  total_var <- sum(res_combined^2)

  if(is.null(x = k_select)){
    k_max <- min(c(k_max, min(dim(res_combined)) - 1))

    # Two-stage SVD: run a small pilot first to estimate k
    k_pilot <- min(k_max, 100)
    k_pilot <- min(k_pilot, min(dim(res_combined)) - 1)

    if(isTRUE(verbose)) {
      message("Performing pilot SVD with k_pilot = ", k_pilot)
    }

    pilot_svds <- safe_svds(res_combined, k = k_pilot, verbose = verbose)
    DD1 <- pilot_svds$d

    # Calculate cumulative variance explained from pilot
    positive_sv <- DD1[DD1 > 0]
    if(length(positive_sv) == 0) {
      stop("No positive singular values found in SVD decomposition")
    }

    variance <- cumsum(positive_sv^2) / total_var
    k <- which(variance >= var_cutoff)[1]

    if (is.na(k) || k > k_pilot) {
      # Pilot was not enough — fall back to full k_max SVD
      k_max <- min(k_max, min(dim(res_combined)) - 1)
      if(isTRUE(verbose)) {
        message("Pilot SVD insufficient, performing full SVD with k_max = ", k_max)
      }
      final_svds <- safe_svds(res_combined, k = k_max, verbose = verbose)
      DD1 <- final_svds$d
      positive_sv <- DD1[DD1 > 0]
      variance <- cumsum(positive_sv^2) / total_var
      k <- which(variance >= var_cutoff)[1]
      if (is.na(k)) {
        k <- length(positive_sv)
        warning("Could not find k components explaining ", var_cutoff * 100, "% variance. Using all ", k, " components.")
      }
    } else {
      final_svds <- pilot_svds
    }

    if(isTRUE(verbose)) {
      message("Selected k = ", k, " components explaining ", round(variance[k] * 100, 2), "% variance")
    }
  } else {
    k <- k_select
    if(k > min(dim(res_combined)) - 1) {
      k <- min(dim(res_combined)) - 1
      warning("k_select too large, reduced to ", k)
    }

    final_svds <- safe_svds(res_combined, k = k, verbose = verbose)

    if(isTRUE(verbose)) {
      message("Using user-specified k = ", k, " components")
    }
  }
  #

  # Extract selected components
  DD1 <- final_svds$d[1:k]
  VV1T <- final_svds$v[, 1:k, drop = FALSE]

  # Store intermediate results in object tools slot for debugging/inspection
  slot(object = object, name = "tools")[["M_overall"]] <- M_overall
  slot(object = object, name = "tools")[["DD1"]] <- DD1
  slot(object = object, name = "tools")[["VV1T"]] <- VV1T

  if(isTRUE(verbose)) {
    message("Computing batch effect correction...")
  }

  if (ondisk) {
    # =========================================================================
    # On-disk path: lazy correction via BPCells TransformLinearResidual
    # =========================================================================
    # corrected = GEX_full - VV1T %*% (t(VV1T) %*% (GEX_full - M_overall %*% t(LL)))
    #
    # Rearranged as:  corrected = GEX_full - VV1T %*% correction_coef
    # where correction_coef = t(VV1T) %*% GEX_full - t(VV1T) %*% M_overall %*% t(LL)
    #
    # This is exactly X - t(row_params) %*% col_params, which BPCells'
    # TransformLinearResidual computes lazily (no full matrix materialised).

    tVV1T <- t(VV1T)  # (k x genes)

    # proj = t(VV1T) %*% GEX_full  -->  (k x cells), computed by streaming GEX_full once
    # BPCells: matrix %*% IterableMatrix returns a dense matrix
    if (isTRUE(verbose)) message("  Computing projection (streaming on-disk data)...")
    proj <- tVV1T %*% GEX_full  # (k x cells)

    # small_term = t(VV1T) %*% M_overall %*% t(LL)  -->  (k x cells) via tiny intermediates
    small_term <- (tVV1T %*% M_overall) %*% base::t(LL)  # (k x n_pcs) %*% (n_pcs x cells)

    # correction_coef = proj - small_term  -->  (k x cells)
    correction_coef <- proj - small_term
    rm(proj, small_term, tVV1T)

    # Build lazy IterableMatrix: corrected = GEX_full - t(row_params) %*% col_params
    row_params <- base::t(VV1T)   # (k x genes)
    col_params <- correction_coef  # (k x cells)

    # BPCells' TransformLinearResidual applies row_params/col_params in the
    # internal (pre-transpose) coordinate system.  When the matrix is stored
    # transposed (storage_order == "row"), we must swap the params — exactly
    # as BPCells' own regress_out() does (see BPCells/R/transforms.R).
    if (GEX_full@transpose) {
      tmp <- row_params
      row_params <- col_params
      col_params <- tmp
    }

    if(isTRUE(verbose)) {
      message("Creating corrected assay (lazy on-disk transform)...")
    }

    corrected <- new("TransformLinearResidual",
      matrix = BPCells::convert_matrix_type(GEX_full, "double"),
      transpose = GEX_full@transpose,
      dim = GEX_full@dim,
      dimnames = dimnames(GEX_full),
      row_params = row_params,
      col_params = col_params,
      global_params = numeric(0),
      vars_to_regress = character(0)
    )

    new.assay <- CreateAssay5Object(data = corrected)
    object[[new.assay.name]] <- new.assay

  } else {
    # =========================================================================
    # Dense (in-memory) path — original logic
    # =========================================================================
    # Compute residuals after removing overall mean effect
    part1 <- GEX_full - tcrossprod(M_overall, LL)

    # Project onto batch effect subspace and back to compute batch effect
    # The intermediate crossprod(VV1T, part1) is (k x cells), which is small
    be <- VV1T %*% crossprod(VV1T, part1)
    rm(part1)

    # Apply correction by subtracting estimated batch effect
    corrected <- GEX_full - be
    rm(be, GEX_full)

    if(isTRUE(verbose)) {
      message("Creating corrected assay...")
    }

    # Create new assay with corrected data
    new.assay <- suppressWarnings(
      expr = CreateAssayObject(
        data = corrected,
        min.cells = -Inf,
        min.features = -Inf,
        check.matrix = FALSE
      )
    )

    # Add corrected assay to object
    object[[new.assay.name]] <- new.assay
  }

  if(isTRUE(verbose)) {
    message("Batch correction completed. New assay '", new.assay.name, "' added to object.")
  }

  return(object)
}
