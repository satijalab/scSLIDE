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
  # Compute true total variance (Frobenius norm squared = sum of ALL singular values squared)
  # This is needed because truncated SVD only returns the top-k singular values,
  # so sum(top_k_sv^2) / sum(top_k_sv^2) always reaches 1.0 at the last component.
  total_var <- sum(res_combined^2)

  if(is.null(x = k_select)){
    k_max <- min(c(k_max, dim(res_combined) - 1))

    # Two-stage SVD: run a small pilot first to estimate k
    k_pilot <- min(k_max, 100)

    if(isTRUE(verbose)) {
      message("Performing pilot SVD with k_pilot = ", k_pilot)
    }

    pilot_svds <- RSpectra::svds(res_combined, k = k_pilot)
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
      if(isTRUE(verbose)) {
        message("Pilot SVD insufficient, performing full SVD with k_max = ", k_max)
      }
      final_svds <- RSpectra::svds(res_combined, k = k_max)
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

    final_svds <- RSpectra::svds(res_combined, k = k)

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
