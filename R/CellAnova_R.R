#' @include generics.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
#' @param verbose display progress + messages
#' @return Returns a Seurat object with a new assay added containing the batch-corrected expression matrix
#'
#' @importFrom irlba irlba
#'

cellanova_calc_BE <- function(object = NULL, assay = NULL, layer = "scale.data", integrate_key = NULL,
                              features = NULL, control_dict = NULL, reduction = NULL, var_cutoff = 0.9, k_max = 1500, k_select = NULL,
                              new.assay.name = "CORRECTED", verbose = TRUE){

  assay <- assay %||% DefaultAssay(object = object )
  if (is.null(x = reduction)) {
    stop('Please provide dimensionality reduction name')
  }
  if (is.null(x = integrate_key)) {
    stop("Please specify the colname of the integration unit")
  }
  if (is.null(x = control_dict)) {
    stop("Please provide a list or a vector indicating the control samples")
  }
  features <- features %||% VariableFeatures(object = object[[assay]])
  if (length(x = features) == 0) {
    features <- rownames(x = LayerData(object = object[[assay]], layer = layer))
  }
  if(is.vector(control_dict)){
    control_dict <- list(g1 = control_dict)
  }
  control_groups <- names(control_dict)

  # get the integration embeddings and scale each one
  LL <- Embeddings(object = object, reduction = reduction)

  # need to make sure the columns that contain all zeros are removed from LL
  rm_idx <- vector()
  for(i in 1:ncol(LL)){
    if(all(abs(LL[, i]) <= .Machine$double.eps)){
      rm_idx <- c(rm_idx, i)
    }
  }
  if(length(rm_idx) > 0){
    LL <- LL[, -rm_idx]
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
  overall_bloc_coef <- list()
  counter <- 0
  for(d in batch_list){
    if(length(cell_donor_list[[d]]) < 1) next()

    LL_tmp <- LL[cell_donor_list[[d]], ]
    rownames(LL_tmp) <- cell_donor_list[[d]]
    #
    qr_LL <- qr(LL_tmp)
    # get the coef for all the m genes
    GEX_mat <- LayerData(object = object[[assay]], layer = layer, cells = cell_donor_list[[d]])

    coef <- qr.coef(qr_LL, t(GEX_mat))
    overall_bloc_coef[[d]] <- t(coef)
    counter <- counter + 1
    if(isTRUE(verbose)) {
      message("Done with processing ", d, " in 'integrate_key'")
    }
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

  # clean up a bit
  rm(overall_bloc_coef, bloc_coef, res, res_temp)

  ## next perform SVD decomposition for res_combined
  if(is.null(x = k_select)){
    k_max <- min(c(k_max, dim(res_combined)-1))
    #
    # initial_svds <- irlba::irlba(res_combined, nv = k_max)
    initial_svds <- RSpectra::svds(res_combined, k = k_max)
    DD1 <- initial_svds$d
    # Sort positive singular values in descending order
    DD1_rev <- sort(sqrt(DD1[DD1 > 0]), decreasing = TRUE)
    variance <- cumsum(DD1_rev^2) / sum(DD1_rev^2)
    k <- which(variance >= var_cutoff)[1]
    #
    if(is.na(x = k)) k <- length(DD1_rev)
    #
    final_svds <- initial_svds
  } else {
    k <- k_select
    # final_svds <- irlba::irlba(res_combined, nv = k)
    final_svds <- RSpectra::svds(res_combined, k = k)
  }
  #

  DD1 <- final_svds$d[1:k]
  VV1T <- final_svds$v[, 1:k]

  slot(object = object, name = "tools")[["M_overall"]] <- M_overall
  slot(object = object, name = "tools")[["DD1"]] <- DD1
  slot(object = object, name = "tools")[["VV1T"]] <- VV1T

  # Calculate batch effect
  # Using %*% for matrix multiplication (equivalent to @ in Python)
  part1 <- LayerData(object = object[[assay]], layer = layer) - tcrossprod(M_overall, LL)
  # part2 <- tcrossprod(VV1T)
  be <- VV1T %*% crossprod(VV1T, part1)

  # Apply correction
  corrected <- LayerData(object = object[[assay]], layer = layer) - be
  #
  new.assay <- suppressWarnings(
    expr = CreateAssayObject(
      data =  corrected,
      min.cells = -Inf,
      min.features = -Inf,
      check.matrix = FALSE
    )
  )
  object[[new.assay.name]] <- new.assay
  return(object)
}






