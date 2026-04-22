#' Perform preprocessing procedures for a Seurat object to prepare for the sample-level analysis
#'
#' This is a wrapper function to apply landmark sketching, training data sketching (optional), PLS learning,
#' and weighted nearest neighbor (WNN) workflow to a single-cell Seurat object to get the necessary components for
#' sample-level analyses later. The output object can be directly used by \code{\link{GenerateSampleObject}} to obtain
#' a sample-level matrix.
#'
#' @param object a Seurat object.
#' @param assay the assay to perform the analyses on
#' @param add.hvg a boolen to specify whether the function runs additional correlation test between features and response (Y) to
#'  obtain more response-relevant features to be included in VariableFeatures and used for PLS learning.
#'  See \code{\link{QuickCorTest}} for details. Default is TRUE.
#' @param group.by.CorTest A metadata column name to group cells by before QuickCorTest. If NULL,
#'   falls back to QuickCorTest without grouping. Default is NULL.
#' @param Y A metadata column name of responses. Will also be used as the responses for PLS learning.
#' @param sketch.training a boolen to specify whether the function performs sketching for the Seurat object to get a training subset for PLS learning
#'  (the PLS learned from subset will be projected to the full data). For efficiency, it is recommended to do so when the Seurat object is very
#'  large. If the 'assay' being specified for this Seurat object is an on-disk assay, this parameter must be set to TRUE. Default is FALSE.
#' @param group.by.Sketch A metadata column name to group cells by before sketching. If NULL,
#'   falls back to standard sketching without grouping. Default is NULL.
#' @param ncells.per.group A positive integer, named vector, or list specifying the number of cells to sample.
#'   Default is 1000 per category in 'group.by.Sketch'. See \code{\link{SketchDataByGroup}} for detailed usage.
#' @param training.assay.name The assay name of the training data.
#' @param training.sketch.method Sketching method to use for the training subset. Can be 'LeverageScore' or 'Uniform'.
#'   Default is 'Uniform' because for training we want to obtain an unbiased subset of the original data.
#'   See \code{\link{LeverageScore}} for details.
#' @param ncells.landmark A positive integer specifying the number of 'landmark cells'. The function will generate a new
#'   assay for these 'landmark cells' (please do not confuse it with the training subset). The landmark cells should be a
#'   subset that covers the diverse cell types and cell states of the data. Default is 2000.
#' @param landmark.assay.name The assay name of the landmark cells. Default is "LANDMARK".
#' @param landmark.sketch.method Sketching method to use for the landmark cells. Can be 'LeverageScore' or 'Uniform'.
#'   Default is 'LeverageScore' because we want the landmark cells to be as diverse/representative as possible.
#'   See \code{\link{LeverageScore}} for details.
#' @param ncomp Number of components to compute
#' @param pls.function PLS function from pls package to run (options: plsr, spls, cppls)
#' @param pls.reduction.name PLS dimensional reduction name
#'
#' @param k.nn The number of nearest neighbors to compute for each modality
#' @param weighted.nn.name Multimodal neighbor object name
#' @param fix.wnn.weights Pre-specified modality weights. If provided, skips the calculation and uses these weights directly.
#'   Should be a list with the same length as reduction.list.
#' @param name.reduction.1 The name of the DimReduc to use as the 1st embedding (the 2nd is the PLS embedding) in the WNN process.
#' @param dims.reduction.1 The dimensions for reduction.1 to use during the WNN process.
#' @param rm.training.assay Whether to remove the training assay after running PrepareSampleObject(). This is used to reduce the
#'   memory usage when the seurat object is large (e.g., a large on-disk object). Default is FALSE.
#' @param max_core The number of cores to use for parallelization when running the WNN process (but not other processes). 
#'   Note that if the user has already set the "future::plan()" before running the function, 
#'   it will ignore this parameter and respect user's future plan (and the plan()
#'   will be applied to other functions being called, including ScaleData()). 
#'   Default is 1 for sequential processing.
#' @param future.memory.per.core The memory allocation per core for options(future.globals.maxSize = ...), and the calculation 
#'   is future.globals.maxSize = max_core × future.memory.per.core × 1024 × 1024 bytes. Default is 2000 (unit in MB). 
#' @param verbose Print progress and diagnostic messages
#' @param ...	Arguments passed to other methods
#' @export
#' @concept scSLIDE
#' @return return a Seurat object that contains a weighted.nn Neighbor object between the landmark cells and all the other cells.
#'   the Seurat object will also contain a landmark assay (and a training assay if sketch.training == TRUE).
#'
#' @importFrom SeuratObject DefaultAssay DefaultAssay<- Assays VariableFeatures VariableFeatures<- CreateDimReducObject
#' @importFrom Seurat ScaleData ProjectCellEmbeddings SketchData
#' @importFrom future plan
#'
PrepareSampleObject <- function(
    object = NULL,
    assay = NULL,
    add.hvg = TRUE,
    group.by.CorTest = NULL,
    Y = NULL,
    sketch.training = TRUE,
    group.by.Sketch = NULL,
    ncells.per.group = 1000,
    training.assay.name = "TRAINING",
    training.sketch.method = "Uniform",
    ncells.landmark = 2000,
    landmark.assay.name = "LANDMARK",
    landmark.sketch.method = "LeverageScore",
    ncomp = 10,
    pls.function = c("plsr", "spls", "cppls"),
    pls.reduction.name = "pls",
    k.nn = 5,
    name.reduction.1 = "pca",
    dims.reduction.1 = 1:30,
    weighted.nn.name = "weighted.nn",
    fix.wnn.weights = c(0.5, 0.5),
    rm.training.assay = FALSE,
    max_core = 1,
    future.memory.per.core = 2000, 
    verbose = TRUE,
    ...
){
  # house-keeping checks
  assay <- assay[1L] %||% DefaultAssay(object = object)
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  name.reduction.1 <- match.arg(arg = name.reduction.1, choices = names(object@reductions))
  pls.function <- match.arg(arg = pls.function)

  # whether to perform sketching to get a subset of the full data for PLS learning
  if(isTRUE(x = sketch.training)){
    object <- SketchDataByGroup(object,
                                assay = assay,
                                ncells = ncells.per.group,
                                sketched.assay = training.assay.name,
                                method = "Uniform",
                                group.by = group.by.Sketch,
                                verbose = verbose,
                                ...)
  } else {
    training.assay.name <- assay
  }
  
  # an extra step to make sure variableFeatures is set 
  var_feature <- VariableFeatures(object)
  if(length(var_feature) <= 0) {
      stop("Please run FindVariableFeatures() or manually set VariableFeatures for this object.")
  }
  
  # whether to perform additional cor.test to get more relevant genes for PLS training
  if(isTRUE(x = add.hvg)){
    deg_list <- QuickCorTest(object,
                             assay = training.assay.name,
                             layer = "data",
                             group.by = group.by.CorTest,
                             Y = Y,
                             verbose = verbose,
                             ...)
    all_HVG <- union(VariableFeatures(object), deg_list)
    VariableFeatures(object) <- all_HVG
    #
    if(verbose){
      message("The number of HVGs used for PLS learning is ", length(all_HVG))
    }
  } else {
    all_HVG <- VariableFeatures(object)
    if(verbose){
      message("Using the default VariableFeatures for PLS learning, n = ", length(all_HVG))
    }
  }
  if(length(all_HVG) < 1) stop("The number of highly variable genes (HVGs) for running PLS is 0, something is wrong.")
  
  # to get a assay for the landmark cells
  object <- SketchData(object,
                       assay = training.assay.name,
                       ncells = ncells.landmark,
                       sketched.assay = landmark.assay.name,
                       method = landmark.sketch.method,
                       features = var_feature, 
                       verbose = verbose,
                       ...)

  # perform PLS learning for the training.assay
  DefaultAssay(object) <- training.assay.name
  object <- ScaleData(object, assay = training.assay.name, features = all_HVG)
  message("\nRunning PLS learning for the data.")
  object <- RunPLS(object = object,
                   assay = training.assay.name,
                   features = all_HVG,
                   ncomp = ncomp,
                   Y = Y,
                   pls.function = pls.function,
                   reduction.name = pls.reduction.name,
                   verbose = verbose,
                   ...)

  # to project the PLS back to the full data (if we are using a subset of training data)
  if(isTRUE(x = sketch.training)){
    proj.pls <- ProjectCellEmbeddings(
      query = object,
      query.assay = assay,
      reference = object,
      reference.assay = training.assay.name,
      reduction = pls.reduction.name,
      dims = 1:ncomp,
      scale = TRUE,
      normalization.method = "LogNormalize",
      verbose = verbose,
      ...
    )
    # append the projected PLS embeddings to the full data
    object[[paste0("proj.", pls.reduction.name)]] <- CreateDimReducObject(embeddings = proj.pls, key = "pPLS_", assay = assay)
  }

  if(rm.training.assay && sketch.training){
    if(verbose) message("removing sketch.training.assay to reduce memory usage.")
    DefaultAssay(object) <- assay
    object[[training.assay.name]] <- NULL
    gc()
  }

  # to perform multi-modal nearest neighbors between the landmarks and the other cells
  message("Constructing WNN graph between the landmark cells and all other cells.")
  DefaultAssay(object) <- assay
  # get the reduction.list ready
  if(isTRUE(x = sketch.training)){
    reduction.list <- list(name.reduction.1, paste0("proj.", pls.reduction.name))
  } else {
    reduction.list <- list(name.reduction.1, pls.reduction.name)
  }
  # get the dims.list ready
  if(dims.reduction.1[length(dims.reduction.1)] > dim(object[[name.reduction.1]])[2]){
    dims.reduction.1 <- dims.reduction.1[1]:dim(object[[name.reduction.1]])[2]
    warning("The dims.reduction.1 specified exceeds the number of components in ", name.reduction.1,
            ", so switching to ", dims.reduction.1[1], ":", dims.reduction.1[length(dims.reduction.1)])
  }
  dims.list <- list(dims.reduction.1, 1:ncomp)
  
  # ============================================================================
  # Handle future parallelization settings
  # ============================================================================
  # Check if user has already set a future plan (not sequential)
  original_plan <- future::plan()
  user_has_custom_plan <- !inherits(original_plan, "sequential")
  
  # Flag to track if we modified future settings
  future_modified <- FALSE
  
  if (!user_has_custom_plan) {
      # User has not set custom future plan
      if (max_core > 1) {
          # Set up multicore processing
          if (verbose) {
              message("Setting up future multicore with ", max_core, " workers")
              message("Setting future.globals.maxSize to ", 
                      max_core * future.memory.per.core, " MB")
          }
          
          future::plan(future::multicore, workers = max_core)
          options(future.globals.maxSize = max_core * future.memory.per.core * 1024^2)
          future_modified <- TRUE
          
      } else {
          # max_core = 1 or not specified, use sequential
          if (verbose) {
              message("Using sequential future plan (no parallelization)")
          }
          future::plan("sequential")
          future_modified <- TRUE
      }
  } else {
      # User has custom plan, respect it
      if (verbose) {
          message("Using user-specified future plan")
      }
  }
  
  # Ensure future settings are restored on function exit (if we modified them)
  if (future_modified) {
      on.exit({
          future::plan(original_plan)
          if (verbose) message("Restored original future plan")
      }, add = TRUE)
  }
  
  # run FindmmNN
  object = FindmmNN(object,
                    sketch.assay = landmark.assay.name,
                    reduction.list = reduction.list,
                    k.nn = k.nn,
                    weighted.nn.name = weighted.nn.name,
                    dims.list = dims.list,
                    fix.wnn.weights = fix.wnn.weights,
                    verbose = verbose
                   )
  # return the object
  return(object)
}


#' Generate a sample-level count matrix based on landmark assay and its weighted.nn object
#'
#' Returns summed counts of weighted.nn for each landmark cell within each sample.
#'
#' @param object Seurat object
#' @param nn.name Name of the Neighbor object to use for the calculation
#' @param return.seurat Whether to return the data as a Seurat object. Default is TRUE
#' @param k.nn the number of nearest neighbors to perform the summing
#' @param sketch.assay the name of the sketch.assay you used to perform the FindmmNN()
#' @param group.by Category (or vector of categories) for grouping (e.g, Donor ID); 'ident' by default
#' To use multiple categories, specify a vector, such as c('batch', 'replicate')
#' @param normalization.method Method for normalization. Supports LogNormalize and ChiSquared. see details at
#' \code{\link{NormalizeData}} and \code{\link{NormalizeChiSquared}}
#' @param scale.factor Scale factor for Log-Normalization, see \code{\link{NormalizeData}}
#' @param rename.group.by if rename.group.by is NULL, the rownames of the landmark matrix will used the original cell IDs.
#'  But, user can indicate which meta-data column to use to rename the rows.
#'  A suffix of "_LM" + order number will be added automatically.
#' @param add.meta.data if TRUE, the function will automatically detect sample-level meta-data (based on 'group.by') and append it
#'  to the sample-level object; if FALSE, it will not do so.
#' @param remove.sketch.cell.from.col if TRUE, the function will detect if the columns of the NN object and remove the cells that
#'  have been used as the landmark cells.
#' @param new_assay_name Name for the new assay containing landmark counts
#' @param cells.use An optional character vector of cell names to subset before building
#'   the sample-level matrix. When non-NULL, only cells in this vector (that are also present
#'   in the NN object) are retained. This is used by \code{\link{GenerateSampleObject_v2}} for
#'   cell-type proportion harmonization. Default is NULL (use all cells).
#' @param verbose Print progress and diagnostic messages
#' @param ... Arguments to be passed to methods such as \code{\link{CreateSeuratObject}}
#'
#' @export
#' @concept scSLIDE
#' @return return a Seurat object that contains a count matrix with number of landmark as rows, sample as columns
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom SeuratObject FetchData CreateSeuratObject AddMetaData
#' @importFrom Seurat NormalizeData
#'
GenerateSampleObject <- function(
    object,
    nn.name = NULL,
    k.nn = 5,
    sketch.assay = "LANDMARK",
    return.seurat = TRUE,
    new_assay_name = "LMC",
    group.by = 'ident',
    normalization.method = "ChiSquared",
    scale.factor = 10000,
    rename.group.by = NULL,
    add.meta.data = TRUE,
    remove.sketch.cell.from.col = TRUE,
    cells.use = NULL,
    verbose = TRUE,
    ...
) {
  # Get internal function from Seurat
  CreateCategoryMatrix <- getFromNamespace("CreateCategoryMatrix", "Seurat")
  
  # quick check on the group.by
  data <- FetchData(object = object, vars = rev(x = group.by))
  group.by <- intersect(group.by, colnames(data))
  if (length(group.by) < 1) stop("Please specify the correct meta-data column names.")

  # get the cells in the nn object
  if (!nn.name %in% names(object@neighbors)) stop("Please specify the correct name of the NN object.")
  slct_cells <- Cells(object[[nn.name]])

  # now extract the cell identity based on group.by
  data <- data[slct_cells, , drop = F]
  data <- data[which(rowSums(x = is.na(x = data)) == 0), , drop = F]
  for (i in 1:ncol(x = data)) {
    data[, i] <- as.factor(x = data[, i])
  }
  num.levels <- sapply(
    X = 1:ncol(x = data),
    FUN = function(i) {
      length(x = levels(x = data[, i]))
    }
  )
  if (any(num.levels == 1)) {
    message(
      paste0(
        "The following grouping variables have 1 value and will be ignored: ",
        paste0(colnames(x = data)[which(num.levels <= 1)], collapse = ", ")
      )
    )
    group.by <- rev(colnames(x = data)[which(num.levels > 1)])
    data <- data[, which(num.levels > 1), drop = F]
  }
  category.matrix <- CreateCategoryMatrix(labels = data, method = "aggregate")

  # now starts to generate cell-level count matrix
  rows <- as.vector(t(object[[nn.name]]@nn.idx[, 1:k.nn, drop = FALSE]))  # No copy if nn.idx is integer
  cols <- rep(1:length(object[[nn.name]]@cell.names), each = k.nn)
  raw_ct_mat <- sparseMatrix(i = rows, j = cols,
                             dims = c(ncol(object[[sketch.assay]]), length(object[[nn.name]]@cell.names)),
                             x = 1)
  rownames(raw_ct_mat) <- colnames(object[[sketch.assay]])
  colnames(raw_ct_mat) <- object[[nn.name]]@cell.names

  # optionally subset to user-specified cells
  if (!is.null(cells.use)) {
    keep <- which(colnames(raw_ct_mat) %in% cells.use)
    if (length(keep) == 0) stop("None of the cells in 'cells.use' found in the NN object.")
    raw_ct_mat <- raw_ct_mat[, keep]
    if (verbose) message("Subsetting to ", length(keep), " cells from 'cells.use'.")
  }

  # an extra step to remove the cells that have been used as landmarks
  if (isTRUE(x = remove.sketch.cell.from.col)){
    rm_col <- which(colnames(raw_ct_mat) %in% rownames(raw_ct_mat))
    if (length(x = rm_col) > 0){
      raw_ct_mat <- raw_ct_mat[, -rm_col]
    }
  }

  # collapse the cell-level matrix into sample-level matrix, based on the category.matrix
  landmark_ct_mat <- matrix(nrow = nrow(raw_ct_mat),
                            ncol = ncol(category.matrix))
  i <- 1
  valid_cols <- colnames(raw_ct_mat)
  for(SAMPLE_ID in colnames(category.matrix)){
    cell_idx <- rownames(category.matrix)[which(category.matrix[, SAMPLE_ID] == 1)]
    cell_idx <- cell_idx[! cell_idx %in% rownames(raw_ct_mat)]
    cell_idx <- cell_idx[cell_idx %in% valid_cols]
    if (length(cell_idx) == 0) {
      landmark_ct_mat[, i] <- 0
      i <- i + 1
      next
    }
    landmark_ct_mat[, i] <- SeuratObject::rowSums(raw_ct_mat[, cell_idx])
    i <- i + 1
  }

  # assign row and col names to the count matrix
  if (!is.null(rename.group.by)) {
    mdata <- FetchData(object = object, vars = rev(x = rename.group.by))
    rename.group.by <- intersect(rename.group.by, colnames(mdata))
    if(rename.group.by > 1) {
      rename.group.by <- rename.group.by[1]
    }
    #
    if (length(rename.group.by) < 1) {
      warning("Cannot find the correct meta-data columns to rename the landmark matrix.")
      rownames(landmark_ct_mat) <- paste0(colnames(object[[sketch.assay]]), "_LM", 1:nrow(landmark_ct_mat))
      rownames(landmark_ct_mat) <- gsub(pattern = "\\s", replacement = "_", rownames(landmark_ct_mat))
      colnames(landmark_ct_mat) <- colnames(category.matrix)
    } else {
      rownames(landmark_ct_mat) <- paste0(mdata[colnames(object[[sketch.assay]]), ], "_LM", 1:nrow(landmark_ct_mat))
      rownames(landmark_ct_mat) <- gsub(pattern = "\\s", replacement = "_", rownames(landmark_ct_mat))
      colnames(landmark_ct_mat) <- colnames(category.matrix)
    }
  } else {
    rownames(landmark_ct_mat) <- paste0(colnames(object[[sketch.assay]]), "_LM", 1:nrow(landmark_ct_mat))
    rownames(landmark_ct_mat) <- gsub(pattern = "\\s", replacement = "_", rownames(landmark_ct_mat))
    colnames(landmark_ct_mat) <- colnames(category.matrix)
  }

  if(return.seurat == FALSE) {
    return(landmark_ct_mat)
  }

  # generate a new Seurat object for the sample-level matrix
  landmark_obj <- CreateSeuratObject(counts = landmark_ct_mat, assay = new_assay_name)

  # we can now add donor-level meta data to this sample-level object and perform some basic QC and processings
  if(isTRUE(x = add.meta.data)){
    meta_data <- object@meta.data
    find_donor_specific_cols <- function(meta_data, donor_col = group.by) {
      # Check each column for within-donor consistency
      donor_specific <- sapply(names(meta_data), function(col_name) {
        if (col_name == donor_col) return(TRUE)  # Always include donor column

        # Split by donor and check if each donor has only 1 unique value
        by_donor <- split(meta_data[[col_name]], meta_data[[donor_col]])
        all(sapply(by_donor, function(x) length(unique(x)) == 1))
      })
      return(names(donor_specific)[donor_specific])
    }

    meta_data <- meta_data[, find_donor_specific_cols(meta_data, group.by)]
    meta_data <- meta_data[!duplicated(meta_data), ]
    # rownames(meta_data) <- meta_data[[group.by]]
    rownames(meta_data) <- gsub(pattern = "\\_", replacement = "-", meta_data[[group.by]])
    landmark_obj <- AddMetaData(landmark_obj, metadata = meta_data)
  }

  # normalization
  if(normalization.method == "LogNormalize"){
    landmark_obj <- NormalizeData(landmark_obj, normalization.method = normalization.method, scale.factor = scale.factor)
  } else if(normalization.method == "ChiSquared"){
    landmark_obj <- NormalizeChiSquared(object = landmark_obj, assay = new_assay_name)
  }
  #
  return(landmark_obj)
}


#' Generate a sample-level count matrix with cell-type proportion harmonization across batches
#'
#' Wraps \code{\link{GenerateSampleObject}} with a pre-processing step that
#' downsamples cells so that cell-type proportions are harmonized across batches
#' within each condition. For each condition, the target proportion for each cell
#' type is the median across batches. A bottleneck-based downsampling is then
#' applied per (condition, batch) to match those target proportions exactly.
#' Cell types whose target proportion falls below \code{min.prop} are excluded
#' from harmonization entirely and pass through unchanged; raising
#' \code{min.prop} removes minor cell types from the bottleneck calculation
#' and increases overall cell retention.
#'
#' @param object Seurat object (must already contain a Neighbor object, e.g.
#'   from \code{\link{PrepareSampleObject}})
#' @param integrate_key Metadata column identifying the batch (technical
#'   replicate / sequencing run)
#' @param condition_key Metadata column identifying the biological condition
#'   (e.g. disease vs. control)
#' @param group.by Metadata column identifying sample identity (e.g. donor).
#'   Default is \code{"ident"}.
#' @param cell_type_key Metadata column identifying the cell-type annotation
#' @param min.prop Minimum target proportion threshold. Cell types whose
#'   median (target) proportion is below this value are excluded from
#'   harmonization and all their cells are kept as-is. Set to 0 to
#'   harmonize every cell type. Default is 0.001.
#' @param seed Random seed for reproducible downsampling. Default is 42.
#' @param verbose Print progress and diagnostic messages. Default is TRUE.
#' @param nn.name Name of the Neighbor object. If NULL, inferred from the
#'   Seurat object.
#' @param k.nn Number of nearest neighbours used for aggregation. Default is 5.
#' @param sketch.assay Name of the sketch assay. Default is \code{"LANDMARK"}.
#' @param return.seurat Whether to return a Seurat object. Default is TRUE.
#' @param new_assay_name Name for the new assay. Default is \code{"LMC"}.
#' @param normalization.method Normalization method. Default is
#'   \code{"ChiSquared"}.
#' @param scale.factor Scale factor for Log-Normalization. Default is 10 000.
#' @param rename.group.by Metadata column used to rename landmark rows.
#' @param add.meta.data Attach sample-level metadata. Default is TRUE.
#' @param remove.sketch.cell.from.col Remove landmark cells from columns.
#'   Default is TRUE.
#' @param ... Additional arguments passed to
#'   \code{\link{GenerateSampleObject}}.
#'
#' @return The return value of \code{\link{GenerateSampleObject}} (a Seurat
#'   object or matrix), built from the harmonized subset of cells.
#'
#' @export
#' @concept scSLIDE
#'
#' @importFrom SeuratObject FetchData Cells
#'
GenerateSampleObject_v2 <- function(
    object,
    integrate_key,
    condition_key,
    group.by = "ident",
    cell_type_key,
    min.prop = 0.001,
    seed = 42,
    verbose = TRUE,
    # --- forwarded to GenerateSampleObject ---
    nn.name = NULL,
    k.nn = 5,
    sketch.assay = "LANDMARK",
    return.seurat = TRUE,
    new_assay_name = "LMC",
    normalization.method = "ChiSquared",
    scale.factor = 10000,
    rename.group.by = NULL,
    add.meta.data = TRUE,
    remove.sketch.cell.from.col = TRUE,
    ...
) {
  # ------------------------------------------------------------------
  # Step 1 — Validate inputs
  # ------------------------------------------------------------------
  required_keys <- c(integrate_key, condition_key, cell_type_key)
  missing_keys <- setdiff(required_keys, colnames(object@meta.data))
  if (length(missing_keys) > 0) {
    stop("The following metadata columns were not found: ",
         paste(missing_keys, collapse = ", "))
  }

  # Resolve the nn.name (same default logic as GenerateSampleObject)
  if (is.null(nn.name)) {
    nn.name <- names(object@neighbors)[1]
    if (is.null(nn.name)) stop("No Neighbor object found in the Seurat object.")
  }
  if (!nn.name %in% names(object@neighbors)) {
    stop("Neighbor object '", nn.name, "' not found in the Seurat object.")
  }

  # Get the cells that are in the NN object
  nn_cells <- Cells(object[[nn.name]])
  meta <- object@meta.data[nn_cells, , drop = FALSE]

  batch_vec     <- meta[[integrate_key]]
  condition_vec <- meta[[condition_key]]
  celltype_vec  <- meta[[cell_type_key]]

  # ------------------------------------------------------------------
  # Step 2 — Per-(condition, batch) cell-type proportions
  # ------------------------------------------------------------------
  conditions <- unique(condition_vec)
  batches    <- unique(batch_vec)
  celltypes  <- unique(celltype_vec)

  # Build a 3-way count table: condition × batch × celltype
  # Note: table() alphabetizes dimension names, so always use named indexing
  # with [cond, b, celltypes] to retrieve values in our desired order.
  count_tab <- table(condition = condition_vec,
                     batch     = batch_vec,
                     celltype  = celltype_vec)

  if (verbose) {
    message("=== Cell-type proportion harmonization ===")
    message("Conditions : ", paste(conditions, collapse = ", "))
    message("Batches    : ", paste(batches, collapse = ", "))
    message("Cell types : ", paste(celltypes, collapse = ", "))
  }

  # ------------------------------------------------------------------
  # Step 3 — Median target proportions per condition
  # ------------------------------------------------------------------
  # For each condition, compute proportion matrix (batches × celltypes),
  # then take the column-wise median.
  target_props <- list()
  for (cond in conditions) {
    # subset to batches that actually contain cells for this condition
    cond_batches <- batches[sapply(batches, function(b) sum(count_tab[cond, b, celltypes]) > 0)]
    prop_mat <- matrix(0, nrow = length(cond_batches), ncol = length(celltypes),
                       dimnames = list(cond_batches, celltypes))
    for (b in cond_batches) {
      ct_counts <- count_tab[cond, b, celltypes]
      total <- sum(ct_counts)
      if (total > 0) {
        prop_mat[b, ] <- ct_counts / total
      }
    }
    target_props[[cond]] <- apply(prop_mat, 2, stats::median)

    if (verbose) {
      message("\n--- Condition: ", cond, " ---")
      message("  Batches with cells: ", paste(cond_batches, collapse = ", "))
      for (b in cond_batches) {
        message("  Batch '", b, "' counts : ",
                paste(celltypes, "=", count_tab[cond, b, celltypes], collapse = ", "),
                " (total=", sum(count_tab[cond, b, celltypes]), ")")
        message("  Batch '", b, "' props  : ",
                paste(celltypes, "=", round(prop_mat[b, ], 4), collapse = ", "))
      }
      message("  Target proportions  : ",
              paste(celltypes, "=", round(target_props[[cond]], 4), collapse = ", "))
    }
  }

  # ------------------------------------------------------------------
  # Step 4 & 5 — Downsample and collect retained cells
  #
  # Strategy: bottleneck with min.prop filtering.
  #   - Cell types with target proportion < min.prop (or == 0) bypass
  #     harmonization entirely (all cells kept).
  #   - Among harmonized types, target proportions are renormalized to
  #     sum to 1. The bottleneck type (min of N[t] / tp_h[t]) determines
  #     the total cells to keep, and each harmonized type is set to
  #     floor(total_keep * tp_h[t]).  This guarantees exact proportion
  #     matching across batches for harmonized types.
  # ------------------------------------------------------------------
  set.seed(seed)
  retained_cells <- character(0)

  for (cond in conditions) {
    cond_batches <- batches[sapply(batches, function(b) sum(count_tab[cond, b, celltypes]) > 0)]

    # If only one batch for this condition, no harmonization needed
    if (length(cond_batches) <= 1) {
      idx <- which(condition_vec == cond & batch_vec %in% cond_batches)
      retained_cells <- c(retained_cells, nn_cells[idx])
      if (verbose) {
        message("\nCondition '", cond, "': single batch, keeping all ",
                length(idx), " cells.")
      }
      next
    }

    tp <- target_props[[cond]]

    for (b in cond_batches) {
      # Cells in this (condition, batch) group
      idx <- which(condition_vec == cond & batch_vec == b)
      cells_cb <- nn_cells[idx]
      ct_cb    <- celltype_vec[idx]

      N <- table(factor(ct_cb, levels = celltypes))

      # Types present in this batch
      types_present <- names(N)[N > 0]

      # Classify types into harmonized vs passthrough
      #   harmonized : present AND target >= min.prop
      #   passthrough: present AND target < min.prop  (includes target == 0)
      types_harmonized  <- types_present[tp[types_present] >= min.prop]
      types_passthrough <- types_present[tp[types_present] <  min.prop]

      if (length(types_harmonized) == 0) {
        # Nothing to harmonize — keep everything
        retained_cells <- c(retained_cells, cells_cb)
        if (verbose) {
          message("\nCondition '", cond, "', Batch '", b,
                  "': no harmonized types, keeping all ", length(cells_cb), " cells.")
        }
        next
      }

      # Renormalize target proportions so harmonized types sum to 1
      tp_h <- tp[types_harmonized]
      tp_h <- tp_h / sum(tp_h)

      # Bottleneck: the type with the smallest N[t] / tp_h[t] ratio
      # determines the total number of harmonized cells to keep
      ratios <- as.numeric(N[types_harmonized]) / tp_h
      total_keep <- floor(min(ratios))

      n_keep <- setNames(rep(0L, length(celltypes)), celltypes)
      for (ct in types_harmonized) {
        n_keep[ct] <- floor(total_keep * tp_h[ct])
      }
      # Passthrough types: keep all their cells
      for (ct in types_passthrough) {
        n_keep[ct] <- as.integer(N[ct])
      }

      if (verbose) {
        orig_total <- length(cells_cb)
        kept_total <- sum(n_keep)
        n_pass     <- sum(N[types_passthrough])
        # Identify the bottleneck type
        bottleneck_type <- types_harmonized[which.min(ratios)]
        message("\nCondition '", cond, "', Batch '", b,
                "': keeping ", kept_total, " / ", orig_total, " cells",
                " (", length(types_passthrough), " rare types pass through, ",
                n_pass, " cells; bottleneck: ", bottleneck_type, ")")
        message("  Per-type keep: ",
                paste(celltypes, "=", n_keep, "/", as.integer(N), collapse = ", "))
      }

      # Sample cells for each type
      for (ct in celltypes) {
        if (n_keep[ct] == 0) next
        ct_cells <- cells_cb[ct_cb == ct]
        if (length(ct_cells) <= n_keep[ct]) {
          retained_cells <- c(retained_cells, ct_cells)
        } else {
          retained_cells <- c(retained_cells, sample(ct_cells, n_keep[ct]))
        }
      }
    }
  }

  if (verbose) {
    message("\n=== Total cells retained: ", length(retained_cells),
            " / ", length(nn_cells), " ===\n")
  }

  # ------------------------------------------------------------------
  # Step 6 — Call GenerateSampleObject with cells.use
  # ------------------------------------------------------------------
  GenerateSampleObject(
    object     = object,
    nn.name    = nn.name,
    k.nn       = k.nn,
    sketch.assay = sketch.assay,
    return.seurat = return.seurat,
    new_assay_name = new_assay_name,
    group.by   = group.by,
    normalization.method = normalization.method,
    scale.factor = scale.factor,
    rename.group.by = rename.group.by,
    add.meta.data = add.meta.data,
    remove.sketch.cell.from.col = remove.sketch.cell.from.col,
    cells.use  = retained_cells,
    verbose    = verbose,
    ...
  )
}
