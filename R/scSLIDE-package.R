#' scSLIDE: Single Cell Sketching and Landmark Integrated Dimensional Embedding
#'
#' scSLIDE extends Seurat with advanced functionality for single-cell RNA sequencing analysis,
#' including partial least squares (PLS) dimensionality reduction, enhanced sketching methods,
#' trajectory analysis, and sample-level aggregation methods for large-scale single-cell datasets.
#'
#' @section Key Functions:
#' \describe{
#'   \item{RunPLS}{Perform Partial Least Squares dimensionality reduction}
#'   \item{SketchDataByGroup}{Enhanced sketching with group-aware sampling}
#'   \item{FindmmNN}{Find multi-modal nearest neighbors}
#'   \item{PrepareSampleObject}{Prepare data for sample-level analysis}
#'   \item{GenerateSampleObject}{Generate sample-level count matrices}
#'   \item{TrajDETest}{Trajectory-based differential expression analysis}
#'   \item{QuickCorTest}{Fast correlation testing between genes and responses}
#'   \item{NormalizeChiSquared}{Chi-squared normalization for sample-level data}
#'   \item{RunDiffusionMap}{Diffusion map dimensionality reduction}
#' }
#'
#' @section Visualization Functions:
#' \describe{
#'   \item{BuildLandmarkObject}{Build landmark objects for visualization}
#'   \item{PlotLandmarkObject}{Plot landmark-trajectory correlations}
#'   \item{SampleLevelDimPlot}{Sample-level visualization}
#'   \item{RunAndProjectUMAP}{Generate and project UMAP embeddings}
#' }
#'
#' @useDynLib scSLIDE, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"

# Operator for setting default values
`%||%` <- function(lhs, rhs) {
  if (!is.null(lhs)) {
    lhs
  } else {
    rhs
  }
}