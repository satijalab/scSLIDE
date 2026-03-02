# Test script for cellanova_calc_BE: dense path vs on-disk BPCells path
# Run with: Rscript tests/test_cellanova.R

library(Seurat)
library(SeuratObject)
library(BPCells)
devtools::load_all(".")

set.seed(42)

# ============================================================================
# 1. Create a small test Seurat object
# ============================================================================
cat("=== Creating test Seurat object ===\n")

n_genes <- 200
n_cells <- 500
n_batches <- 5
cells_per_batch <- n_cells / n_batches

# Simulate count matrix
counts <- matrix(rpois(n_genes * n_cells, lambda = 5),
                 nrow = n_genes, ncol = n_cells)
rownames(counts) <- paste0("Gene", seq_len(n_genes))
colnames(counts) <- paste0("Cell", seq_len(n_cells))

# Create Seurat object
obj <- CreateSeuratObject(counts = counts)

# Add batch metadata
batch_ids <- rep(paste0("Batch", seq_len(n_batches)), each = cells_per_batch)
obj$batch <- batch_ids

# Standard preprocessing
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = n_genes, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 20, verbose = FALSE)

# Define control batches (use first 3 batches as controls)
control_dict <- list(g1 = paste0("Batch", 1:3))

cat("Object created:", n_genes, "genes x", n_cells, "cells,", n_batches, "batches\n")
cat("Control batches:", paste(control_dict$g1, collapse = ", "), "\n")

# ============================================================================
# 2. Run dense path
# ============================================================================
cat("\n=== Running dense path ===\n")

obj_dense <- cellanova_calc_BE(
  object = obj,
  assay = "RNA",
  layer = "scale.data",
  integrate_key = "batch",
  control_dict = control_dict,
  reduction = "pca",
  var_cutoff = 0.9,
  k_max = 50,
  new.assay.name = "CORRECTED_DENSE",
  verbose = TRUE
)

cat("Dense path completed.\n")
cat("  M_overall dims:", paste(dim(obj_dense@tools$M_overall), collapse = " x "), "\n")
cat("  DD1 length:", length(obj_dense@tools$DD1), "\n")
cat("  VV1T dims:", paste(dim(obj_dense@tools$VV1T), collapse = " x "), "\n")
cat("  Corrected assay class:", class(obj_dense[["CORRECTED_DENSE"]]), "\n")

corrected_dense <- LayerData(obj_dense[["CORRECTED_DENSE"]], layer = "data")
cat("  Corrected data class:", class(corrected_dense)[1], "\n")
cat("  Corrected data dims:", paste(dim(corrected_dense), collapse = " x "), "\n")
cat("  Corrected data range:", round(range(corrected_dense), 4), "\n")

# ============================================================================
# 3. Create on-disk version of the same object
# ============================================================================
cat("\n=== Setting up on-disk object ===\n")

ondisk_dir <- file.path(tempdir(), "cellanova_test_ondisk")
if (dir.exists(ondisk_dir)) unlink(ondisk_dir, recursive = TRUE)
dir.create(ondisk_dir, recursive = TRUE)

# Write counts to BPCells on-disk format
counts_dir <- file.path(ondisk_dir, "counts")
bpcells_counts <- BPCells::write_matrix_dir(
  mat = as(counts, "dgCMatrix"),
  dir = counts_dir
)

# Create a new Seurat v5 object with on-disk counts
obj_ondisk <- CreateSeuratObject(counts = bpcells_counts)
obj_ondisk$batch <- batch_ids

# Run the same preprocessing — ScaleData on IterableMatrix produces a lazy chain
obj_ondisk <- NormalizeData(obj_ondisk, verbose = FALSE)
obj_ondisk <- FindVariableFeatures(obj_ondisk, nfeatures = n_genes, verbose = FALSE)
obj_ondisk <- ScaleData(obj_ondisk, verbose = FALSE)
obj_ondisk <- RunPCA(obj_ondisk, npcs = 20, verbose = FALSE)

# Check the layer type
scale_data <- LayerData(obj_ondisk[["RNA"]], layer = "scale.data")
cat("scale.data class:", class(scale_data)[1], "\n")
cat("Is IterableMatrix:", inherits(scale_data, "IterableMatrix"), "\n")

# ============================================================================
# 4. Run on-disk path
# ============================================================================
cat("\n=== Running on-disk path ===\n")

obj_ondisk <- cellanova_calc_BE(
  object = obj_ondisk,
  assay = "RNA",
  layer = "scale.data",
  integrate_key = "batch",
  control_dict = control_dict,
  reduction = "pca",
  var_cutoff = 0.9,
  k_max = 50,
  new.assay.name = "CORRECTED_ONDISK",
  verbose = TRUE
)

cat("On-disk path completed.\n")
cat("  M_overall dims:", paste(dim(obj_ondisk@tools$M_overall), collapse = " x "), "\n")
cat("  DD1 length:", length(obj_ondisk@tools$DD1), "\n")
cat("  VV1T dims:", paste(dim(obj_ondisk@tools$VV1T), collapse = " x "), "\n")
cat("  Corrected assay class:", class(obj_ondisk[["CORRECTED_ONDISK"]]), "\n")

corrected_ondisk <- LayerData(obj_ondisk[["CORRECTED_ONDISK"]], layer = "data")
cat("  Corrected data class:", class(corrected_ondisk)[1], "\n")
cat("  Corrected data dims:", paste(dim(corrected_ondisk), collapse = " x "), "\n")

# Materialise for comparison
corrected_ondisk_mat <- as.matrix(corrected_ondisk)
cat("  Corrected data range:", round(range(corrected_ondisk_mat), 4), "\n")

# ============================================================================
# 5. Compare results
# ============================================================================
cat("\n=== Comparing dense vs on-disk results ===\n")

# Ensure same cell and gene ordering
common_genes <- intersect(rownames(corrected_dense), rownames(corrected_ondisk_mat))
common_cells <- intersect(colnames(corrected_dense), colnames(corrected_ondisk_mat))
cat("Common genes:", length(common_genes), " Common cells:", length(common_cells), "\n")

dense_sub <- as.matrix(corrected_dense[common_genes, common_cells])
ondisk_sub <- corrected_ondisk_mat[common_genes, common_cells]

# M_overall comparison
M_dense <- obj_dense@tools$M_overall
M_ondisk <- obj_ondisk@tools$M_overall
m_diff <- max(abs(M_dense - M_ondisk))
cat("M_overall max abs diff:", format(m_diff, scientific = TRUE), "\n")

# VV1T comparison (sign may flip — compare absolute values or use abs correlation)
VV1T_dense <- obj_dense@tools$VV1T
VV1T_ondisk <- obj_ondisk@tools$VV1T
cat("VV1T dims dense:", paste(dim(VV1T_dense), collapse = "x"),
    " ondisk:", paste(dim(VV1T_ondisk), collapse = "x"), "\n")

# Corrected data comparison
max_diff <- max(abs(dense_sub - ondisk_sub))
mean_diff <- mean(abs(dense_sub - ondisk_sub))
cor_val <- cor(as.vector(dense_sub), as.vector(ondisk_sub))

cat("Corrected data max abs diff:", format(max_diff, scientific = TRUE), "\n")
cat("Corrected data mean abs diff:", format(mean_diff, scientific = TRUE), "\n")
cat("Corrected data correlation:", format(cor_val, digits = 10), "\n")

# ============================================================================
# 6. Verdict
# ============================================================================
cat("\n=== Test verdict ===\n")
tol <- 1e-8
if (max_diff < tol) {
  cat("PASS: Dense and on-disk results are numerically identical (max diff < ", tol, ")\n")
} else if (cor_val > 0.9999) {
  cat("PASS (with tolerance): Results are highly correlated (r =", format(cor_val, digits = 8),
      ") but max diff =", format(max_diff, scientific = TRUE), "\n")
} else {
  fail_msg <- paste0(
    "FAIL: Results differ substantially. max diff = ",
    format(max_diff, scientific = TRUE),
    ", correlation = ",
    format(cor_val, digits = 8)
  )
  cat(fail_msg, "\n")
  stop(fail_msg)
}

# Clean up
unlink(ondisk_dir, recursive = TRUE)
cat("\nDone. Temp directory cleaned up.\n")
