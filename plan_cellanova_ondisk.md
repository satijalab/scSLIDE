# Plan: On-Disk `cellanova_calc_BE()` for BPCells-Backed Seurat Objects

## Context

The current `cellanova_calc_BE()` operates entirely on dense in-memory matrices.
With ~1M cells and 3,000 features the full expression matrix alone is ~24 GB, and
the correction step creates several such matrices simultaneously. An on-disk
implementation would remove the cell-count ceiling by never materialising the full
(genes × cells) matrix into RAM.

This plan is informed by the following BPCells capabilities (verified from source):

| Operation | BPCells support | Notes |
|---|---|---|
| `LayerData(layer = "scale.data")` | Returns `IterableMatrix` when the assay is on-disk | Lazy transformation chain from `ScaleData.IterableMatrix` |
| `t()` | Lazy | Flips the `@transpose` flag |
| `mat[i, j]` | Lazy | Row/column subsetting |
| `as.matrix(mat[i, j])` | Materialises subset to dense | Needed for `qr.coef` |
| `dense %*% IterableMatrix` | Returns dense | Streams through on-disk columns |
| `IterableMatrix %*% dense` | Returns dense | Same |
| `TransformLinearResidual` | Lazy IterableMatrix | Represents `X − R^T C` without materialising; also overrides `%*%` to exploit low-rank structure |
| `BPCells::svds()` | Iterative SVD on IterableMatrix | Same return format as `RSpectra::svds` |
| `write_matrix_dir()` | Writes IterableMatrix to on-disk directory | Returns `MatrixDir` (IterableMatrix) |
| `cbind(IterableMatrix, IterableMatrix)` | Lazy `ColBindMatrices` | No data copied; stays on-disk |
| `crossprod()` / `tcrossprod()` | **NOT supported** | Must be rewritten via `t()` + `%*%` |
| `qr()` / `qr.coef()` | **NOT supported on IterableMatrix** | Must materialise the batch subset to dense first |
| `IterableMatrix − IterableMatrix` | **NOT supported** | Only `IterableMatrix − numeric` |
| `CreateAssay5Object(data = IterableMatrix)` | Natively supported | Stored without conversion |

---

## Architecture Overview

```
cellanova_calc_BE()
  │
  ├── detect whether scale.data layer is IterableMatrix
  │     ├── YES → on-disk path (new)
  │     └── NO  → existing dense path (unchanged)
  │
  └── on-disk path
        ├─ Stage 1  per-batch regression (materialise per-batch subset)
        ├─ Stage 2  residual matrix (small dense, unchanged)
        ├─ Stage 3  SVD via BPCells::svds() or RSpectra
        ├─ Stage 4  correction via TransformLinearResidual + chunked write
        └─ Stage 5  store result as on-disk assay layer
```

A single function with an internal branch, rather than a separate S3 method,
keeps the public API identical (`cellanova_calc_BE(object, ...)`). The user does
not need to know which path runs — it is selected automatically.

---

## Step 0 — Prerequisites

### 0a. Guard and detection

After the existing input-validation block (and after features / assay are resolved),
add:

```r
GEX_full <- LayerData(object = object[[assay]], layer = layer)
ondisk <- inherits(GEX_full, "IterableMatrix")

if (ondisk && !requireNamespace("BPCells", quietly = TRUE)) {
  stop("Package 'BPCells' is required for on-disk CellAnova correction.")
}
```

The rest of the function branches on `ondisk`.

### 0b. New parameters

| Parameter | Default | Purpose |
|---|---|---|
| `ondisk.path` | `NULL` | Directory for on-disk output. Required when input is on-disk. If `NULL` and on-disk is detected, a `tempdir()` sub-directory is created with a warning. |
| `chunk_size` | `NULL` | Number of cells per chunk in the correction step. Default: `floor(.Machine$integer.max / n_features)` (ensures each chunk stays under 2^31 elements). |

These parameters are ignored on the dense path.

### 0c. Remove the 2^31 early-stop for on-disk path

The recently-added size check (`n_elements > .Machine$integer.max → stop(...)`)
should only fire on the dense path. On the on-disk path, chunked processing
handles arbitrarily large matrices.

---

## Step 1 — Per-Batch Regression (on-disk path)

### What changes

The regression needs `qr.coef(qr_LL, tGEX_mat)` where `tGEX_mat` is
(cells_in_batch × genes). BPCells does not support `qr.coef` on
IterableMatrix. However, each batch is typically thousands to tens-of-thousands
of cells — well under 2^31 — so materialising one batch at a time is safe.

### Implementation

```r
# On the dense path, t(GEX_full) is pre-computed in one shot.
# On the on-disk path, we skip the full transpose and materialise per-batch.

overall_bloc_coef <- future_lapply(batch_list, function(d) {
  cells <- cell_donor_list[[d]]
  if (length(cells) < 1) return(NULL)

  LL_tmp <- LL[cells, , drop = FALSE]
  qr_LL  <- qr(LL_tmp)

  if (ondisk) {
    # Materialise only this batch's subset (cells x genes)
    tGEX_mat <- as.matrix(t(GEX_full[, cells, drop = FALSE]))
  } else {
    tGEX_mat <- tGEX_full[cells, , drop = FALSE]
  }

  t(qr.coef(qr_LL, tGEX_mat))
})
```

Memory cost per batch: O(cells_in_batch × genes). With a typical batch of
10K–50K cells and 3K genes this is 30M–150M elements (~240 MB–1.2 GB), entirely
manageable.

**Note on parallelism**: When `ondisk = TRUE` and `max_core > 1`, each worker
would independently stream its batch from disk. BPCells on-disk reads are
thread-safe (separate file handles), so `future::multicore` works. However,
`future::multisession` would require serialising the IterableMatrix pointer,
which is not supported. We should document that `future::multicore` (or
`future::sequential`) must be used, or fall back to sequential if the platform
does not support forking (Windows).

---

## Step 2 — Residual Computation

No change. The `overall_bloc_coef` elements are small dense matrices
(genes × embedding_dims). All operations (`Reduce`, `cbind`, subtraction) work
on these small matrices regardless of the input path.

`res_combined` is (n_control_batches × dims) × genes — always small.

---

## Step 3 — SVD

### What changes

`res_combined` is a small dense matrix. Both `RSpectra::svds` and the two-stage
pilot approach work identically. No change required for the on-disk path.

(If in the future `res_combined` becomes large enough to warrant on-disk SVD,
`BPCells::svds(as(res_combined, "IterableMatrix"), k)` could be used, but this
is not expected.)

---

## Step 4 — Batch-Effect Correction (the critical stage)

### The algebra

```
part1     = GEX_full − M_overall %*% t(LL)          # (genes × cells)
be        = VV1T %*% crossprod(VV1T, part1)          # (genes × cells)
corrected = GEX_full − be                            # (genes × cells)
```

### Key insight: `TransformLinearResidual`

BPCells provides `TransformLinearResidual`, a lazy IterableMatrix that computes
`X − row_params^T %*% col_params` on-the-fly without materialising the full
result. It also overrides `%*%` to exploit the low-rank structure:

```
(X − R^T C) %*% B  =  X %*% B  −  R^T (C %*% B)
```

This means `t(VV1T) %*% part1` can be computed efficiently: BPCells will stream
through the sparse on-disk matrix once, and the low-rank correction is applied
via small dense multiplies.

### Implementation

```r
# ---- part1 as a lazy IterableMatrix ----
# part1 = GEX_full − M_overall %*% t(LL)
#
# TransformLinearResidual represents: X − row_params^T %*% col_params
#   row_params = t(M_overall)   →  (dims × genes)
#   col_params = t(LL)          →  (dims × cells)

row_params <- t(M_overall)   # (dims × genes)
col_params <- t(LL)          # (dims × cells)

# Handle BPCells storage order
if (BPCells::storage_order(GEX_full) == "row") {
  tmp <- row_params
  row_params <- col_params
  col_params <- tmp
}

part1 <- methods::new(
  "TransformLinearResidual",
  matrix       = BPCells::convert_matrix_type(GEX_full, "double"),
  transpose    = GEX_full@transpose,
  dim          = GEX_full@dim,
  dimnames     = dimnames(GEX_full),
  row_params   = as.matrix(row_params),
  col_params   = as.matrix(col_params),
  global_params    = numeric(0),
  vars_to_regress  = character(0)
)
# part1 is now a lazy IterableMatrix — O(dims × genes + dims × cells) memory

# ---- Compute VV1T^T %*% part1 using BPCells %*% ----
# t(VV1T) is (k × genes);  part1 is (genes × cells) lazy IterableMatrix
# Result is (k × cells) dense — k is typically 10–50, so this is small.
proj <- t(VV1T) %*% part1   # dispatches to dense_multiply_left on IterableMatrix
                              # exploits low-rank override: O(stream) + O(k × dims × cells)

# ---- be = VV1T %*% proj ----
# VV1T is (genes × k);  proj is (k × cells)
# Result is (genes × cells) dense — THIS is the one full-size dense matrix.
be <- VV1T %*% proj          # plain dense multiply, (genes × cells)
rm(proj)
```

At this point `be` is a dense (genes × cells) matrix. This is unavoidable:
the batch effect `be` is a dense low-rank matrix with no sparsity to exploit.

### Avoiding a second full-size dense matrix: chunked `corrected = GEX_full − be`

`IterableMatrix − matrix` is not supported. We cannot represent `corrected`
as a single lazy object. Instead, compute and write in chunks:

```r
all_cells   <- colnames(GEX_full)
n_cells     <- length(all_cells)
chunk_size  <- chunk_size %||% floor(.Machine$integer.max / n_features)
chunk_starts <- seq(1, n_cells, by = chunk_size)

chunk_dirs <- file.path(ondisk.path, paste0("chunk_", seq_along(chunk_starts)))
chunk_mats <- vector("list", length(chunk_starts))

for (i in seq_along(chunk_starts)) {
  idx <- chunk_starts[i]:min(chunk_starts[i] + chunk_size - 1, n_cells)
  cells_i <- all_cells[idx]

  # Materialise this chunk of GEX_full
  GEX_chunk <- as.matrix(GEX_full[, cells_i, drop = FALSE])

  # Subtract the corresponding columns of be
  corrected_chunk <- GEX_chunk - be[, idx, drop = FALSE]
  rm(GEX_chunk)

  # Write chunk to disk as BPCells on-disk matrix
  corrected_chunk <- as(corrected_chunk, "dgCMatrix")  # or keep dense
  chunk_mats[[i]] <- BPCells::write_matrix_dir(
    mat       = corrected_chunk,
    dir       = chunk_dirs[i],
    compress  = FALSE
  )
  rm(corrected_chunk)
}

rm(be)

# Lazy cbind — still fully on-disk
corrected <- Reduce(cbind, chunk_mats)

# Optionally consolidate into a single directory
corrected <- BPCells::write_matrix_dir(corrected, file.path(ondisk.path, "corrected"),
                                        compress = FALSE, overwrite = TRUE)
# Clean up chunk directories
lapply(chunk_dirs, unlink, recursive = TRUE)
```

**Peak memory**: one chunk of GEX_full (dense, ≤ 2^31 elements) + the full `be`
matrix. With 3K genes × 1M cells and chunk_size = 715K:

- `be`: 3K × 1M × 8 bytes = ~24 GB (unavoidable — it is the batch effect)
- chunk: 3K × 715K × 8 bytes = ~17 GB
- Peak: ~41 GB (vs. ~96 GB for the original 4-matrix dense approach)

### Alternative: avoid materialising `be` entirely (fully streaming)

If even 24 GB for `be` is too much, the chunk loop can compute `be` per-chunk
on the fly from `VV1T` and `proj`:

```r
# Instead of pre-computing the full be:
for (i in seq_along(chunk_starts)) {
  idx <- chunk_starts[i]:min(chunk_starts[i] + chunk_size - 1, n_cells)
  cells_i <- all_cells[idx]

  GEX_chunk <- as.matrix(GEX_full[, cells_i, drop = FALSE])

  # Compute part1 for this chunk
  part1_chunk <- GEX_chunk - M_overall %*% t(LL[cells_i, , drop = FALSE])

  # Project + correct
  be_chunk <- VV1T %*% (t(VV1T) %*% part1_chunk)
  corrected_chunk <- GEX_chunk - be_chunk

  rm(GEX_chunk, part1_chunk, be_chunk)
  # ... write_matrix_dir as above ...
}
```

This recomputes `part1` and `be` per-chunk but avoids ever holding a full
(genes × cells) dense matrix. **Peak memory**: 2 × chunk_size × genes × 8 bytes
(one for `GEX_chunk`, one for `part1_chunk` or `be_chunk` or `corrected_chunk`).
With 3K genes and 715K chunk size: ~34 GB peak. With a smaller chunk (e.g. 200K):
~10 GB.

**Trade-off**: the fully-streaming approach reads `GEX_full` from disk twice per
chunk (once for `GEX_chunk`, and the `part1` lazy IterableMatrix already streams
it once for `t(VV1T) %*% part1`). However, since we're materialising `GEX_chunk`
directly and computing `part1_chunk` from it in-memory, there is actually only
one disk read per chunk. The cost is recomputing `M_overall %*% t(LL[chunk])` per
chunk — a cheap (genes × dims) × (dims × chunk_size) multiply.

**Recommendation**: Use the fully-streaming approach (no pre-materialised `be`).
It is simpler, has strictly lower peak memory, and the per-chunk overhead is
negligible.

---

## Step 5 — Store Result

### On-disk path

```r
# corrected is an on-disk IterableMatrix (MatrixDir or ColBindMatrices)
# Store it directly in a Seurat v5 assay layer
LayerData(object[[new.assay.name]], layer = "data") <- corrected
```

Or, using `CreateAssay5Object`:

```r
new.assay <- CreateAssay5Object(data = corrected)
object[[new.assay.name]] <- new.assay
```

`CreateAssay5Object` natively accepts `IterableMatrix`. The corrected data
remains on-disk.

### Dense path

Unchanged — uses `CreateAssayObject(data = corrected)` as today.

---

## Step 6 — Update Package Infrastructure

### DESCRIPTION

No new hard dependency. BPCells is already in `Suggests`. Add a note if desired:

```
Suggests:
    BPCells (>= 0.2.0),
    ...
```

### Roxygen documentation

Add to `cellanova_calc_BE`:

```r
#' @param ondisk.path Character. Directory path for on-disk output when the
#'   input assay uses a BPCells on-disk layer. Each correction chunk is written
#'   to a subdirectory. If NULL (default) and the input is on-disk, a
#'   temporary directory is used with a warning.
#' @param chunk_size Integer. Number of cells per chunk in the on-disk
#'   correction step. Default: floor(.Machine$integer.max / n_features).
```

---

## Summary of Memory Profile

| Scenario | Peak RAM | Notes |
|---|---|---|
| **Current dense** (3K × 1M) | ~48 GB | After Step 3 optimisation (2 full matrices) |
| **On-disk, pre-computed `be`** | ~41 GB | 24 GB for `be` + 17 GB for one chunk |
| **On-disk, fully streaming** (default chunk 715K) | ~34 GB | 2 chunks in flight |
| **On-disk, fully streaming** (chunk 200K) | ~10 GB | User-tunable via `chunk_size` |
| **On-disk, fully streaming** (chunk 50K) | ~2.4 GB | Minimal memory, more disk I/O passes |

The fully-streaming approach with a user-configurable `chunk_size` gives users
direct control over the memory-vs-speed trade-off.

---

## Risks and Mitigations

| Risk | Mitigation |
|---|---|
| `TransformLinearResidual` is not exported by BPCells — constructing it via `methods::new()` uses internal API | Guard with `requireNamespace("BPCells")`, pin minimum BPCells version in Suggests, add a unit test that catches breakage. Also: the fully-streaming approach (recommended) does **not** need `TransformLinearResidual` at all — it materialises `part1_chunk` in-memory per chunk. |
| `future::multisession` cannot serialise IterableMatrix pointers | Document that `future::multicore` is required; detect and warn on Windows (no fork support → fall back to sequential). |
| Column-name ordering after chunked cbind | Preserve original `colnames(GEX_full)` ordering when constructing chunks; BPCells `cbind` concatenates names in order. |
| Temporary chunk directories left behind on error | Wrap the chunk loop in `on.exit(unlink(chunk_dirs, recursive = TRUE), add = TRUE)` to clean up on any exit path. |
| BPCells not installed | `requireNamespace` check at entry; clear error message. On the dense path, BPCells is never touched. |
| `write_matrix_dir` of double data with `compress = TRUE` warns | Use `compress = FALSE` for correction output (scaled data is double-precision). |

---

## Implementation Order

1. Add `ondisk.path` and `chunk_size` parameters to the function signature and
   roxygen block.
2. Add the on-disk detection logic (`ondisk <- inherits(GEX_full, "IterableMatrix")`).
3. Gate the existing 2^31 size check to fire only when `!ondisk`.
4. Branch Stage 1 (regression loop): on the on-disk path, materialise per-batch
   via `as.matrix(t(GEX_full[, cells]))` instead of pre-transposing the full
   matrix.
5. Implement Stage 4 (correction) as the fully-streaming chunked approach.
6. Implement Stage 5 (store result): use `CreateAssay5Object` on the on-disk
   path.
7. Add cleanup logic (`on.exit` for chunk directories).
8. Add tests: a small BPCells-backed Seurat object through the full pipeline,
   comparing results against the dense path within floating-point tolerance.
