# Tests for on-disk Sparse PLS (SPLS) implementation

# Helper: convert n x p dense matrix to BPCells IterableMatrix in features x samples layout.
to_bp <- function(X_np) {
  as(as(t(X_np), "dgCMatrix"), "IterableMatrix")
}

# Helper: compare on-disk SPLS result to in-memory spls::spls() reference.
# Because on-disk SPLS centers but does not scale (unlike spls::spls which
# scales by default), we call spls::spls(..., scale.x = FALSE, scale.y = FALSE)
# for a fair comparison.
compare_spls_results <- function(bp, ref, tol = 1e-6) {
  # Active set should match
  expect_equal(sort(bp$A), sort(ref$A))

  # betahat should match (p x q)
  expect_equal(bp$betahat, ref$betahat, tolerance = tol)

  # Projection dimensions should match: |A| x ncomp_actual
  expect_equal(dim(bp$projection), dim(ref$projection))

  # Projection may have per-component sign flips; strip names to avoid
  # mismatches from differing formula variable names (e.g. "t(X_A_centered)1"
  # vs "xA1").
  bp_proj <- unname(bp$projection)
  ref_proj <- unname(ref$projection)
  for (a in seq_len(ncol(bp_proj))) {
    bp_col <- bp_proj[, a]
    ref_col <- ref_proj[, a]
    idx <- which.max(abs(ref_col))
    s <- sign(bp_col[idx]) * sign(ref_col[idx])
    expect_equal(s * bp_col, ref_col, tolerance = tol,
                 label = sprintf("projection component %d", a))
  }
}


test_that("on-disk SPLS matches in-memory spls: single response", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("spls")
  skip_if_not_installed("pls")

  set.seed(42)
  n <- 80
  p <- 15
  ncomp <- 5

  X <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(2, 5), rep(0, 10))
  Y <- X %*% beta + rnorm(n, sd = 0.5)
  Y <- as.matrix(Y)

  # In-memory reference (no scaling, to match on-disk behavior)
  ref <- spls::spls(x = X, y = Y, K = ncomp, eta = 0.5,
                     scale.x = FALSE, scale.y = FALSE, fit = "kernelpls")

  X_bp <- to_bp(X)
  bp <- spls_ondisk(X_bp, Y, ncomp = ncomp, eta = 0.5, fit = "kernelpls")

  compare_spls_results(bp, ref)
})


test_that("on-disk SPLS matches in-memory spls: multiple responses", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("spls")
  skip_if_not_installed("pls")

  set.seed(123)
  n <- 100
  p <- 20
  q <- 3
  ncomp <- 5

  X <- matrix(rnorm(n * p), n, p)
  beta <- matrix(0, p, q)
  beta[1:8, ] <- matrix(rnorm(8 * q), 8, q)
  Y <- X %*% beta + matrix(rnorm(n * q, sd = 0.3), n, q)

  ref <- spls::spls(x = X, y = Y, K = ncomp, eta = 0.5,
                     scale.x = FALSE, scale.y = FALSE, fit = "kernelpls")

  X_bp <- to_bp(X)
  bp <- spls_ondisk(X_bp, Y, ncomp = ncomp, eta = 0.5, fit = "kernelpls")

  compare_spls_results(bp, ref)
})


test_that("on-disk SPLS matches in-memory spls with high sparsity (eta=0.9)", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("spls")
  skip_if_not_installed("pls")

  set.seed(7)
  n <- 60
  p <- 20
  ncomp <- 3

  X <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(3, 3), rep(0, 17))
  Y <- X %*% beta + rnorm(n, sd = 0.3)
  Y <- as.matrix(Y)

  ref <- spls::spls(x = X, y = Y, K = ncomp, eta = 0.9,
                     scale.x = FALSE, scale.y = FALSE, fit = "kernelpls")

  X_bp <- to_bp(X)
  bp <- spls_ondisk(X_bp, Y, ncomp = ncomp, eta = 0.9, fit = "kernelpls")

  compare_spls_results(bp, ref)
})


test_that("on-disk SPLS matches in-memory spls with low sparsity (eta=0.1)", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("spls")
  skip_if_not_installed("pls")

  set.seed(99)
  n <- 60
  p <- 12
  ncomp <- 3

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)

  ref <- spls::spls(x = X, y = Y, K = ncomp, eta = 0.1,
                     scale.x = FALSE, scale.y = FALSE, fit = "kernelpls")

  X_bp <- to_bp(X)
  bp <- spls_ondisk(X_bp, Y, ncomp = ncomp, eta = 0.1, fit = "kernelpls")

  compare_spls_results(bp, ref)
})


test_that("on-disk SPLS works with ncomp = 1", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("spls")
  skip_if_not_installed("pls")

  set.seed(314)
  n <- 50
  p <- 10

  X <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(2, 3), rep(0, 7))
  Y <- X %*% beta + rnorm(n, sd = 0.5)
  Y <- as.matrix(Y)

  ref <- spls::spls(x = X, y = Y, K = 1, eta = 0.5,
                     scale.x = FALSE, scale.y = FALSE, fit = "kernelpls")

  X_bp <- to_bp(X)
  bp <- spls_ondisk(X_bp, Y, ncomp = 1, eta = 0.5, fit = "kernelpls")

  compare_spls_results(bp, ref)
})


test_that("on-disk SPLS input validation: bad eta", {
  skip_if_not_installed("BPCells")

  X <- matrix(rnorm(30), 10, 3)
  X_bp <- to_bp(X)
  Y <- matrix(rnorm(10), 10, 1)

  expect_error(spls_ondisk(X_bp, Y, ncomp = 1, eta = -0.1), "eta")
  expect_error(spls_ondisk(X_bp, Y, ncomp = 1, eta = 1.0), "eta")
  expect_error(spls_ondisk(X_bp, Y, ncomp = 1, eta = 1.5), "eta")
})


test_that("on-disk SPLS input validation: bad kappa", {
  skip_if_not_installed("BPCells")

  X <- matrix(rnorm(30), 10, 3)
  X_bp <- to_bp(X)
  Y <- matrix(rnorm(10), 10, 1)

  expect_error(spls_ondisk(X_bp, Y, ncomp = 1, kappa = 0), "kappa")
  expect_error(spls_ondisk(X_bp, Y, ncomp = 1, kappa = 0.6), "kappa")
})


test_that("on-disk SPLS input validation: bad ncomp", {
  skip_if_not_installed("BPCells")

  X <- matrix(rnorm(30), 10, 3)
  X_bp <- to_bp(X)
  Y <- matrix(rnorm(10), 10, 1)

  expect_error(spls_ondisk(X_bp, Y, ncomp = 0), "ncomp")
  expect_error(spls_ondisk(X_bp, Y, ncomp = -1), "ncomp")
})


test_that("on-disk SPLS input validation: dimension mismatch", {
  skip_if_not_installed("BPCells")

  X <- matrix(rnorm(30), 10, 3)
  X_bp <- to_bp(X)
  Y_wrong <- matrix(rnorm(5), 5, 1)

  expect_error(spls_ondisk(X_bp, Y_wrong, ncomp = 1), "must match")
})


test_that("on-disk SPLS output dimensions are correct", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("pls")

  set.seed(555)
  n <- 50
  p <- 15
  q <- 2
  ncomp <- 3

  X <- matrix(rnorm(n * p), n, p)
  beta <- matrix(0, p, q)
  beta[1:5, ] <- matrix(rnorm(5 * q), 5, q)
  Y <- X %*% beta + matrix(rnorm(n * q, sd = 0.3), n, q)

  X_bp <- to_bp(X)
  result <- spls_ondisk(X_bp, Y, ncomp = ncomp, eta = 0.5, fit = "kernelpls")

  # betahat: p x q
  expect_equal(dim(result$betahat), c(p, q))

  # A is a set of indices in 1:p
  expect_true(all(result$A >= 1 & result$A <= p))

  # projection: |A| x ncomp_actual
  expect_equal(nrow(result$projection), length(result$A))

  # Xmeans: length p
  expect_equal(length(result$Xmeans), p)

  # Ymeans: length q
  expect_equal(length(result$Ymeans), q)

  # fitted.values: n x q x ncomp
  expect_equal(dim(result$fitted.values), c(n, q, ncomp))
})


test_that("on-disk SPLS works with sparse BPCells matrix", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("spls")
  skip_if_not_installed("pls")

  set.seed(456)
  n <- 100
  p <- 20
  ncomp <- 3

  X <- matrix(rnorm(n * p), n, p)
  X[abs(X) < 1] <- 0  # sparsify
  beta <- c(rep(2, 5), rep(0, 15))
  Y <- X %*% beta + rnorm(n, sd = 0.5)
  Y <- as.matrix(Y)

  ref <- spls::spls(x = X, y = Y, K = ncomp, eta = 0.5,
                     scale.x = FALSE, scale.y = FALSE, fit = "kernelpls")

  X_bp <- to_bp(X)
  bp <- spls_ondisk(X_bp, Y, ncomp = ncomp, eta = 0.5, fit = "kernelpls")

  compare_spls_results(bp, ref)
})


test_that("on-disk SPLS betahat is sparse", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("pls")

  set.seed(42)
  n <- 80
  p <- 30
  ncomp <- 3

  X <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(2, 3), rep(0, 27))
  Y <- X %*% beta + rnorm(n, sd = 0.5)
  Y <- as.matrix(Y)

  X_bp <- to_bp(X)
  result <- spls_ondisk(X_bp, Y, ncomp = ncomp, eta = 0.8, fit = "kernelpls")

  # With high eta, most betahat rows should be zero
  n_nonzero_rows <- sum(rowSums(abs(result$betahat)) > 0)
  expect_true(n_nonzero_rows < p)
  expect_true(n_nonzero_rows > 0)
})


test_that("ust helper works correctly", {
  # Basic soft thresholding
  b <- c(0.8, 0.3, -0.5, 0.1, -0.9)
  eta <- 0.5

  result <- scSLIDE:::.ust(b, eta)

  # max(|b|) = 0.9, threshold = 0.45
  # 0.8 - 0.45 = 0.35 > 0, result[1] = 0.35
  # 0.3 - 0.45 < 0, result[2] = 0
  # 0.5 - 0.45 = 0.05 > 0, result[3] = -0.05
  # 0.1 - 0.45 < 0, result[4] = 0
  # 0.9 - 0.45 = 0.45 > 0, result[5] = -0.45
  expect_equal(result[1, 1], 0.35, tolerance = 1e-10)
  expect_equal(result[2, 1], 0, tolerance = 1e-10)
  expect_equal(result[3, 1], -0.05, tolerance = 1e-10)
  expect_equal(result[4, 1], 0, tolerance = 1e-10)
  expect_equal(result[5, 1], -0.45, tolerance = 1e-10)
})


test_that("ust helper: eta = 0 returns all elements", {
  b <- c(1, 2, 3)
  result <- scSLIDE:::.ust(b, 0)
  # eta = 0: threshold = 0, all elements survive unchanged
  expect_equal(result[, 1], b)
})


test_that("spls_dv helper: single response is just ust", {
  Z <- matrix(c(0.5, -0.3, 0.8, 0.1), ncol = 1)
  result <- scSLIDE:::.spls_dv(Z, eta = 0.5, kappa = 0.5, eps = 1e-4, maxstep = 100)
  # For q=1, spls_dv normalizes Z by median(|Z|) and applies ust
  # median(|Z/median(|Z|)|) ... just check it's a valid sparse vector

  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 1)
  # Some elements should be zero with eta=0.5
  expect_true(any(result == 0))
})


test_that("on-disk SPLS message about no scaling is printed", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("pls")

  set.seed(42)
  n <- 30
  p <- 8

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n), n, 1)
  X_bp <- to_bp(X)

  expect_message(
    spls_ondisk(X_bp, Y, ncomp = 2, eta = 0.5),
    "scaling is not"
  )
})
