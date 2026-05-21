library(matclust)

# Well-separated 2-cluster fixture; no MixSim dependency
make_sep_data <- function() {
  set.seed(42)
  n_per <- 50L; p <- 2L; R <- 3L; n <- 2L * n_per
  x <- array(NA_real_, dim = c(R, p, n))
  x[, , 1:n_per]               <- rnorm(n_per * p * R, mean =  10, sd = 0.3)
  x[, , (n_per + 1):(2*n_per)] <- rnorm(n_per * p * R, mean = -10, sd = 0.3)
  list(x = x, true_class = c(rep(1L, n_per), rep(2L, n_per)))
}

sep    <- make_sep_data()
x_sep  <- sep$x
tc_sep <- sep$true_class
n_sep  <- dim(x_sep)[3]  # 100
p_sep  <- dim(x_sep)[2]  # 2
K_sep  <- 2L

set.seed(42)
res <- repclust(x_sep, nclusters = K_sep, iter_max = 30, tol = 1e-4, init = "kmeans")

# ---- Output structure ----

test_that("repclust returns a named list with all expected elements", {
  expect_named(res, c("z", "pi", "mu", "Sigma", "class", "ll", "bic"),
               ignore.order = FALSE)
})

test_that("z is an n x K numeric matrix with values in [0, 1]", {
  z <- res$z
  expect_true(is.matrix(z))
  expect_equal(dim(z), c(n_sep, K_sep))
  expect_true(all(z >= -1e-10) && all(z <= 1 + 1e-10))
})

test_that("rows of z sum to 1 (posterior probabilities)", {
  expect_true(all(abs(rowSums(res$z) - 1) < 1e-10))
})

test_that("pi is a length-K numeric vector summing to 1 with values in (0, 1)", {
  pi_vec <- res$pi
  expect_length(pi_vec, K_sep)
  expect_equal(sum(pi_vec), 1, tolerance = 1e-10)
  expect_true(all(pi_vec > 0) && all(pi_vec < 1))
})

test_that("mu is a K x p numeric matrix", {
  expect_true(is.matrix(res$mu))
  expect_equal(dim(res$mu), c(K_sep, p_sep))
})

test_that("Sigma is a p x p x K array with symmetric positive-definite slices", {
  Sigma <- res$Sigma
  expect_true(is.array(Sigma))
  expect_equal(dim(Sigma), c(p_sep, p_sep, K_sep))
  for (k in 1:K_sep) {
    S <- Sigma[, , k]
    expect_equal(S, t(S), tolerance = 1e-10,
                 info = paste("cluster", k, "Sigma not symmetric"))
    expect_true(all(eigen(S, only.values = TRUE)$values > 0),
                info = paste("cluster", k, "Sigma not positive definite"))
  }
})

test_that("class is a length-n vector with values in 1:K", {
  cl <- res$class
  expect_length(cl, n_sep)
  expect_true(all(cl %in% 1:K_sep))
})

test_that("class[i] equals which.max(z[i, ]) for all i", {
  expect_equal(as.integer(res$class), apply(res$z, 1, which.max))
})

test_that("ll is numeric with length >= 1 and bic has the same length", {
  expect_true(is.numeric(res$ll))
  expect_gte(length(res$ll), 1L)
  expect_equal(length(res$bic), length(res$ll))
})

test_that("ll is non-decreasing (EM convergence property)", {
  ll <- res$ll
  if (length(ll) > 1) {
    tol_mono <- 1e-6 * abs(ll[1])
    expect_true(all(diff(ll) >= -tol_mono),
                info = paste("largest decrease:", min(diff(ll))))
  }
})

# ---- Cluster recovery ----

test_that("repclust recovers well-separated clusters with >= 95% accuracy", {
  cl  <- as.integer(res$class)
  tc  <- tc_sep
  acc <- max(mean(cl == tc), mean(cl == (3L - tc)))  # try both label permutations
  expect_gte(acc, 0.95, label = sprintf("best accuracy = %.3f", acc))
})

test_that("recovered cluster means are close to true means (+/-10)", {
  mu <- res$mu
  row_pos <- which.max(mu[, 1])
  row_neg <- 3L - row_pos
  expect_equal(mu[row_pos, 1], 10, tolerance = 0.5)
  expect_equal(mu[row_neg, 1], -10, tolerance = 0.5)
})

# ---- Initialization variants ----

test_that("init='random' produces a valid result", {
  set.seed(99)
  r <- repclust(x_sep, nclusters = 2L, iter_max = 20, init = "random")
  expect_named(r, c("z", "pi", "mu", "Sigma", "class", "ll", "bic"), ignore.order = FALSE)
  expect_equal(sum(r$pi), 1, tolerance = 1e-10)
  expect_true(all(abs(rowSums(r$z) - 1) < 1e-10))
})

test_that("init='given' continues from a provided parameter set", {
  set.seed(7)
  r_short <- repclust(x_sep, nclusters = 2L, iter_max = 5, init = "kmeans")
  set.seed(7)
  r_given <- repclust(x_sep, nclusters = 2L, iter_max = 10, init = "given",
                      params = r_short)
  expect_named(r_given, c("z", "pi", "mu", "Sigma", "class", "ll", "bic"), ignore.order = FALSE)
  expect_equal(sum(r_given$pi), 1, tolerance = 1e-10)
  expect_gte(r_given$ll[length(r_given$ll)],
             r_short$ll[length(r_short$ll)] - 1e-4)
})

# ---- Covariance constraints ----

test_that("sph=TRUE (VII) produces diagonal Sigma slices", {
  set.seed(11)
  r <- repclust(x_sep, nclusters = 2L, iter_max = 20, sph = TRUE, hom = FALSE)
  for (k in 1:2L) {
    S <- r$Sigma[, , k]
    expect_true(all(abs(S[!diag(nrow = p_sep, ncol = p_sep)]) < 1e-10),
                info = paste("cluster", k, "Sigma not diagonal"))
  }
})

test_that("sph=TRUE, hom=TRUE (EII) produces equal Sigma slices", {
  set.seed(12)
  r <- repclust(x_sep, nclusters = 2L, iter_max = 20, sph = TRUE, hom = TRUE)
  expect_equal(r$Sigma[, , 1], r$Sigma[, , 2], tolerance = 1e-10)
})

test_that("sph=FALSE, hom=TRUE (EEE) produces equal Sigma slices", {
  set.seed(13)
  r <- repclust(x_sep, nclusters = 2L, iter_max = 20, sph = FALSE, hom = TRUE)
  expect_equal(r$Sigma[, , 1], r$Sigma[, , 2], tolerance = 1e-10)
})

# ---- Edge cases ----

test_that("repclust with K=1 assigns all observations to cluster 1", {
  set.seed(55)
  x1 <- array(rnorm(3 * 2 * 30), dim = c(3L, 2L, 30L))
  r <- repclust(x1, nclusters = 1L, iter_max = 10)
  expect_true(all(r$class == 1L))
  expect_equal(as.numeric(r$pi), 1, tolerance = 1e-10)
  expect_equal(dim(r$z), c(30L, 1L))
})

test_that("repclust on complete data (no NAs) runs without error", {
  set.seed(33)
  x_comp <- array(rnorm(4 * 3 * 40), dim = c(4L, 3L, 40L))
  expect_no_error(repclust(x_comp, nclusters = 2L, iter_max = 15))
})

test_that("repclust with some missing-replicate rows runs without error", {
  x_miss <- x_sep
  x_miss[2, , 1:5] <- NA
  x_miss[3, , 1:5] <- NA
  expect_no_error(repclust(x_miss, nclusters = 2L, iter_max = 10))
})

test_that("repclust with tol=0 runs all iter_max steps", {
  set.seed(22)
  r <- repclust(x_sep, nclusters = 2L, iter_max = 5L, tol = 0, init = "kmeans")
  expect_equal(length(r$ll), 5L)
})
