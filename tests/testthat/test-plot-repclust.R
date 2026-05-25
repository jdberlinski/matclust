library(matclust)

make_fit <- function(K = 2L, seed = 42L) {
  set.seed(seed)
  n_per <- 40L; p <- 2L; R <- 3L; n <- K * n_per
  x <- array(NA_real_, dim = c(R, p, n))
  for (k in seq_len(K))
    x[, , ((k - 1) * n_per + 1):(k * n_per)] <-
      rnorm(n_per * p * R, mean = (k - 1) * 20, sd = 0.5)
  repclust(x, nclusters = K, iter_max = 20, tol = 1e-4, init = "kmeans")
}

res2 <- make_fit(K = 2L)
res3 <- make_fit(K = 3L, seed = 7L)

test_that("repclust result has class 'repclust'", {
  expect_s3_class(res2, "repclust")
  expect_s3_class(res3, "repclust")
})

test_that("plot.repclust runs without error for K=2", {
  expect_no_error(plot(res2))
})

test_that("plot.repclust runs without error for K=3", {
  expect_no_error(plot(res3))
})

test_that("plot.repclust returns the input invisibly", {
  out <- plot(res2)
  expect_identical(out, res2)
})

test_that("plot.repclust resets par on exit", {
  old_mfrow <- par("mfrow")
  plot(res2)
  expect_equal(par("mfrow"), old_mfrow)
})

test_that("plot.repclust runs without error for K=1", {
  set.seed(1)
  x1 <- array(rnorm(3 * 2 * 30), dim = c(3L, 2L, 30L))
  res1 <- repclust(x1, nclusters = 1L, iter_max = 10)
  expect_s3_class(res1, "repclust")
  expect_no_error(plot(res1))
})
