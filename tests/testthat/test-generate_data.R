library(matclust)

set.seed(1)
gd <- generate_data(n = 60, K = 2, p = 2, gom = 0.05, min_reps = 2, max_reps = 4)

test_that("generate_data returns a named list with three elements", {
  expect_named(gd, c("data", "class", "params"), ignore.order = FALSE)
})

test_that("data array has correct type and three dimensions", {
  expect_true(is.array(gd$data))
  expect_equal(length(dim(gd$data)), 3L)
})

test_that("data array third dimension equals n", {
  expect_equal(dim(gd$data)[3], 60L)
})

test_that("data array second dimension equals p", {
  expect_equal(dim(gd$data)[2], 2L)
})

test_that("data array first dimension is between min_reps and max_reps", {
  R <- dim(gd$data)[1]
  expect_gte(R, 2L)
  expect_lte(R, 4L)
})

test_that("class vector has length n and values in 1:K", {
  cl <- gd$class
  expect_length(cl, 60L)
  expect_true(all(cl %in% 1:2))
  expect_equal(length(unique(cl)), 2L)
})

test_that("params contains Mu (K x p) and S (p x p x K)", {
  pr <- gd$params
  expect_named(pr, c("Mu", "S"), ignore.order = FALSE)
  expect_equal(dim(pr$Mu), c(2L, 2L))
  expect_equal(dim(pr$S),  c(2L, 2L, 2L))
})

test_that("NA values have no observed data after them within each feature column", {
  dat <- gd$data
  R <- dim(dat)[1]; p <- dim(dat)[2]; n <- dim(dat)[3]
  for (i in 1:n) {
    for (j in 1:p) {
      col_vals <- dat[, j, i]
      na_pos   <- which(is.na(col_vals))
      if (length(na_pos) > 0)
        expect_true(all(is.na(col_vals[min(na_pos):R])))
    }
  }
})

test_that("at least one observation has missing replicates when min_reps < max_reps", {
  expect_true(anyNA(gd$data))
})

test_that("generate_data with min_reps == max_reps produces no NAs", {
  set.seed(2)
  gd_complete <- generate_data(n = 30, K = 2, p = 2, gom = 0.05,
                               min_reps = 3, max_reps = 3)
  expect_false(anyNA(gd_complete$data))
  expect_equal(dim(gd_complete$data)[1], 3L)
})

test_that("generate_data respects K=3", {
  set.seed(3)
  gd3 <- generate_data(n = 90, K = 3, p = 2, gom = 0.05)
  expect_equal(length(unique(gd3$class)), 3L)
  expect_equal(dim(gd3$params$Mu)[1], 3L)
  expect_equal(dim(gd3$params$S)[3], 3L)
})

test_that("generate_data is reproducible with the same seed", {
  set.seed(77)
  gd_a <- generate_data(n = 40, K = 2, p = 2, gom = 0.05)
  set.seed(77)
  gd_b <- generate_data(n = 40, K = 2, p = 2, gom = 0.05)
  expect_identical(gd_a$data,  gd_b$data)
  expect_identical(gd_a$class, gd_b$class)
})
