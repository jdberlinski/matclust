library(matclust)

make_ll_inputs <- function(n = 10, R = 2, p = 2, K = 2, seed = 42) {
  set.seed(seed)
  x   <- array(rnorm(R * p * n), dim = c(R, p, n))
  A   <- array(FALSE, dim = c(R, p, n))
  mu  <- matrix(c(0, 0, 1, 1), nrow = K, ncol = p, byrow = TRUE)
  Sig <- array(0, dim = c(p, p, K))
  for (k in 1:K) Sig[, , k] <- diag(p)
  z   <- matrix(1 / K, nrow = n, ncol = K)
  list(x = x, A = A, mu = mu, Sig = Sig, z = z, R = R, p = p)
}

# ---- Scalar output ----

test_that("get_ll returns a single finite negative number", {
  inp <- make_ll_inputs()
  ll  <- get_ll(inp$x, inp$mu, inp$Sig, inp$A, inp$R, inp$p, inp$z)
  expect_length(ll, 1L)
  expect_true(is.finite(ll))
  expect_lt(ll, 0)
})

# ---- Analytical cross-checks ----

test_that("get_ll matches analytic log-density: R=1, p=1, xi=mu=0, Sigma=1", {
  x   <- array(0.0, dim = c(1L, 1L, 1L))
  A   <- array(FALSE, dim = c(1L, 1L, 1L))
  mu  <- matrix(0, nrow = 1L, ncol = 1L)
  Sig <- array(1.0, dim = c(1L, 1L, 1L))
  z   <- matrix(1.0, nrow = 1L, ncol = 1L)
  ll  <- get_ll(x, mu, Sig, A, 1L, 1L, z)
  expect_equal(ll, dnorm(0, log = TRUE), tolerance = 1e-7)
})

test_that("get_ll matches analytic log-density: R=2, p=1, both replicates at mu", {
  mu_val <- 3.0
  x   <- array(mu_val, dim = c(2L, 1L, 1L))
  A   <- array(FALSE,  dim = c(2L, 1L, 1L))
  mu  <- matrix(mu_val, nrow = 1L, ncol = 1L)
  Sig <- array(1.0, dim = c(1L, 1L, 1L))
  z   <- matrix(1.0, nrow = 1L, ncol = 1L)
  ll  <- get_ll(x, mu, Sig, A, 2L, 1L, z)
  expect_equal(ll, 2 * dnorm(0, log = TRUE), tolerance = 1e-7)
})

test_that("get_ll ignores replicate rows where A is TRUE (missing)", {
  # obs 1: replicate 1 observed (x=0, mu=0), replicate 2 missing (x=99, masked by A)
  # obs 2: both replicates observed (x=0, mu=0)
  # Expected total: dnorm(0,log=T) + 2*dnorm(0,log=T) = 3*dnorm(0,log=T)
  x          <- array(0.0,   dim = c(2L, 1L, 2L))
  x[2, 1, 1] <- 99.0
  A          <- array(FALSE, dim = c(2L, 1L, 2L))
  A[2, 1, 1] <- TRUE
  mu  <- matrix(0, nrow = 1L, ncol = 1L)
  Sig <- array(1.0, dim = c(1L, 1L, 1L))
  z   <- matrix(1.0, nrow = 2L, ncol = 1L)
  ll  <- get_ll(x, mu, Sig, A, 2L, 1L, z)
  expect_equal(ll, 3 * dnorm(0, log = TRUE), tolerance = 1e-7)
})

test_that("get_ll with full obs matches manual multivariate log-density computation", {
  set.seed(5)
  R <- 1L; p <- 2L; n <- 5L; K <- 2L
  x   <- array(rnorm(R * p * n), dim = c(R, p, n))
  A   <- array(FALSE, dim = c(R, p, n))
  mu  <- matrix(c(0, 0, 1, 1), nrow = K, ncol = p, byrow = TRUE)
  Sig <- array(0, dim = c(p, p, K))
  for (k in 1:K) Sig[, , k] <- diag(p)
  z   <- matrix(c(rep(c(1, 0), 3), rep(c(0, 1), 2)), nrow = n, ncol = K, byrow = TRUE)

  ll_val <- get_ll(x, mu, Sig, A, R, p, z)

  manual_ll <- 0
  for (i in 1:n) {
    k_used <- which.max(z[i, ])
    diff   <- matrix(x[1, , i], nrow = 1) - mu[k_used, , drop = FALSE]
    S_k    <- Sig[, , k_used]
    manual_ll <- manual_ll +
      -0.5 * p * log(2 * pi) -
      0.5 * log(det(S_k)) -
      0.5 * as.numeric(diff %*% solve(S_k) %*% t(diff))
  }
  expect_equal(ll_val, manual_ll, tolerance = 1e-7)
})

# ---- log_f_k (non-exported, accessed via :::) ----

test_that("log_f_k for R=1, p=1, xi=mu=0, sig=1 equals dnorm(0, log=TRUE)", {
  result <- matclust:::log_f_k(
    xi  = matrix(0, 1, 1),
    mu  = matrix(0, 1, 1),
    sig = matrix(1, 1, 1),
    R = 1L, p = 1L
  )
  expect_equal(result, dnorm(0, log = TRUE), tolerance = 1e-7)
})

test_that("log_f_k for R=2, p=2 returns a single finite number < 0", {
  set.seed(6)
  result <- matclust:::log_f_k(
    xi  = matrix(rnorm(4), 2, 2),
    mu  = matrix(c(0, 0), 1, 2),
    sig = diag(2),
    R = 2L, p = 2L
  )
  expect_length(result, 1L)
  expect_true(is.finite(result))
  expect_lt(result, 0)
})

test_that("log_f_k is maximized when xi equals the mean matrix", {
  mu  <- matrix(c(1.5, -0.5), 1, 2)
  sig <- diag(2)
  xi_at_mean <- matrix(rep(as.numeric(mu), 2), nrow = 2, byrow = TRUE)
  xi_off      <- xi_at_mean + 1.0
  expect_gt(
    matclust:::log_f_k(xi_at_mean, mu, sig, R = 2L, p = 2L),
    matclust:::log_f_k(xi_off,     mu, sig, R = 2L, p = 2L)
  )
})

test_that("log_f_k increases as sigma grows for an off-center observation", {
  xi  <- matrix(c(2, 2), 1, 2)
  mu  <- matrix(c(0, 0), 1, 2)
  ll_small <- matclust:::log_f_k(xi, mu, 0.5 * diag(2), R = 1L, p = 2L)
  ll_large <- matclust:::log_f_k(xi, mu, 5.0 * diag(2), R = 1L, p = 2L)
  expect_gt(ll_large, ll_small)
})
