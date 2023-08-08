#' run_em
#'
#' A function to run the EM algorithm on matrix-variate data
#'
#' @return A list
#' @export
#'
#' @examples
#' run_em()
run_em <- function(x, nclusters, iter_max = 101, tol = .001, init = "kmeans") {
  R <- dim(x)[1]
  p <- dim(x)[2]
  n <- dim(x)[3]
  K <- nclusters

  A <- apply(x, c(2, 3), is.na)

  # replace with 0 for cpp
  for (i in 1:n)
    for (j in 1:p)
      x[which(is.na(x[,j,i])),j,i] <- 0

  # obtain initial solution using K-means on first replicate for each data point
  firstx <- t(apply(x, 3, function(sl) sl[1,]))
  if (init == "kmeans") {
    res <- kmeans(firstx, nclusters, nstart = 10)
    mu <- res$centers
    n_k <- res$size
    pr <- res$size / sum(res$size)
  } else if (init == "random") {
    res <- list()
    res$cluster <- sample(1:nclusters, n, replace = T)
    mu <- matrix(0, nrow = K, ncol = p)
    for (k in 1:K)
      mu[k,] <- colMeans(firstx[res$cluster == k,])
    n_k <- table(res$cluster)
    pr <- n_k / sum(n_k)
  }

  z <- matrix(0, nrow = n, ncol = K)
  # each row contains a cluster mean
  # mu <- matrix(0, nrow = K, ncol = p)
  # each slice a pxp covariance matrix
  Sigma <- array(0, dim = c(p, p, K))

  for (k in 1:K) {
    for (i in 1:n) {
      if (res$cluster[i] == k) {
        Sigma[,,k] <- Sigma[,,k] + t(x[,,i] - rep(1, R) %*% t(mu[k,])) %*% (x[,,i] - rep(1, R) %*% t(mu[k,]))
      }
    }
    Sigma[,,k] <- Sigma[,,k] / (n_k[k] * R)
  }

  ll <- numeric(iter_max + 1)
  bic <- numeric(iter_max + 1)
  cl <- apply(z, 1, which.max) - 1
  # ll[1] <- get_ll(x, mu, Sigma, R, p, res$cluster - 1)
  # TODO: bandaid
  z <- model.matrix(~ 0 + x, data.frame(x = factor(res$cluster))) |> as.data.frame() |> as.matrix()
  ll[1] <- get_ll(x, mu, Sigma, R, p, z)

  # for testing initialization
  if (iter_max == 0) {
    return(
      list(z = z, pi = pr, mu = mu, Sigma = Sigma)
    )
  }

  # EM loop
  for (iter in 1:iter_max) {
    params <- em_step(x, mu, Sigma, z, pr, cl, A, n, K, R, p, iter)

    mu <- params$mu
    Sigma <- params$Sigma
    z <- params$z
    pr <- params$pr
    cl <- params$cl
    x <- params$x

    ll[iter + 1] <- params$ll
    bic[iter + 1] <- params$bic

    if(abs((ll[iter + 1] - ll[iter]) / ll[iter]) < tol)
      break
  }
  ll <- ll[1:iter + 1]
  bic <- bic[1:iter + 1]
  return(
    list(z = z, pi = pr, mu = mu, Sigma = Sigma, class = as.numeric(cl + 1), ll = ll, bic = bic)
  )
}

