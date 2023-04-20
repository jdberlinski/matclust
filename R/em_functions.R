#' run_em
#'
#' A function to run the EM algorithm on matrix-variate data
#'
#' @return A list
#' @export
#'
#' @examples
#' run_em()
run_em <- function(x, nclusters, iter_max = 100, tol = .001, init = "kmeans") {
  R <- dim(x)[1]
  p <- dim(x)[2]
  n <- dim(x)[3]
  K <- nclusters

  # replace with column means
  for (i in 1:n) {
    cm <- colMeans(x[,,i], na.rm = TRUE)
    for (j in 1:p) {
      x[which(is.na(x[,j,i])),j,i] <- cm[j]
    }
  }

  # obtain initial solution using K-means on first replicate for each data point
  firstx <- t(apply(x, 3, function(sl) sl[1,]))
  if (init == "kmeans") {
    res <- kmeans(firstx, nclusters, nstart = 10)
    mu <- res$centers
    n_k <- res$size
    pr <- res$size / sum(res$size)
  } else if (init == "random") {
    res <- list()
    res$cluster <- sample(1:4, n, replace = T)
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
  cl <- apply(z, 1, which.max)
  ll[1] <- get_ll(x, pr, mu, Sigma, R, p, res$cluster)

  # EM loop
  for (iter in 1:iter_max) {

    # E-step
    for (k in 1:K) {
      for (i in 1:n) {
        z[i, k] <- pr[k] * f_k(x[,,i], mu[k,], Sigma[,,k], R, p) / f(x[,,i], pr, mu, Sigma, R, p, K)
      }
    }

    # M-step
    for (k in 1:K) {
      acc <- 0
      sacc <- matrix(0, p, p)
      den <- sum(z[,k])
      for (i in 1:n) {
        acc <- acc + z[i,k] * colMeans(x[,,i])
        sacc <- sacc + z[i,k] * t(x[,,i] - rep(1, R) %*% t(mu[k,])) %*% (x[,,i] - rep(1, R) %*% t(mu[k,]))
      }
      mu[k,] <- acc / den
      Sigma[,,k] <- sacc / (R * den)
      pr[k] <- den / n
    }
    cl <- apply(z, 1, which.max)

    ll[iter + 1] <- get_ll(x, pr, mu, Sigma, R, p, cl)

    # print((ll[iter + 1] - ll[iter]) / ll[iter])
    if(abs((ll[iter + 1] - ll[iter]) / ll[iter]) < tol)
      break
  }
  ll <- ll[1:iter + 1]

  max_ll <- ll[length(ll)]
  bic <- -2*max_ll + log(n) * (K + K*p + K*p*(p+1)/2)

  # cl <- apply(z, 1, which.max)
  return(
    list(z = z, pi = pr, mu = mu, Sigma = Sigma, class = cl, ll = ll, bic = bic)
  )
}

f_k <- function(xi, mu, sig, R, p) {
  tryCatch(invisible(solve(sig)), error = function(e) print(sig))
  M <- rep(1, R) %*% t(mu)
  (2*pi)^(-R*p/2) * det(sig)^(-p/2) * exp(-.5 * sum(diag((xi - M) %*% solve(sig) %*% t(xi - M))))
}

f <- function(xi, pr, mu, sig, R, p, K) {
  acc <- 0
  for (k in 1:K)
    acc <- acc + pr[k] * f_k(xi, mu[k,], sig[,,k], R, p)
  return(acc)
}

get_ll <- function(x, pr, mu, sig, R, p, cl) {
  n <- dim(x)[3]
  ll <- 0
  for (i in 1:n)
    ll <- ll + log(f_k(x[,,i], mu[cl[i],], sig[,,cl[i]], R, p))
  return(ll)
}

