#' repclust
#'
#' A function to fit a finite mixture model of Gaussian distributions on
#' replicated data.
#'
#' @return A list
#' @export
#'
#' @examples
#' repclust()
repclust <- function(x, nclusters, iter_max = 101, tol = .001, init = "kmeans",
                     params = NULL,
                     emEM_args = list(nstarts = round(sqrt(nclusters*prod(dim(x)))),
                                      em_iter = 1, nbest = 10)) {
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
  } else if (init == "emEM") {
    if (emEM_args$nbest > emEM_args$nstarts) {
      message("Specified number of starts is less than `nbest`, setting `nstarts` to `nbest`.")
      emEM_args$nstarts <- emEM_args$nbest
    }

    message(
      paste0(
        "Determining best ", emEM_args$nbest, " values from ",
        emEM_args$nstarts, " short em runs assuming ", nclusters, " clusters."
      )
    )
    # best_ll <- -Inf
    best_lls <- rep(-Inf, emEM_args$nbest)
    # best_res <- NULL
    best_res <- rep(list(NA), emEM_args$nbest)
    cutoff <- -Inf
    for (i in seq_len(emEM_args$nstarts)) {
      tmp <- repclust(
        x, nclusters, iter_max = emEM_args$em_iter, tol,
        init = "kmeans"
      )

      if (i <= emEM_args$nbest) {
        best_res[[i]] <- tmp
        best_lls[i] <- tmp$ll[length(tmp$ll)]
      } else if (tmp$ll[length(tmp$ll)] > cutoff) {
        best_res[[emEM_args$nbest]] <- tmp
        best_lls[emEM_args$nbest] <- tmp$ll[length(tmp$ll)]
        ll_order <- order(best_lls)
        best_lls <- best_lls[ll_order]
        best_res <- best_res[ll_order]
      }

      cutoff <- min(best_lls, na.rm = TRUE)
      # NOTE: possibly change this to have best_lls and best_res store all of
      # the `nstarts` results, and only select the `nbest` ones after
    }

    best_long_ll <- -Inf
    best_long_res <- NULL

    for (i in seq_len(emEM_args$nbest)) {
      message(
        paste0(
          "Running long EM from best em runs: ", i, " of ",
          emEM_args$nbest, " total assuming ", nclusters,
          " clusters..."
        )
      )
      tmp <- repclust(
        x, nclusters, iter_max = iter_max, tol = tol,
        init = "given", params = best_res[[i]]
      )

      if (tmp$ll[length(tmp$ll)] > best_long_ll) {
        best_long_ll <- tmp$ll[length(tmp$ll)]
        best_long_res <- tmp
      }
    }

    return(best_long_res)
  }

  ll <- numeric(iter_max + 1)
  bic <- numeric(iter_max + 1)
  if (init != "given") {
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

    cl <- apply(z, 1, which.max) - 1
    # ll[1] <- get_ll(x, mu, Sigma, R, p, res$cluster - 1)
    # TODO: bandaid
    z <- model.matrix(~ 0 + x, data.frame(x = factor(res$cluster))) |> as.data.frame() |> as.matrix()
    ll[1] <- get_ll(x, mu, Sigma, R, p, z)
  } else {
    mu <- params$mu
    Sigma <- params$Sigma
    z <- params$z
    pr <- params$pi
    cl <- params$class - 1

    ll[1] <- params$ll[length(params$ll)]
    bic[1] <- params$bic[length(params$bic)]
  }

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

