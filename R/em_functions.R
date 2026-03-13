#' Gaussian mixture models for replicated data
#'
#' A function to fit a finite mixture model of Gaussian distributions on
#' replicated data. The data are not required to have an equal number of
#' replicates for each observation, or each feature within each observation.
#'
#' A single observation may be of the form:
#' \preformatted{
#'   x_i11, x_i21, x_i31, x_i41
#'   x_i12,   -  , x_i31, x_i41
#'   x_i13,   -  , x_i31,   -
#'     -  ,   -  , x_i31,   -
#' }
#'
#' The finite mixture model assumes that if the the observations were fully
#' observed (meaning all features had R available replicates), they follow the distribution
#'
#' \deqn{
#' X_i^* \sim \sum_{k=1}^K \pi_k f(X_i^*;\; \mu_k, \Sigma_k).
#' }{X_i^* ~ sum_{k=1}^k pi_k * f(X_i^*; mu_k, Sigma_k).}
#'
#' where \eqn{f}{f} is a matrix-variate normal distribution. The mean matrix has
#' equal columns equal to \eqn{\mu_k}{mu_k}, column covariance equal to
#' \eqn{\Sigma_k}{Sigma_k}, and identity row covariance.
#'
#' This inherently assumes that replicates are independent within each subject.
#'
#' @param x The array of data to be clustered. Should have dimensions (R, p, n),
#' where R is the maximum ' number of replicates, p is the dimensionality of the
#' data, and n is the number of observations. See the output of
#' `generate_data()` for examples
#' @param nclusters numeric value, the desired number of clusters.
#' @param iter_max numeric maximum number of EM steps to conduct
#' @param tol numeric tolerance for difference in likelihood at successive EM
#' steps. When two successive steps differ in less than this value the algorithm
#' terminates
#' @param init string corresponding to the initialization method. Options are
#' "kmeans," ' which calls `kmeans()` on the first row of each observation,
#' "random," which randomly ' assigns observations to clusters, or "emEM," which
#' performs a specified number of short  initial EM runs specified in the
#' `emEM_args` parameter, or "given," which uses the inital ' parameters specified
#' in `params`.
#' @param params list of initial parameter values if `init == "given"`. The list should have values:
#' \describe{
#'   \item{`mu`}{A matrix with `nclusters` rows and p columns representing the mean vectors for each cluster}
#'   \item{`Sigma`}{An array with dimensions (p, p, `nclusters`) corresponding to the covariance matrices for each cluster}
#'   \item{`z`}{A matrix with n rows and `nclusters` columns that corresponds to the probabilty of each point belonging in each cluster}
#'   \item{`pi`}{An `nclusters` dimensional vector with values corresponding to the mixing proportions of each cluster.}
#'   \item{`class`}{An n dimensional vector with the assigned cluster for each observation.}
#' }
#' @param emEM_args a list of settings for the "emEM" initialization option. The list should have values:
#' \describe{
#'   \item{`nstarts`}{The number of short EM runs to conduct.}
#'   \item{`em_iter`}{The number of EM steps to conduct in each run}
#'   \item{`nbest`}{The number of short runs to fully run. The returned solution will be the best of these runs.}
#' }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{`z`}{The matrix of class probabilities for each observation}
#'   \item{`pi`}{The mixing proportions for each cluster}
#'   \item{`mu`}{The mean vectors for each cluster}
#'   \item{`Sigma`}{The covariance matrices for each cluster}
#'   \item{`class`}{The assigned class for each observation (with the highest value in `z`)}
#'   \item{`ll`}{The sequence evaluted log-likelihood values. The length is equal to the number of EM runs that were conducted.}
#'   \item{`bic`}{The Bayesian information critereon associated with each of the values in `ll`}
#' }
#'
#' @export
#'
#' @examples
#' simulated_data <- generate_data(1000, 4, 10, 0.2)
#' res <- repclust(simulated_data$data, 4)
repclust <- function(
  x,
  nclusters,
  iter_max = 101,
  tol = .001,
  init = "kmeans",
  params = NULL,
  sigma_constr = "VVV",
  emEM_args = list(
    nstarts = round(sqrt(nclusters * prod(dim(x)))),
    em_iter = 1,
    nbest = 10
  )
) {
  R <- dim(x)[1]
  p <- dim(x)[2]
  n <- dim(x)[3]
  K <- nclusters

  A <- apply(x, c(2, 3), is.na)

  # replace with 0 for cpp
  for (i in 1:n) {
    for (j in 1:p) {
      x[which(is.na(x[, j, i])), j, i] <- 0
    }
  }

  # obtain initial solution using K-means on first replicate for each data point
  firstx <- t(apply(x, 3, function(sl) sl[1, ]))
  if (init == "kmeans") {
    res <- kmeans(firstx, nclusters, nstart = 10)
    mu <- res$centers
    n_k <- res$size
    pr <- res$size / sum(res$size)
  } else if (init == "random") {
    res <- list()
    res$cluster <- sample(1:nclusters, n, replace = T)
    mu <- matrix(0, nrow = K, ncol = p)
    for (k in 1:K) {
      mu[k, ] <- colMeans(firstx[res$cluster == k, ])
    }
    n_k <- table(res$cluster)
    pr <- n_k / sum(n_k)
  } else if (init == "emEM") {
    if (emEM_args$nbest > emEM_args$nstarts) {
      message(
        "Specified number of starts is less than `nbest`, setting `nstarts` to `nbest`."
      )
      emEM_args$nstarts <- emEM_args$nbest
    }

    message(
      paste0(
        "Determining best ",
        emEM_args$nbest,
        " values from ",
        emEM_args$nstarts,
        " short em runs assuming ",
        nclusters,
        " clusters."
      )
    )
    # best_ll <- -Inf
    best_lls <- rep(-Inf, emEM_args$nbest)
    # best_res <- NULL
    best_res <- rep(list(NA), emEM_args$nbest)
    cutoff <- -Inf
    for (i in seq_len(emEM_args$nstarts)) {
      # tmp <- repclust(
      #   x,
      #   nclusters,
      #   iter_max = emEM_args$em_iter,
      #   tol,
      #   init = "kmeans"
      # )

      tmp <-
        tryCatch({
          repclust(
            x,
            nclusters,
            iter_max = emEM_args$em_iter,
            tol = tol,
            init = "kmeans"
          )},
          error = function(e) return(NA)
        )

      if (!is.list(tmp)) next

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
          "Running long EM from best em runs: ",
          i,
          " of ",
          emEM_args$nbest,
          " total assuming ",
          nclusters,
          " clusters..."
        )
      )
      tmp <-
        tryCatch({
          repclust(
            x,
            nclusters,
            iter_max = iter_max,
            tol = tol,
            init = "given",
            params = best_res[[i]]
          )},
          error = function(e) return(NA)
        )

      if (!is.list(tmp)) {
        message("\tEncountered error, skipping run.")
        next
      }

      if (tmp$ll[length(tmp$ll)] > best_long_ll) {
        best_long_ll <- tmp$ll[length(tmp$ll)]
        best_long_res <- tmp
        best_long_res$ll <- c(best_res[[i]]$ll, tmp$ll)
        best_long_res$bic <- c(best_res[[i]]$bic, tmp$bic)
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
          Sigma[,, k] <- Sigma[,, k] +
            t(x[,, i] - rep(1, R) %*% t(mu[k, ])) %*%
              (x[,, i] - rep(1, R) %*% t(mu[k, ]))
        }
      }
      Sigma[,, k] <- Sigma[,, k] / (n_k[k] * R)
    }

    cl <- apply(z, 1, which.max) - 1
    # ll[1] <- get_ll(x, mu, Sigma, R, p, res$cluster - 1)
    # TODO: bandaid
    z <- model.matrix(~ 0 + x, data.frame(x = factor(res$cluster))) |>
      as.data.frame() |>
      as.matrix()
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
    params <- em_step(x, mu, Sigma, z, pr, cl, A, n, K, R, p, iter, sigma_constr)

    mu <- params$mu
    Sigma <- params$Sigma
    z <- params$z
    pr <- params$pr
    cl <- params$cl
    x <- params$x

    ll[iter + 1] <- params$ll
    bic[iter + 1] <- params$bic

    if (is.na(abs((ll[iter + 1] - ll[iter]) / ll[iter])))
      print(ll)

    if (abs((ll[iter + 1] - ll[iter]) / ll[iter]) < tol) {
      break
    }
  }
  ll <- ll[1:iter + 1]
  bic <- bic[1:iter + 1]
  return(
    list(
      z = z,
      pi = pr,
      mu = mu,
      Sigma = Sigma,
      class = as.numeric(cl + 1),
      ll = ll,
      bic = bic
    )
  )
}

