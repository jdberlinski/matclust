#' Generate matrix-variate data
#'
#' @param n Number of observations to generate
#' @param K Number of clusters
#' @param p Number of dimensions
#' @param gom Generalized overlap of Maitra
#' @param min_reps Minimum number of repetitions per observation
#' @param max_reps Maximum number of repetitions per observation
#' @param int Vector containing lower and upper bounds for the centers of each distribution
#' @param mix_prob Vector of mixing proportions for each cluster
#' @param complete Proportion of fully observed cases (with `max_reps` rows of
#' fully observed data) to be present for each cluster
#' @param rep_prob Vector of probabilities for each value of min_reps:max_reps
#' to occur in each observation
#' @return A list, containing an array of data with observations on the 3rd axis
#' and a vector of class memberships
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom MixSim MixGOM
#' @examples
#' data <- generate_data(1000, 4, 10, 0.2)
#' res <- run_em(data$data, 4)
generate_data <- function(n, K, p, gom, min_reps = 1, max_reps = 6, int = c(-5, 5), mix_prob = rep(1, K),
  complete = 0.05, rep_prob = rep(1, length(min_reps:max_reps))) {
  mix_prob <- mix_prob / sum(mix_prob)

  params <- MixSim::MixGOM(goMega = gom, K = K, p = p, int = int)

  rep_prob <- rep_prob / sum(rep_prob)

  if (min_reps == max_reps) {
    n_reps <- matrix(max_reps, nrow = n, ncol = p)
  } else {
    n_reps <- matrix(sample(min_reps:max_reps, n*p, prob = rep_prob, replace = T), nrow = n)
  }
  classes <- sample(1:K, n, prob = mix_prob, replace = T)
  n_each <- table(classes)
  X <- Reduce(c, Map(
     function(i) {
       n_complete <- ceiling(complete * n_each[i])
       data <- list()
       if (i == 1) {
         start_ind <- 0
       } else {
         start_ind <- sum(n_each[1:(i-1)])
       }
       for (j in 1:n_each[i]) {
         nij <- n_reps[start_ind + j, ]
         if (j <= n_complete)
           nij <- rep(max_reps, p)
         # nij <- ifelse(j <= n_complete, rep(max_reps, p), nij)
         nmax <- max(nij)
         x <- MASS::mvrnorm(nmax, params$Mu[i, ], params$S[,,i])
         for (k in 1:p) {
           if (nij[k] < nmax)
             x[(nij[k] + 1):nmax, k] <- NA
         }
         data[[length(data) + 1]] <- list(data = x, class = i)
       }
       return(data)
     },
     1:K
  ))

  R <- max(sapply(X, function(e) ifelse(is.matrix(e$data), nrow(e$data), 1))) # |> max()
  cl <- sapply(X, function(e) e$class)
  X <- lapply(
    X,
    function(x) {
      if(!is.matrix(x$data))
        x$data <- matrix(x$data, nrow = 1)
      # n <- ifelse(is.matrix(x$data), nrow(data), 1)
      n <- nrow(x$data)
      if (n < R) {
        nw <- matrix(NA, nrow = R, ncol = p)
        nw[1:n, ] <- x$data
        x$data <- nw
      }
      return(x$data)
    }
  )
  X <- array(unlist(X), dim = c(R, p, n))
  X <- list(data = X, class = cl)

  return(X)
}

#' An Efficient Kernel K-Means Algorithm
#' @name kkmeans
#'
#' @description Performs kernel k-means with the specified kernel using an
#' optimal-transfer quick-transfer algorithm.
#'
#' @param data Numeric data to cluster. This will be converted to a matrix using `as.matrix`.
#' @param k Number of clusters.
#' @param kern Kernel to use, one of ('gaussian', 'poly').
#' @param param value of parameter to pass to kernel function.(eg sigma in
#' gaussian kernel). The Gaussian kernel is K(x, y) = exp(- ||x - y||^2 / (2*`param`))),
#' and the polynomial kernel is K(x, y) = (x'y + 1) ^ `param`
#' @param nstart Number of times to run the algorithm. The run with the lowest
#' total within cluster SSE (in feature space) will be returned
#' @param iter_max The maximum number of iterations to allow.
#' @param estimate If using the Gaussian kernel, specifying `estimate = "mknn"` will use
#' an `nn`-nearest neighbor method for estimating `param`.
#' @param nn How many neighbors to consider for mknn estimation.
#' @param init_centers The initial values for cluster membership. If `nstart` is greater
#' than 1, any start beyond the first iteration will use randomized centers.
#' @param method Which method to use for kernel k-means iteration. One of ("otqt", "macqueen", "lloyd").
#' "otqt" is a method using optimal-transfer and quick-transfer heuristics similar to the Hartigan and
#' Wong algorithm for k-means clustering.
#' @param trueest Whether or not the within-cluster sum of squares should be
#' recomputed in R after clustering is finished
#' @param kmat kernel matrix, if using a custom kernel
#' @return A list containing the following useful information
#' \describe{
#'   \item{cluster}{The final cluster membership.}
#'   \item{centers}{A k x p matrix, the rows of which contain the centers of the clusters in R^n (not to be confused
#' with the clusters in feature space)}
#'   \item{wss}{The within-cluster sum of squares for each cluster in feature space.}
#'   \item{param}{The parameter value used.}
#' }
#' @export
#' @examples
#' data <- as.matrix(iris[, 1:4])
