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

