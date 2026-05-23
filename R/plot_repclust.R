#' Diagnostic plots for a repclust fit
#'
#' Displays a 2×2 grid of model-fit diagnostics that do not require access to
#' the original (potentially multi-dimensional) data:
#' \enumerate{
#'   \item EM log-likelihood trace
#'   \item Estimated cluster mixing proportions
#'   \item Histogram of per-observation maximum posterior probability
#'   \item Heatmap of posterior probabilities sorted by cluster assignment
#' }
#'
#' @param x A `repclust` object returned by \code{\link{repclust}}.
#' @param ... Currently unused; included for S3 consistency.
#'
#' @return Invisibly returns `x`.
#'
#' @method plot repclust
#' @export
#'
#' @examples
#' sim <- generate_data(200, 3, 4, 0.2)
#' res <- repclust(sim$data, nclusters = 3)
#' plot(res)
plot.repclust <- function(x, ...) {
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit(par(op))

  if (any(par("pin") <= 0))
    stop("Plot window is too small for a 2×2 layout. Please enlarge it and call plot() again.", call. = FALSE)

  K <- ncol(x$z)
  n <- nrow(x$z)

  # 1. EM convergence
  plot(
    seq_along(x$ll), x$ll,
    type = "b", pch = 16, cex = 0.6,
    xlab = "EM iteration", ylab = "Observed log-likelihood",
    main = "EM Convergence"
  )

  # 2. Cluster mixing proportions
  pi_vec <- as.vector(x$pi)
  obs_counts <- tabulate(x$class, nbins = K)
  bp <- barplot(
    pi_vec,
    names.arg = paste0("C", seq_len(K)),
    xlab = "Cluster", ylab = "Estimated mixing proportion",
    main = "Cluster Proportions",
    ylim = c(0, max(pi_vec) * 1.15)
  )
  text(bp, pi_vec + max(pi_vec) * 0.04,
       labels = paste0("n=", obs_counts), cex = 0.8)

  # 3. Assignment ambiguity (max posterior per observation)
  max_post <- apply(x$z, 1, max)
  hist(
    max_post,
    breaks = 20,
    xlab = "Max posterior probability",
    main = "Assignment Ambiguity",
    col = "grey80", border = "white",
    xlim = c(0, 1)
  )
  if (K > 1) {
    abline(v = 1 / K, lty = 2, col = "#0072B2")
    text(1 / K, par("usr")[4], "1/K", adj = c(-0.15, 1), cex = 0.8, col = "#0072B2")
  }

  # 4. Posterior probability heatmap (rows sorted by hard assignment)
  z_sorted <- x$z[order(x$class), , drop = FALSE]
  image(
    seq_len(K), seq_len(n), t(z_sorted),
    xlab = "Cluster", ylab = "Observation (sorted by assignment)",
    main = "Posterior Probabilities",
    axes = FALSE,
    col = grey.colors(100, start = 1, end = 0)
  )
  axis(1, at = seq_len(K), labels = paste0("C", seq_len(K)))
  box()

  invisible(x)
}
