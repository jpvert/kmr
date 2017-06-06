#' Gaussian RBF kernel
#' 
#' Compute the Gaussian RBF kernel
#' 
#' @param x Input matrix of covariates
#' @param sigma Bandwidth of the Gaussian RBF kernel (default=1)
#' @param x2 (optional) A second matrix of covariates
#' @return If \code{x2} is provided, returns the kernel matrix between samples
#'   in x and samples in x2. If \code{x2} is not provided, returns the kernel
#'   Gram matrix of x versus itself.

gausskernel <- function(x, sigma=1, x2) {
  if (missing(x2)) {
    x2 = x
  }
  
  n = nrow(x)
  m = nrow(x2)
  i = rep(seq_len(n), times=m)
  j = rep(seq_len(m), each=n)
  u = x[i, , drop = FALSE] - x2[j, , drop = FALSE]
  kernel <- matrix(exp(- rowSums(u ^ 2) / (2*sigma^2)), n, m)
  kernel
}