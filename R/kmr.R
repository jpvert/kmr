#' Kernel multitask regression
#' 
#' Trains a kernel multitask regression model
#' 
#' @param x Input matrix of covariates, of dimension \code{nobs x nvars}; each 
#'   row is an observation vector. If a precomputed kernel is used, then 
#'   \code{x} is the square Gram matrix.
#' @param y Output matrix of responses. \code{y} should be an \code{nobs x 
#'   ntasks} matrix, where each row corresponds to an observation and each 
#'   column to a task.
#' @param kx_type Kernel for observations.  \code{kx_type="linear"} is the
#'   linear kernel (default). \code{kx_type="gaussian"} is the Gaussian RBF
#'   kernel with bandwidth \code{sigma=1} by default, or any other valued is
#'   passed as an element of the \code{kx_option} list.
#'   \code{kx_type="precomputed"} assumes that the \code{x} matrix provided is
#'   the kernel Gram matrix.
#' @param kx_option An optional list of parameters for the observation kernel.
#' @param kt_type Kernel for tasks. \code{kt_type="multitask"} (default) is the 
#'   multitask kernel with parameter \eqn{0 \le \alpha \le 1}, which 
#'   interpolates between the Dirac kernel for \eqn{\alpha=1} (default) and the 
#'   constant kernel for \eqn{\alpha=0}. The parameter \eqn{alpha} can be passed
#'   to the \code{kt_option} list as a field \code{alpha}. 
#'   \code{kt_type="empirical"} takes the empirical correlation between outputs 
#'   as kernel between the tasks. \code{kt_type="precomputed"} allows to provide
#'   a precomputed kernel as a field \code{kt} in the \code{kt_type} list.
#' @param kt_option An optional list of parameters for the task kernel.
#'   
#' @return An object of class \code{"kmr"}, which can then be used to make 
#'   predictions for the different tasks on new observations.
#' @export
kmr <- function(x, y, lambda=1, kx_type=c("linear", "gaussian", "precomputed"), kx_option=list(sigma=1), kt_type=c("multitask", "empirical", "precomputed"), kt_option=list(alpha=1)) {

  this.fcall=match.call()
  kx_type=match.arg(kx_type)
  kt_type=match.arg(kt_type)
  
  # Eigenvectors and eigenvalues of the covariate kernel
  if (kx_type=="linear") {
    # In the case of the linear kernel we just make a SVD of the data matrix
    s=svd(x)
    Ux=s$u
    Dx=(s$d)^2
    rm(s)
  } else {
    # Otherwise we make an eigendecomposition of the kernel Gram matrix
    if (kx_type=="precomputed") {
      Kx = x
    } else {
      # Gaussian kernel
      Kx = exp(-as.matrix(dist(x))^2/(2*kx_option[['sigma']]^2))
    }
    s=eigen(Kx,symmetric=TRUE)
    Ux=s$vectors
    Dx=s$values
    rm(s)
  }
  
  # Eigenvectors and eigenvalues of the task kernel
  ntasks=ncol(y)
  Kt = switch(kt_type,
              multitask=kt_option[['alpha']]*diag(ntasks)+matrix(1-kt_option[['alpha']],nrow=ntasks,ncol=ntasks),
              empirical=cor(y),
              precomputed=kt_option[['kt']]
              )
  s=eigen(Kt,symmetric = TRUE)
  Ut=s$vectors
  Dt=s$values
  rm(s)
  
  # Return
  res=list(Ux=Ux,Dx=Dx,Ut=Ut,Dt=Dt,kx_type=kx_type,kx_option=kx_option,call=this.fcall)
  class(res)="kmr"
  return(res)
}