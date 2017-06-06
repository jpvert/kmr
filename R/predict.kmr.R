#' Make prediction with a kmr object
#' 
#' Similar to other \code{predict} methods, this function predicts fitted values
#' from a fitted "\code{kmr}" object
#' @param object Fitted "\code{kmr}" model object
#' @param newx Matrix of new values for \code{x} at which predictions are to be
#'   made.
#' @param lambda Value(s) of the regularization parameter, a single scalar or a sequence of values. Default is 1. 


predict.kmr <- function(object, newx, lambda=1) {
  
  # Kernel matrix between test and train sets
  K = switch(object$kx_type,
             precomputed=newx,
             linear=newx %*% object$x,
             gaussian=gausskernel(newx, object$x, object$kx_option$sigma))
  
  
}