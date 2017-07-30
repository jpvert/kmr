#' Plot method for a "\code{cv.kmr}" object
#' 
#' Similar to other \code{plot} methods, this function plots pred perf over regularization parameter for a "\code{cv.kmr}" model object.
#' @param x Fitted "\code{cv.kmr}" model object
#' 
#' @export
plot.cv.kmr <- function(x) {
  
  cvobj = x
  matplot(log(cvobj$lambda), t(cvobj$meanCV), type = "l", lwd = 2, xlab = "log(lambda)", ylab = cvobj$type.measure)
  points(log(cvobj$bestlambda), cvobj$bestCV, lwd = 2)
  invisible()
}
