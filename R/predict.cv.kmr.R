#' Make prediction with a "\code{cv.kmr}" object
#' 
#' Similar to other \code{predict} methods, this function predicts fitted values
#' from a fitted "\code{cv.kmr}" object.
#' 
#' @param object Fitted "\code{cv.kmr}" model object.
#' @param newx Matrix of new values for \code{x} or kernel matrix for new values
#'  crossing old values for \code{x}, at which predictions are to be made.
#' @param lambda Value(s) of the regularization parameter, a single scalar or a
#'   sequence of values, or "\code{lambda.opt}" which allows different lambdas 
#'   for each task tuned by cross-validation. Default is "\code{lambda.opt}".
#' @param ... Not used. Other arguments to predict.
#'   
#' @return A matrix of predicted values for the new samples (in rows) and all
#'   tasks (in columns), corresponding to the regularization parameter
#'   \code{lambda}. If \code{lambda} is a list of specified values, the function 
#'   returns a list of matrices, corresponding to the predictions for the different 
#'   values of \code{lambda} in the list. If \code{lambda} is "\code{lambda.opt}",
#'   a single matrix returned too.
#'   
#' @note Note that "\code{lambda.opt}" will allow different lambdas for each task tuned
#'   by cross-validation; otherwise, specified value(s) of lambda will be fixed for all tasks.
#'   
#' @export
#' 
#' @seealso \code{\link{cv.kmr}}
#' 

predict.cv.kmr <- function(object, 
                           newx, 
                           lambda = "lambda.opt", 
                           ...)
{
  if (missing(lambda)) {
    lambda <- object$bestlambda
    i.lambda <- match(lambda, object$lambda)
  }
  
  # predict along lambda
  pred <- predict.kmr(object, newx, lambda)
  
  # Return the list of prediction matrices if there are several specified lambda's,
  # or just the prediction matrix if there is a single lambda or "lambda.opt"
  if (length(pred) == 1) {
    return(pred[[1]])
  } else if (exists("i.lambda")) {
    pred <- lapply(seq_along(i.lambda), function(i) pred[[i.lambda[i]]][ ,i,drop=T])
    return(do.call("cbind", pred))
  } else {
    return(pred)
  }
}