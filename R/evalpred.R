#' Evaluate KMR
#' 
#' Evaluate prediction performance of KMR.
#' 
#' @param ypred Matrix or list of matrices (eg. each corresp to a lambda) of predicted reponses of cell lines x tasks.
#' @param ytest Matrix of true responses of cell lines x tasks.
#' @param type.measure Measure type of evaluation. The default is \code{type.measure="ci"} the concordance index. Other options are \code{type.measure="mse"} the mean squared error.
#' @importFrom Hmisc rcorr.cens
#' @export
#' 
#' 

evalpred <- function(ypred, ytest, type.measure = c("ci","mse")) {
  
  type.measure <- match.arg(type.measure)
  errfun <- switch(type.measure,
              "mse"= function(u) {apply(((u-ytest)^2),2,mean)},
              "ci" = function(u) {apply(rbind(u,ytest), 2, function(v) {Hmisc::rcorr.cens(v[1:nrow(u)], v[-(1:nrow(u))], outx=FALSE)[[1]]})})
  if (!is.list(ypred))
    ypred <- list(ypred)
  return(sapply(ypred, errfun))
}
