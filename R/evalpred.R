#' Evaluate KMR
#' 
#' Evaluate prediction performance of KMR.
#' 
#' @param ypred Matrix or list of matrices (eg. each corresp to a lambda) of predicted reponses of cell lines x tasks.
#' @param ytest Matrix of true responses of cell lines x tasks.
#' @param type.measure Character indicating the measure type of evaluation. Possible options are \code{"ci"} (concordance index, default), \code{"mse"} (mean squared error), \code{"cor"} (pearson correlation).
#' @useDynLib kmr cidx
#' @export
#' 
evalpred <- function(ypred, ytest, type.measure = c("ci","mse","cor")) {
  
  type.measure <- match.arg(type.measure)
  
  errfun <- switch(type.measure,
              "ci" = 
                function(u) {
                  apply(rbind(u,ytest), 2, function(v) {
                    if (isTRUE(requireNamespace("Hmisc", quietly = TRUE))) {
                      Hmisc::rcorr.cens(v[1:nrow(u)], v[-(1:nrow(u))], outx = FALSE)[[1]]
                    } else {
                      .Call(cidx, v[1:nrow(u)], v[-(1:nrow(u))])
                    }
                  })
                },
              "mse" = 
                function(u) {
                  apply( ((u-ytest)^2), 2, mean )
                },
              "cor" = function(u) {
                diag(cor(u, ytest))
              }
  )
  
  if (!is.list(ypred))
    ypred <- list(ypred)
  
  return(sapply(ypred, errfun))
}
