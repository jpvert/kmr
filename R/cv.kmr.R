#' Cross-validation for KMR
#' 
#' Does a k-fold cross-validation for \code{kmr}, and returns performance and optimal values for the regularization parameter \code{lambda}.
#' 
#' @param x \code{x} matrix as in \code{kmr}.
#' @param y Reponse matrix \code{y} as in \code{kmr}.
#' @param kx_type Kernel type for observations as in \code{kmr}.
#' @param kx_option Optional list of parameters for the observation kernel as in \code{kmr}.
#' @param kt_type Kernel type for tasks as in \code{kmr}.
#' @param kt_option Optional list of parameters for the task kernel as in \code{kmr}.
#' @param lambda Sequence of values for lambda that must be tested. Default is 10^(-5:5).
#' @param type.measure Loss to use for cross-validation. The default is \code{type.measure="ci"} the concordance index. Other options are \code{type.measure="mse"} the mean squared error.
#' @param nfolds Number of folds for cross-validation. Default is 5.
#' @param nrepeats Number of times the k-fold cross-validation is performed Default is 1.
#' @param seed A seed number for the random number generator (useful to have the same CV splits).
#' @return An object of class \code{"cv.kmr"}, which can then be used to make predictions for the different tasks on new observations as being a list containing the following useful slots:
#' \item{meanCV}{A matrix of CV performance scores of dim ntask x nlambda.}
#' \item{bestlambda}{A vector of lambdas of length ntask, each corresp to the underlying min CV score.}
#' \item{lambda}{Lambda grid to tune over.}
#' \item{type.measure}{Measure type.}
#' @export
#' 
#' 
cv.kmr <- function(x, y, kx_type=c("linear", "gaussian", "precomputed"), kx_option=list(sigma=1), kt_type=c("multitask", "empirical", "precomputed"), kt_option=list(alpha=1), lambda=10^(-5:5), type.measure = c("ci","mse"), nfolds=5, nrepeats=1, seed=9182456, mc.cores=1) {
  
  kx_type=match.arg(kx_type)
  kt_type=match.arg(kt_type)
  type.measure = match.arg(type.measure)
  N = nrow(x)
  Nt = ncol(y)
  Nl = length(lambda)
  stopifnot(Nl > 1)
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  # Set random number generator seed
  set.seed(seed)
  
  ###  Make folds
  folds <- list()
  for (i in seq(nrepeats)) {
    folds <- c(folds, split(sample(seq(N)), rep(1:nfolds, length = N)))
  }
  nexp <- length(folds)
  
  ### Iterate over the folds
  resCV <- mclapply(seq(nexp) , function(iexp) {
    message('.', appendLF = F)
    
    # training samples for the fold
    mytrain <- seq(N)[-folds[[iexp]]]
    
    # test samples for the fold
    mytest <- folds[[iexp]]
    
    # Train the model and make the prediction
    if (kx_type=="precomputed") {
      xtrain = x[mytrain,mytrain]
      xtest = x[mytest,mytrain]
    } else {
      xtrain = x[mytrain,]
      xtest = x[mytest,]
    }
    
    # Train on the training set
    m = kmr(xtrain, y[mytrain,,drop=F], kx_type, kx_option, kt_type, kt_option)
    # Predict on the test set
    ypred = predict(m, xtest, lambda=lambda)
    ytest = y[mytest,,drop=F]
    return(evalpred(ypred, ytest, type.measure))
  }, mc.cores = mc.cores)
  
  meanCV = Reduce("+", resCV) / length(resCV)
  which.lambda <- switch(type.measure,
                         "mse"= which.min,
                         "ci" = which.max)
  ilambda <- apply(meanCV,1,which.lambda)
  bestlambda <- lambda[ilambda]
  bestCV <- meanCV[cbind(seq_along(ilambda),ilambda)]
  
  ### Train model on full data
  res <- kmr(x, y, kx_type, kx_option, kt_type, kt_option)
  
  res[['meanCV']] <- meanCV
  res[['bestlambda']] <- bestlambda
  res[['bestCV']] <- bestCV
  res[['lambda']] <- lambda
  res[['type.measure']] <- type.measure
  
  class(res) <- "cv.kmr"
  return(res)
}
  
