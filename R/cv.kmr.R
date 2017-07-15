#' Cross-validation for KMR
#' 
#' Does a k-fold cross-validation for \code{kmr}, and returns performance and optimal values for the regularization parameter \code{lambda}.
#' 
#' @param x \code{x} matrix as in \code{kmr}
#' @param y Reponse matrix \code{y} as in \code{kmr}
#' @param kx_type Kernel type for observations as in \code{kmr}
#' @param kx_option Optional list of parameters for the observation kernel as in \code{kmr}
#' @param kt_type Kernel type for tasks as in \code{kmr}
#' @param kt_option Optional list of parameters for the task kernel as in \code{kmr}
#' @param lambda Sequence of values for lambda that must be tested. Default is 10^(-5:5)
#' @param type.measure Loss to use for cross-validation. The default is \code{type.measure="mse"} the mean squared error. \code{type.measure="ci"} measures the concordance index.
#' @param nfolds Number of folds for cross-validation. Default is 5.
#' @param nrepeats Number of times the k-fold cross-validation is performed Default is 10.
#' @param seed A seed number for the random number generator (useful to have the same CV splits)
#' @export
#' 
#' 
cv.kmr <- function(x, y, kx_type=c("linear", "gaussian", "precomputed"), kx_option=list(sigma=1), kt_type=c("multitask", "empirical", "precomputed"), kt_option=list(alpha=1), lambda=10^(1:10), type.measure = c("mse","ci"), nfolds=5, nrepeats=10, seed=9182456) {
  
  kx_type=match.arg(kx_type)
  kt_type=match.arg(kt_type)
  type.measure = match.arg(type.measure)
  N = nrow(x)
  Nt = ncol(y)
  Nl = length(lambda)
  
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
    cat('.')
    
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
    # Assess performance
    errfun = switch(type.measure,
                    "mse"= function(u) {apply(((u-ytest)^2),2,mean)},
                    "ci" = function(u) {apply(rbind(u,ytest), 2, function(v) { rcorr.cens(v[1:nrow(u)], v[-(1:nrow(u))], outx=FALSE)[[1]]})})
    
    return(sapply(ypred, errfun)) } , mc.cores=mc.cores)
  
  meanCV = Reduce("+", resCV) / length(resCV)
  matplot(t(meanCV),type="l",lwd=2)
  bestlambda=lambda[apply(meanCV,1,which.min)]
  return(meanCV,bestlambda)
}
  
