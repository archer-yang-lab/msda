cv.msda <- function(x, y, nfolds = 5, lambda = NULL, lambda.opt = "min", ...) {
    y <- drop(y)
    n <- nrow(x)
    p <- ncol(x)
    K<-length(unique(y))
    prior<-rep(0,K)
    for(i in 1:K){prior[i]<-mean(y==i)}
    ### Fit the model once to get dimensions etc of output
    tmp <- msda(x, y, lambda = lambda, ...)
    lambda <- tmp$lambda
    ### Now fit the nfold models and store them
    foldid <- sample(rep(seq(nfolds), length = n))
    if (nfolds < 3) 
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    if (nfolds > n) 
        stop("The number of folds should be smaller than the sample size.")
    residmat=matrix(NA,nfolds,length(lambda))
    good <- matrix(0, nfolds, length(lambda))
    for (i in seq(nfolds)) {
        which <- foldid == i
		cat("\nFold ", i, " is running.\n")
        fit <- msda(x[!which, , drop = FALSE], y[!which], lambda = lambda, ...)
      	preds <- predict(fit,x[which,,drop=FALSE])
       	nlami <- length(fit$lambda)
        residmat[i,seq(nlami)] <- colMeans(y[which] != preds)
        good[i, seq(nlami)] <- 1
    }
    rN <- colSums(good)
    cv <- colMeans(residmat, na.rm = TRUE)
    cv.error <- sqrt(colMeans(scale(residmat, cv, FALSE)^2, na.rm = TRUE)/(rN - 1))
    if (lambda.opt == "min") {
        bestlambda <- min(lambda[which(cv == min(cv,na.rm=TRUE))])
    } else {
        bestlambda <- max(lambda[which(cv == min(cv,na.rm=TRUE))])
    }
    id.min<-which.min(cv)
    lambda.1se<-max(lambda[cv<min(cv,na.rm=TRUE)+cv.error[id.min]],na.rm=TRUE)
    cv <- na.omit(cv)
    cv.error <- na.omit(cv.error)
	min_len <- min(length(cv), length(cv.error))
    obj <- list(lambda = lambda[1:min_len], cv = cv[1:min_len], cv.error = cv.error[1:min_len], bestlambda = bestlambda,lambda.1se=lambda.1se)
    class(obj) <- "cv.msda"
    obj
} 