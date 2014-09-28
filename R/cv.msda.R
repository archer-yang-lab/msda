cv.msda <- function(x, y, nfolds = 5, foldid, lambda = NULL, lambda.opt = "min", ...) {
    y <- drop(y)
    n <- NROW(x)
    p <- NCOL(x)
    K<-length(unique(y))
    prior<-rep(0,K)
    for(i in 1:K){prior[i]<-mean(y==i)}
    ### Fit the model once to get dimensions etc of output
    tmp <- msda(x, y, lambda = lambda, ...)
    lambda <- tmp$lambda
    nlambda <- length(lambda)
    ### Now fit the nfold models and store them
    if (missing(foldid)) foldid <- sample(rep(seq(nfolds), length = n)) else 			nfolds <- max(foldid)
    if (nfolds < 3) 
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    if (nfolds > n) 
        stop("The number of folds should be smaller than the sample size.")

    nlams=double(nfolds)	
    residmat=matrix(NA,nfolds,length(lambda))
    good <- matrix(0, nfolds, length(lambda))
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- msda(x[!which, , drop = FALSE], y[!which], lambda = lambda)
      	preds <- predict(fitobj,x[which,,drop=FALSE])
       	nlami <- length(fitobj$lambda)
        residmat[i,seq(nlami)] <- colMeans(y[which] != preds)
        good[i, seq(nlami)] <- 1
    }
    real_N <- apply(good, 2, sum)
    # residmat[is.na(residmat)] <- 0.5
    cv <- apply(residmat, 2, mean, na.rm = TRUE)
    cv.error <- sqrt(apply(scale(residmat, cv, FALSE)^2, 2, mean, na.rm = TRUE)/(real_N - 1))
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

