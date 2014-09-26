# cross validation
cv.msda <- function(x, y, nfolds = 10, lambda.opt = "min", ...) {
    y <- drop(y)
    n <- nrow(x)
    p <- ncol(x)
    ### Fit the model once to get dimensions etc of output
    tmp <- msda(x, y)
    prior <- tmp$prior
    lambda <- tmp$lambda
    nlambda <- length(lambda)
    ### Now fit the nfold models and store them
    foldid <- sample(rep(seq(nfolds), length = n))
    if (nfolds < 3) 
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    if (nfolds > n) 
        stop("The number of folds should be smaller than the sample size.")
    residmat <- matrix(0, nlambda, nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fit <- msda(x[!which, , drop = FALSE], y[!which], lambda = lambda)
        mu <- fit$mu
        theta <- list2array(fit$theta)
        pred <- predict.msda(x[which, ], mu, theta, prior)
        for (l in 1:nlambda) {
            residmat[l, i] <- mean(y[which] != pred[, l])
        }
    }
    residmat[is.na(residmat)] <- 0.5
    residmat <- matrix(residmat, nrow = nlambda)
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    if (lambda.opt == "min") {
        bestlambda <- max(lambda[which(cv == min(cv))])
    } else {
        bestlambda <- max(lambda[which(cv == min(cv))])
    }
    obj <- list(lambda = lambda, cv = cv, cv.error = cv.error, bestlambda = bestlambda)
    class(obj) <- "cv.msda"
    obj
} 
