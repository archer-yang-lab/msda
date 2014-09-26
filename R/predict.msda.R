# predict routine
predict.msda <- function(x, mu, theta, prior) {
    mubar <- sweep(mu[, -1], 1, mu[, 1], "+")/2
    n <- nrow(x)
    nclass <- length(prior)
    nlambda <- dim(theta)[3]
    score <- array(0, dim = c(n, nclass, nlambda))
    for (k in 2:nclass) {
        score[, k, ] <- sweep(x, 1, mubar[, k - 1], "-") %*% theta[k - 
            1, , ] + log(prior[k]) - log(prior[1])
    }
    pred <- apply(score, c(1, 3), which.max)
    pred
}
