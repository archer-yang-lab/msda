source("utilities.R")
source("predict.mlda.R")
source("cv.mlda.R")
source("mlda.R")
library(Matrix)
dyn.load("mlda.so")


## generate data
n <- 100
p <- 1000
s <- 3
K <- 3

prior <- rep(1/K, K)

rho <- 0.5
sigma <- matrix(rho, p, p)
diag(sigma) <- 1
sigma.eigen <- eigen(sigma)
sigma.sqrt <- sigma.eigen$vectors %*% diag(sqrt(sigma.eigen$values)) %*% t(sigma.eigen$vectors)

beta <- matrix(0, p, K)
for (i in 1:K) {
    beta[1:s, i] <- i
}
mu <- sigma %*% beta

set.seed(1)
y <- runif(n)
prior.cum <- cumsum(prior)
for (i in 1:n) {
    y[i] <- sum(y[i] < prior.cum)
}

set.seed(1)
x <- matrix(rnorm(n * p), n, p)
x <- x %*% sigma.sqrt
for (i in 1:K) {
    x[y == i, ] <- sweep(x[y == i, ], 2, mu[, i], "+")
}
# mlda.prep(x,y)
pf <- rep(1, p)



    # data setup
    x <- as.matrix(x)
    y <- drop(y)
    nclass <- as.integer(length(unique(y)))
	prior<-rep(0,nclass)
    for (k in 1:nclass) {
        prior[k] <- mean(y == k)
    }
    nvars <- as.integer(ncol(x))
    nobs <- nrow(x)
    nres <- length(y)
    if (nres != nobs) 
        stop("x and y have different number of observations")
    # prepare sigma and delta
    mu <- matrix(0, nvars, nclass)
    sigma <- matrix(0, nvars, nvars)
    for (i in 1:nclass) {
        mu[, i] <- apply(x[y == i, ], 2, mean)
        sigma <- sigma + (sum(y == i) - 1) * cov(x[y == i, ])
    }
    sigma <- sigma/(n - nclass)
    delta <- sweep(mu[, -1], 1, mu[, 1], "-")
    delta <- t(delta)
    outlist <- list(sigma = sigma, delta = delta, mu = mu, prior = prior)
    outlist

## fit the model for a 100-lambda sequence
system.time(obj <- mlda(x = x, y = y, eps = 1e-04, pf = pf, verbose=FALSE, sml = 1e-4))

## KKT check the result
kktchk(obj, pf, thr = 1e-06)
# 
# try cross validation
# cv.mlda(x, y) 
