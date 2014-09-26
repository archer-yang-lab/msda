\name{}
\alias{}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
(x, y, nfolds = 10, lambda.opt = "min", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{nfolds}{
%%     ~~Describe \code{nfolds} here~~
}
  \item{lambda.opt}{
%%     ~~Describe \code{lambda.opt} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, nfolds = 10, lambda.opt = "min", ...) 
{
    y <- drop(y)
    n <- nrow(x)
    p <- ncol(x)
    tmp <- msda(x, y)
    prior <- tmp$prior
    lambda <- tmp$lambda
    nlambda <- length(lambda)
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
    }
    else {
        bestlambda <- max(lambda[which(cv == min(cv))])
    }
    obj <- list(lambda = lambda, cv = cv, cv.error = cv.error, 
        bestlambda = bestlambda)
    class(obj) <- "cv.msda"
    obj
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
