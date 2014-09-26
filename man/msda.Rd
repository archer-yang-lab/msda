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
(x, y, nlambda = 100, lambda.factor = ifelse((nobs - nclass) <= nvars, 0.01, 1e-04), lambda = NULL, dfmax = nobs, pmax = min(dfmax * 1.2, nvars), pf = rep(1, nvars), eps = 1e-04, maxit = 1e+06, sml = 1e-04, verbose = FALSE, perturb = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{nlambda}{
%%     ~~Describe \code{nlambda} here~~
}
  \item{lambda.factor}{
%%     ~~Describe \code{lambda.factor} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{dfmax}{
%%     ~~Describe \code{dfmax} here~~
}
  \item{pmax}{
%%     ~~Describe \code{pmax} here~~
}
  \item{pf}{
%%     ~~Describe \code{pf} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{maxit}{
%%     ~~Describe \code{maxit} here~~
}
  \item{sml}{
%%     ~~Describe \code{sml} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{perturb}{
%%     ~~Describe \code{perturb} here~~
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
function (x, y, nlambda = 100, lambda.factor = ifelse((nobs - 
    nclass) <= nvars, 0.01, 1e-04), lambda = NULL, dfmax = nobs, 
    pmax = min(dfmax * 1.2, nvars), pf = rep(1, nvars), eps = 1e-04, 
    maxit = 1e+06, sml = 1e-04, verbose = FALSE, perturb = NULL) 
{
    this.call <- match.call()
    tmp <- msda.prep(x, y)
    sigma_old <- as.matrix(tmp$sigma)
    sigma <- sigma_old
    if (!is.null(perturb)) 
        diag(sigma) <- diag(sigma) + perturb
    delta <- as.matrix(tmp$delta)
    mu <- as.matrix(tmp$mu)
    prior <- tmp$prior
    nobs <- as.integer(dim(x)[1])
    nvars <- as.integer(dim(x)[2])
    nclass <- as.integer(length(unique(y)))
    vnames <- colnames(x)
    if (is.null(vnames)) 
        vnames <- paste("V", seq(nvars), sep = "")
    nk <- as.integer(dim(delta)[1])
    if (length(pf) != nvars) 
        stop("The size of penalty factor must be same as the number of input variables")
    maxit <- as.integer(maxit)
    verbose <- as.integer(verbose)
    sml <- as.double(sml)
    pf <- as.double(pf)
    eps <- as.double(eps)
    dfmax <- as.integer(dfmax)
    pmax <- as.integer(pmax)
    nlam <- as.integer(nlambda)
    if (missing(lambda)) {
        if (lambda.factor >= 1) 
            stop("lambda.factor should be less than 1")
        flmin <- as.double(lambda.factor)
        ulam <- double(1)
    }
    else {
        flmin <- as.double(1)
        if (any(lambda < 0)) 
            stop("lambdas should be non-negative")
        ulam <- as.double(rev(sort(lambda)))
        nlam <- as.integer(length(lambda))
    }
    fit <- .Fortran("msda", nk, nvars, as.double(sigma), as.double(delta), 
        pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, sml, 
        verbose, nalam = integer(1), theta = double(pmax * nk * 
            nlam), itheta = integer(pmax), ntheta = integer(nlam), 
        alam = double(nlam), npass = integer(1), jerr = integer(1))
    outlist <- formatoutput(fit, maxit, pmax, nvars, vnames, 
        nk)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr, 
        sigma = sigma, sigma_old = sigma_old, delta = delta, 
        mu = mu, prior = prior, call = this.call))
    if (is.null(lambda)) 
        outlist$lambda <- lamfix(outlist$lambda)
    class(outlist) <- c("msda")
    outlist
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
