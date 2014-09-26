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
(x, mu, theta, prior)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{prior}{
%%     ~~Describe \code{prior} here~~
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
function (x, mu, theta, prior) 
{
    mubar <- sweep(mu[, -1], 1, mu[, 1], "+")/2
    n <- nrow(x)
    nclass <- length(prior)
    nlambda <- dim(theta)[3]
    score <- array(0, dim = c(n, nclass, nlambda))
    for (k in 2:nclass) {
        score[, k, ] <- sweep(x, 1, mubar[, k - 1], "-") \%*\% 
            theta[k - 1, , ] + log(prior[k]) - log(prior[1])
    }
    pred <- apply(score, c(1, 3), which.max)
    pred
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
