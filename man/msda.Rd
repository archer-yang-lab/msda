\name{msda}
\alias{msda}
\title{
Fits a regularization path for Multi-class Sparse Discriminant Analysis}
\description{
Fits a regularization path for Multi-class Sparse Discriminant Analysis at a sequence of regularization parameters lambda.
}
\usage{
(x, y, nlambda = 100, 
lambda.factor = ifelse((nobs - nclass) <= nvars, 0.01, 1e-04), 
lambda = NULL, dfmax = nobs, pmax = min(dfmax * 1.2, nvars), 
pf = rep(1, nvars), eps = 1e-04, maxit = 1e+06, sml = 1e-04, 
verbose = FALSE, perturb = NULL)
}
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

Mai, Q.*, Yang, Y.*, and Zou, H. (2014), "Multiclass Sparse Discriminant Analysis." Submitted to \emph{Journal of the American Statistical Association}. (* co-first author)

URL: \url{https://github.com/emeryyi/msda}\cr
}
\author{
Qing Mai <mai@stat.fsu.edu>, Yi Yang <yiyang@umn.edu>,  Hui Zou <hzou@stat.umn.edu>\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{cv.msda}, \code{predict.msda}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{classification}
