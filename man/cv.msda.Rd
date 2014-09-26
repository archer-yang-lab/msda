\name{cv.msda}
\alias{cv.msda}
\title{
Cross-validation for msda
}
\description{
Does k-fold cross-validation for msda, returns a value for \code{lambda}.
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

\seealso{
\code{msda}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{classification}
