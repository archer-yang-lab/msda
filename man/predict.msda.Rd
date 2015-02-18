\name{predict.msda}
\alias{predict.msda}
\title{
make predictions from a "msda" object.
}
\description{This functions predicts class labels from a fitted \code{\link{msda}} object.

}
\usage{
\method{predict}{msda}(object, newx, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{fitted \code{\link{msda}} model object.}
\item{newx}{matrix of new values for \code{x} at which predictions are
to be made. NOTE: \code{newx} must be a matrix, \code{predict} function does not accept a vector or other formats of \code{newx}.}
\item{\dots}{Not used. Other arguments to predict. }
}

\value{
predicted class label(s) at the entire sequence of the penalty parameter \code{lambda} used to create the model.
}
\references{

Mai, Q.*, Yang, Y.*, and Zou, H. (2014), "Multiclass Sparse Discriminant Analysis." Submitted to \emph{Journal of the American Statistical Association}. (* co-first author)

URL: \url{https://github.com/emeryyi/msda}\cr
}
\author{
Qing Mai <mai@stat.fsu.edu>, Yi Yang <yiyang@umn.edu>,  Hui Zou <hzou@stat.umn.edu>\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}

\seealso{
\code{\link{msda}}}
\examples{
data(GDS1615)
x<-GDS1615$x
y<-GDS1615$y
obj <- msda(x = x, y = y)
pred<-predict(obj,x)
}
\keyword{models}
\keyword{classification}
