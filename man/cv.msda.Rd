\name{cv.msda}
\alias{cv.msda}
\title{
Cross-validation for msda
}
\description{
Does k-fold cross-validation for msda, returns a value for \code{lambda}.
}
\usage{
cv.msda(x, y, nfolds = 5, lambda = NULL, lambda.opt = "min", ...)
}
\arguments{
  \item{x}{
matrix of predictors, of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}
  \item{y}{
response variable. This argument should be a factor for classification.}
  \item{nfolds}{
	number of folds - default is 5. Although \code{nfolds}
	can be as large as the sample size (leave-one-out CV), it is not
	recommended for large datasets. Smallest value allowable is \code{nfolds=3}.}
  \item{lambda}{optional user-supplied lambda sequence; default is
	\code{NULL}, and \code{\link{msda}} chooses its own sequence.}
  \item{lambda.opt}{If choose \code{"min"}, 
	the smallest \code{lambda} that gives minimum
	cross validation error \code{cvm} will be returned. If choose \code{"max"}, 
		the largest \code{lambda} that gives minimum
		cross validation error \code{cvm} will be returned.
}
  \item{\dots}{
other arguments that can be passed to msda.}
}
\details{
The function runs \code{\link{msda}} \code{nfolds}+1 times; the
first to get the \code{lambda} sequence, and then the remainder to
compute the fit with each of the folds omitted. The average error and standard deviation over the folds are computed.}
\value{
an object of class \code{\link{cv.msda}} is returned, which is a
list with the ingredients of the cross-validation fit.
		\item{lambda}{the values of \code{lambda} used in the fits.}
		\item{cvm}{the mean cross-validated error - a vector of length
		\code{length(lambda)}.}
		\item{cvsd}{estimate of standard error of \code{cvm}.}
		\item{lambda.min}{the optimal value of \code{lambda} that gives minimum
		cross validation error \code{cvm}.}
		\item{lambda.1se}{the largest value of \code{lambda} such that error is
		within 1 standard error of the minimum.}
		\item{msda.fit}{a fitted \code{\link{msda}} object for the full data.}
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
obj.cv<-cv.msda(x=x[,1:10],y=y,nfolds=5,lambda.opt="max")
lambda.min<-obj.cv$lambda.min
id.min<-which(obj.cv$lambda==lambda.min)
pred<-predict(obj.cv$msda.fit,x[,1:10])[,id.min]
}
\keyword{models}
\keyword{classification}
