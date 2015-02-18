\name{plot.msda}
\alias{plot.msda}
\title{Plot coefficients from a "msda" object}
\description{
Produces a coefficient profile plot of the coefficient paths for a
fitted \code{\link{msda}} object. 
}
\usage{
\method{plot}{msda}(x, xvar = c("norm", "lambda"), ...)
}
\arguments{
		\item{x}{fitted \code{\link{msda}} model}
		\item{xvar}{the variable on the X-axis. The option \code{"norm"} plots the coefficients against the
		L1-norm of the coefficients, and the option \code{"lambda"} plots the coefficient against the log-lambda
		sequence.}
		\item{\dots}{other graphical parameters to plot}
}
\details{
A coefficient profile plot is produced.
}
\references{

Mai, Q.*, Yang, Y.*, and Zou, H. (2014), "Multiclass Sparse Discriminant Analysis." Submitted to \emph{Journal of the American Statistical Association}. (* co-first author)

URL: \url{https://github.com/emeryyi/msda}\cr
}
\author{
Qing Mai <mai@stat.fsu.edu>, Yi Yang <yiyang@umn.edu>,  Hui Zou <hzou@stat.umn.edu>\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}

\examples{
data(GDS1615)
x<-GDS1615$x
y<-GDS1615$y
obj <- msda(x = x, y = y)
plot(obj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{classification}
