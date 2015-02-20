\name{msda}
\alias{msda}
\title{
Fits a regularization path for Multi-Class Sparse Discriminant Analysis}
\description{
Fits a regularization path for Multi-Class Sparse Discriminant Analysis at a sequence of regularization parameters lambda.
}
\usage{
msda(x, y, nlambda = 100, 
lambda.factor = ifelse((nobs - nclass) <= nvars, 0.2, 0.001), 
lambda = NULL, dfmax = nobs, pmax = min(dfmax * 2 + 20, nvars), 
pf = rep(1, nvars), eps = 1e-04, maxit = 1e+06, sml = 1e-06, 
verbose = FALSE, perturb = NULL)
}
\arguments{
  \item{x}{matrix of predictors, of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}
  \item{y}{response variable. This argument should be a factor for classification.
}
  \item{nlambda}{the number of \code{lambda} values - default is 100.}
  \item{lambda.factor}{
	The factor for getting the minimal lambda in \code{lambda} sequence, where \code{min(lambda)} = \code{lambda.factor} * \code{max(lambda)}.  \code{max(lambda)} is the smallest value of \code{lambda} for which all coefficients are zero. The default depends on the relationship between \eqn{N} (the number of rows in the matrix of predictors) and \eqn{p}{p} (the number of predictors). If \eqn{N > p}, the default is \code{0.0001},
	close to zero.  If \eqn{N<p}, the default is \code{0.2}.
	A very small value of \code{lambda.factor} will lead to a saturated fit. It takes no effect if there is user-defined \code{lambda} sequence.}
  \item{lambda}{
	a user supplied \code{lambda} sequence. Typically, by leaving this option unspecified users can have 
	the program compute its own \code{lambda} sequence based on
	\code{nlambda} and \code{lambda.factor}. Supplying a value of
	\code{lambda} overrides this. It is better to supply
	a decreasing sequence of \code{lambda} values than a single (small) value, if not, the program will sort user-defined \code{lambda} sequence in decreasing order automatically.}
  \item{dfmax}{
	limit the maximum number of variables in the
	model. Useful for very large \eqn{p}, if a partial path is desired. Default is \eqn{n}.}
  \item{pmax}{
	limit the maximum number of variables ever to be nonzero. For example once \eqn{\beta} enters the model, no matter how many times it exits or re-enters model through the path, it will be counted only once. Default is \code{min(dfmax*1.2,p)}.}
  \item{pf}{
	L1 penalty factor of length \eqn{p}{p}. Separate L1 penalty weights can be applied to each coefficient of \eqn{\theta}{theta} to allow
	differential L1 shrinkage. Can be 0 for some variables, which implies
	no L1 shrinkage, and results in that variable always being included in the
	model. Default is 1 for all variables (and implicitly infinity for
	variables listed in \code{exclude}).}
  \item{eps}{
	convergence threshold for coordinate descent. Each inner
	coordinate descent loop continues until the relative change in any
	coefficient. Defaults value is \code{1e-8}.}
  \item{maxit}{
maximum number of outer-loop iterations allowed at fixed lambda value. Default is 1e6. If models do not converge, consider increasing \code{maxit}.}
  \item{sml}{
%%     ~~Describe \code{sml} here~~
}
  \item{verbose}{
whether to print out computation progress. The default is \code{FALSE}. }
  \item{perturb}{a scalar number. If it is specified, the number will be added to each diagonal element of the sigma matrix as perturbation. The default is \code{NULL}.
}
}
\details{
Note that for computing speed reason, if models are not converging or running slow, consider increasing \code{eps} and \code{sml}, or decreasing
\code{nlambda}, or increasing \code{lambda.factor} before increasing
\code{maxit}. Users can also reduce \code{dfmax} to limit the maximum number of variables in the model.
}
\value{
An object with S3 class \code{\link{msda}}.
		\item{theta}{a list of length(lambda) for fitted coefficients theta, each one corresponding to one lambda value, each stored as a sparse matrix (\code{dgCMatrix} class, the standard class for sparse numeric matrices in the \code{Matrix} package.). To convert it into normal type matrix use \code{as.matrix()}.}
		\item{df}{the number of nonzero coefficients for each value of
		\code{lambda}.}
		\item{obj}{the fitted value of the objective function for each value of
		\code{lambda}.}
		\item{dim}{dimension of each coefficient matrix at each lambda.}
		\item{lambda}{the actual sequence of \code{lambda} values used.}
		\item{x}{matrix of predictors.}
		\item{y}{response variable.}
		\item{npasses}{total number of iterations (the most inner loop) summed over all lambda values}
		\item{jerr}{error flag, for warnings and errors, 0 if no error.}
		\item{sigma}{estimated sigma matrix.}
		\item{delta}{estimated delta matrix. delta[k] = mu[k]-mu[1].}
		\item{mu}{estimated mu vector.}
		\item{prior}{prior probability that y belong to class k, estimated by mean(y that belong to k).}
		\item{call}{the call that produced this object}
}
\references{

Mai, Q.*, Yang, Y.*, and Zou, H. (2014), "Multiclass Sparse Discriminant Analysis." Submitted to \emph{Journal of the American Statistical Association}. (* co-first author)

URL: \url{https://github.com/emeryyi/msda}\cr
}
\author{
Qing Mai <mai@stat.fsu.edu>, Yi Yang <yiyang@umn.edu>,  Hui Zou <hzou@stat.umn.edu>\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{cv.msda}, \code{predict.msda}
}
\examples{
data(GDS1615)
x<-GDS1615$x
y<-GDS1615$y
obj <- msda(x = x, y = y)
}
\keyword{models}
\keyword{classification}
