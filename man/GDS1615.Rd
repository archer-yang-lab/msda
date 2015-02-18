\name{GDS1615}
\alias{GDS1615}
\alias{x}
\alias{y}
\docType{data}
\title{GDS1615 data introduced in Burczynski et al. (2012).}
\usage{data(GDS1615)}
\description{
The  dataset is a subset of the dataset available on Gene Expression Omnibus with the accession number GDS1615. The original dataset contains 22283 gene expression levels and the disease states of the observed subjects. In Mai, Yang and Zou, the dimension of the original dataset was first reduced to 127 by F-test screening.
}
\value{
  This data frame contains the following:
    \item{x}{gene expression levels.}
    \item{y}{Disease state that is coded as 1,2,3. 1: normal; 2: ulcerative colitis; 3: Crohn's disease.}
}
\references{
M. E. Burczynski, R. L Peterson, N. C. Twine, K. A. Zuberek, B. J. Brodeur, L. Casciotti, V. Maganti, P. S. Reddy, A. Strahs, F. Immermann, W. Spinelli, U. Schwertschlag, A. M. Slager, M. M. Cotreau, and A. J. Dorner.  (2012), "Molecular classification of crohn's disease and ulcerative colitis patients using transcriptional profiles in peripheral blood mononuclear cells". \emph{Journal of Molecular Diagnostics}, 8:51--61.

Mai, Q.*, Yang, Y.*, and Zou, H. (2014), "Multiclass Sparse Discriminant Analysis." Submitted to \emph{Journal of the American Statistical Association}. (* co-first author)

URL: \url{https://github.com/emeryyi/msda}\cr

}

\examples{
data(GDS1615)
}
\keyword{datasets}
