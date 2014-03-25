\name{fat2Lpoly.allSNPs}
\alias{fat2Lpoly.allSNPs}
\docType{data}
\title{
Example results output by the function \code{fat2Lpoly.withinR}
}
\description{
This list is an example of output from the function \code{fat2Lpoly.withinR}. It is provided to test the function \code{get.scores.pvalues} by executing the example code in the \code{get.scores.pvalues} documentation.
}
\usage{data(fat2Lpoly.allSNPs)}
\format{
A list of two objects:
  \describe{
  \item{scores.covs.all.SNPs}{
list of length 'nrow(\code{snp.names.mat})', each element of which contains the estimates of the scores and covariances of all the families.
}
  \item{snp.names.mat}{
(same matrix as provided as argument) matrix of one or two columns giving the names of the SNPs (if one column matrix) or pairs of SNPs (if two columns matrix) to be analyzed.  These SNPs represent all or part of the SNPs in the data files \code{datfilenames}.
}
}
}
\examples{
data(fat2Lpoly.allSNPs)
}
\keyword{datasets}
