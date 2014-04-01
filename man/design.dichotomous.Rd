\name{design.dichotomous}
\alias{design.dichotomous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setting-up the design matrix for a logistic model with two biallelic markers.
}
\description{
This function sets up two identical lists, each containing a design matrix for the linear predictor of the logit of a dichotomous outcome under a full logistic model with main effects and product terms for two biallelic markers. 
}
\usage{
design.dichotomous(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A 2-column matrix of numeric values representing genotypes of biallelic markers, with one column per marker and one row per subject. The two homozygous genotypes must 
be coded 0 and 1, and the heterozygous genotype value depends on the genetic model: 0 (recessive), 1/2 (allelic) or 1 (dominant).
}
  \item{...}{ Additional arguments will be ignored, but must be allowed for compatibility with other design functions.
}
}
\details{
The linear predictors (without intercept) of the logistic function for \eqn{Y = 1}{Y = 1}  against the reference level \eqn{Y = 0}{Y = 0} has the form:

    \eqn{\eta_{11} X_{1} + \eta_{12} X_{2} + \eta_{13} X_{1} X_{2}}{Y[1] = 1 , Y[2] = 0 : \eta[11] X[1] + \eta[12] X[2] + \eta[13] X[1] X[2]}

The design matrix for the above model is constructed by this function.
}
\value{
\item{x.e }{List containing the single design matrix with all terms forming the full 2-locus logistic model}
\item{x.loc.e }{list of character strings containing the indices of the variables in \code{x} involved in each term of the logistic model}
\item{x.l }{identical to \code{x.e}}
\item{x.loc.l }{identical to \code{x.l}}
}
\author{
Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>
}

\seealso{
\code{\link{fat2Lpoly}, \link{design.full}}
}

