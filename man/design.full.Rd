\name{design.full}
\alias{design.full}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Setting-up design matrices for a full polytomous model with two biallelic markers.
}
\description{
This function sets up two identical lists of three design matrices, one for each linear predictor of the logit of the three outcome levels defined by the combination of two dichotomous traits against the reference level (0,0) under a full model with main effects and product terms for two biallelic markers. 
}
\usage{
design.full(x, par.constrained, constraints)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A 2-column matrix of numeric values representing genotypes of biallelic markers, with one column per marker and one row per subject. The two homozygous genotypes must 
be coded 0 and 1, and the heterozygous genotype value depends on the genetic model: 0 (recessive), 1/2 (allelic) or 1 (dominant).
}
  \item{par.constrained}{ Optional matrix of dimensions (\code{n.levels}-1) x nc specifying the parameter (1, 2 or 3) in the linear predictor for each level involved in the nc constraints BETWEEN the logistic models for different levels of the response variable, one constraint per column. 
  This functionality is not yet implemented.
}
  \item{constraints}{
Optional matrix of dimensions (\code{n.levels}-1) x nc specifying the nc linear constraints BETWEEN the logistic models for different levels of the response variable, involving the parameters specified in \code{par.constrained}, one constraint per column.  A 0 means that the corresponding parameter is not involved in the constraint.
  This functionality is not yet implemented.
}
}
\details{
Let \eqn{Y_{1}}{Y[1]} and \eqn{Y_{2}}{Y[2]} be binary variables coding the presence (1) or absence (0) of the two traits (e.g. and endophenotype and a disease trait, respectively).
The linear predictors (without intercept) of the logistic functions between outcome levels and the reference level \eqn{Y_{1} = 0}{Y[1] = 0} and \eqn{Y_{2} = 0}{Y[2] = 0} for the full model are as follows:

    \eqn{Y_{1} = 1 , Y_{2} = 0 : \beta_{11} X_{1} + \beta_{12} X_{2} + \beta_{13} X_{1} X_{2}}{Y[1] = 1 , Y[2] = 0 : \beta[11] X[1] + \beta[12] X[2] + \beta[13] X[1] X[2]}
	
    \eqn{Y_{1} = 0 , Y_{2} = 1 : \beta_{21} X_{1} + \beta_{22} X_{2} + \beta_{23} X_{1} X_{2}}{Y[1] = 0 , Y[2] = 1 : \beta[21] X[1] + \beta[22] X[2] + \beta[23] X[1] X[2]}
	
    \eqn{Y_{1} = 1 , Y_{2} = 1 : \beta_{31} X_{1} + \beta_{32} X_{2} + \beta_{33} X_{1} X_{2}}{Y[1] = 1 , Y[2] = 1 : \beta[31] X[1] + \beta[32] X[2] + \beta[33] X[1] X[2]}
	
The design matrices for the above model are constructed by this function.
}
\value{
\item{x.e }{List of 3 design matrices containing all terms forming the full model}
\item{x.loc.e }{list of character strings containing the indices of the variables in \code{x} involved in each term of the full model}
\item{x.l }{identical to \code{x.e}}
\item{x.loc.l }{identical to \code{x.l}}
}
\author{
Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>
}

\seealso{
\code{\link{fat2Lpoly}, \link{design.endo2disease}}
}

