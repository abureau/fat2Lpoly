\name{design.endo2disease}
\alias{design.endo2disease}
\title{
Setting-up design matrices for the endophenotype-to-disease model.
}
\description{
This function sets up two lists of three design matrices, one for each linear predictor of the logit of the three outcome levels defined by the combination of two dichotomous traits against the reference level (0,0) under the endophenotype-to-disease model of Bureau et al (2014). The design matrices in the first list contain all terms forming the model, and those in the second list contain all main effect and product terms appearing in the model.
}
\usage{
design.endo2disease(x, par.constrained, constraints)
}
\arguments{
  \item{x}{ A matrix of dimensions 2 x n of numeric values representing genotypes of biallelic markers, with one column per marker and one row per subject. The two homozygous genotypes must 
be coded 0 and 1, and the heterozygous genotype value depends on the genetic model: 0 (recessive), 1/2 (allelic) or 1 (dominant).
}
  \item{par.constrained}{ Optional matrix of dimensions (\code{n.levels}-1) x nc specifying the parameter (1 or 2) in the linear predictor for each level involved in the nc constraints BETWEEN the logistic models for different levels of the response variable, one constraint per column. 
  This functionality is not yet implemented.
}
  \item{constraints}{
Optional matrix of dimensions (\code{n.levels}-1) x nc specifying the nc linear constraints BETWEEN the logistic models for different levels of the response variable, involving the parameters specified in \code{par.constrained}, one constraint per column.  A 0 means that the corresponding parameter is not involved in the constraint.
  This functionality is not yet implemented.
}
}
\details{
Let \eqn{Y_{1}}{Y[1]} and \eqn{Y_{2}}{Y[2]} be binary variables coding the presence (1) or absence (0) of the endophenotype and the disease trait, respectively.
The linear predictors (without intercept) of the logistic functions between outcome levels and the reference level \eqn{Y_{1} = 0}{Y[1] = 0} and \eqn{Y_{2} = 0}{Y[2] = 0} specified by the endophenotype-to-disease model are as follows:

    \eqn{Y_{1} = 1 , Y_{2} = 0 : \beta_{11} X_{1} + \beta_{e} X_{1} (1 - X_{2})}{Y[1] = 1 , Y[2] = 0 : \beta[11] X[1] + \beta[e] X[1] (1 - X[2])}
	
    \eqn{Y_{1} = 0 , Y_{2} = 1 : \beta_{21} X_{1}}{Y[1] = 0 , Y[2] = 1 : \beta[21] X[1]}
	
	\eqn{Y_{1} = 1 , Y_{2} = 1 : \beta_{31} X_{1} + \beta_{33} X_{1} X_{2}}{Y[1] = 1 , Y[2] = 1 : \beta[31] X[1] + \beta[33] X[1] X[2]}
   
The design matrices for the above model are constructed by this function.
}
\value{

\item{x.e }{List of 3 design matrices containing all terms forming the endophenotype-to-disease model}
\item{x.loc.e }{list of character strings containing the indices of the variables in \code{x} involved in each term of the endophenotype-to-disease model}
\item{x.l }{List of 3 design matrices containing the terms \eqn{X_{1}}{X[1]}, \eqn{X_{2}}{X[2]} and \eqn{X_{1} X_{2}}{X[1] X[2]} appearing in the endophenotype-to-disease model.}
\item{x.loc.l }{list of character strings containing the indices of the variables in \code{x} involved in each term in the list \code{x.l}}
}
\references{
Bureau A., Croteau J., Chagnon, Y.C., Roy, M.-A. and Maziade, M. Extension of the Generalized Disequilibrium Test to polytomous phenotypes and two locus models. Frontiers in Genetics, 5: Article 258. }
\author{
Alexandre Bureau <alexandre.bureau@msp.ulaval.ca>
}
\seealso{
\code{\link{fat2Lpoly}, \link{design.full}}
}
