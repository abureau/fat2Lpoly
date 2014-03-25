\name{fat2Lpoly.withinR}
\alias{fat2Lpoly.withinR}

\title{
Two-locus Family-based Association Test with Polytomous Outcome (all arguments within R)
}
\description{
Same as \code{fat2Lpoly} except that the first four arguments of \code{fat2Lpoly} are replaced by one object having the format of the objects returned by \code{read.merlin.files}.
}
\usage{
fat2Lpoly.withinR(ped.x.all, snp.names.mat, ibd.loci = NULL, contingency.file = FALSE, 
                  design.constraint, par.constrained, constraints,lc=1,alpha.vec=c(0,0,0))
}
\arguments{
  \item{ped.x.all}{
object returned by the function \code{read.merlin.files} or having the same format.
}
  \item{snp.names.mat}{
matrix of one or two columns giving the names of the SNPs (if one column matrix) or pairs of SNPs (if two columns matrix) to be analyzed.  These SNPs represent all or part of the SNPs in the data files \code{datfilenames}.
}
  \item{ibd.loci}{
matrix of the same dimensions as \code{snp.names.mat}, giving the respective names of the markers (used to obtain the IBD results) nearest to the corresponding SNPs. The marker names must be written exactly the same as in the ibd files \code{ibdfilenames} for extraction of IBD data.  If the IBD data are specified by genetic positions instead of marker names, this matrix must contain the genetic positions of the markers instead of the marker names. If NULL (the default), then we use the kinship coefficients multiplied by two, instead of the expectation of the IBD probabilities, in the computation of the score statistics. The kinship coefficients are obtained using the function \code{\link{kinship}} from the package \code{kinship2}.
}
  \item{contingency.file}{
if 'TRUE' (default is 'FALSE'), then a file called descriptive_statistics'date_and_time'.txt is created and contingency tables with the numbers of subjects per level are progressively added to this file.
}
  \item{design.constraint}{
function building the design matrices WITHIN each category, for constraints specific to each category.  It also returns the design matrices comprising only the loci main effects that are used for computing the covariances.
}
  \item{par.constrained}{ Optional matrix of dimensions (\code{n.levels}-1) x nc specifying the parameter in the linear predictor for each level involved in the nc constraints BETWEEN the logistic models for different levels of the response variable, one constraint per column. 
  This functionality is not yet implemented.
}
  \item{constraints}{
Optional matrix of dimensions (\code{n.levels}-1) x nc specifying the nc linear constraints BETWEEN the logistic models for different levels of the response variable, involving the parameters specified in \code{par.constrained}, one constraint per column.  A 0 means that the corresponding parameter is not involved in the constraint.
  This functionality is not yet implemented.
}
  \item{lc}{
  numerical identifier of the SNP (locus) on which to condition when testing model terms. Default to 1.
  }
  \item{alpha.vec}
  vector of length \code{n.levels}-1 of the coefficients of the polytomous logistic model of association beween the phenotype and the conditionning SNP. Defaults to a vector of zeros.
}
\details{
File "descriptive_statistics'date_and_time'.txt" (will be created if contingency.file='TRUE'): For each tested SNP, it shows contingency tables of the subjects in the 2 or 4 different categories, first for all families together and then for each individual family.

If the argument ibd.loci is left unspecified (or NULL, its default), then we use the kinship coefficients multiplied by two, instead of the expectation of the IBD probabilities, in the computation of the score statistics. The kinship coefficients are obtained using the function \code{\link{kinship}} from the package \code{kinship2}.
}
\value{
  \item{scores.covs.all.SNPs}{
list of length 'nrow(\code{snp.names.mat})', each element of which contains the estimates of the scores and covariances of all the families.
}
  \item{snp.names.mat}{
(same matrix as provided as argument) matrix of one or two columns giving the names of the SNPs (if one column matrix) or pairs of SNPs (if two columns matrix) to be analyzed.  These SNPs represent all or part of the SNPs in the data files \code{datfilenames}.
}
}
\references{
Bureau A., Croteau J., Chagnon, Y.C., Roy, M.-A. and Maziade, M. Extension of the Generalized Disequilibrium Test to polytomous phenotypes and two locus models.
}
\author{
Alexandre Bureau and Jordie Croteau
}
\seealso{
\code{\link{fat2Lpoly}, \link{read.merlin.files}, \link{get.scores.pvalues}}
}
\examples{
data(ped.x.all)

snp.names.mat=cbind(rep("snp4.loc1",5),c("snp2.loc2","snp4.loc2","snp5.loc2",
                                         "snp6.loc2","snp8.loc2"))	
microsat.names.mat=cbind(rep("1_4_mrk:",5),c("2_2_mrk:","2_4_mrk:","2_5_mrk:",
                                             "2_6_mrk:","2_8_mrk:"))	
fat2Lpoly.allSNPs=fat2Lpoly.withinR(ped.x.all,snp.names.mat,ibd.loci=
                                    microsat.names.mat,contingency.file=TRUE,
									design.constraint=design.endo2disease)

joint.tests=list(c(2,5))   
get.scores.pvalues(fat2Lpoly.allSNPs,joint.tests)
#    snp.cond  snp.test global_p params.joint_2-5_p param_1_score param_2_score  
# 1 snp4.loc1 snp2.loc2 5.80e-03           7.12e-01         0.449         0.333  
# 2 snp4.loc1 snp4.loc2 2.14e-07           1.24e-05         0.449         0.890 
# 3 snp4.loc1 snp5.loc2 1.14e-03           1.44e-01         0.449         0.776 
# 4 snp4.loc1 snp6.loc2 5.59e-04           3.84e-02         0.449        -0.082 
# 5 snp4.loc1 snp8.loc2 1.15e-03           1.55e-01         0.449         0.869
# param_3_score param_4_score param_5_score param_1_p param_2_p param_3_p param_4_p
#        -1.427         3.638         0.733     0.653     0.739     0.154  0.000275
#        -1.427         3.638         4.612     0.653     0.373     0.154  0.000275
#        -1.427         3.638         1.785     0.653     0.438     0.154  0.000275		
#        -1.427         3.638         2.553     0.653     0.934     0.154  0.000275
#        -1.427         3.638         1.695     0.653     0.385     0.154  0.000275				
#   param_5_p
# 1  0.464000
# 2  0.000004
# 3  0.074200
# 4  0.010700
# 5  0.090100   
}




