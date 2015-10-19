\name{fat2Lpoly}
\alias{fat2Lpoly}

\title{
Two-locus Family-based Association Test with Polytomous Outcome
}
\description{
Performs family-based association tests with a polytomous outcome under 2-locus and 1-locus models as described in reference [1]. Various functions \code{design.constraint} to create design matrices are provided in this package. When SNP pairs are specified, the tested SNP is the second one of each pair, while the first one is considered the conditioning SNP. The function may also perform one-locus tests if individual SNPs are specified instead of SNP pairs.
}
\usage{
fat2Lpoly(pedfilenames, datfilenames, freq.data, ibdfilenames = NULL, 
          snp.names.mat, ibd.loci = NULL, joint.tests = NULL, 
		  contingency.file = FALSE, design.constraint,
		  par.constrained, constraints, pairweights=calcule.poids.alphafixe,
		  lc = NULL, alpha = NULL)
}
\arguments{
  \item{pedfilenames}{
vector of 1 or 2 (the number of loci involved in the \code{design} function) character strings giving the names of the pedigree files in Merlin format (see Merlin website [2]). Put the full path of the files if they are not in the current working directory. If the phenotype is polytomous with 4 levels created by all combinations of two dichotomous phenotypic variables \eqn{Y_{1}}{Y[1]} and \eqn{Y_{2}}{Y[2]}, then the sixth and seventh columns of each file are respectively for \eqn{Y_{1}}{Y[1]} (e.g. the endophenotype) and \eqn{Y_{2}}{Y[2]} (e.g. the disease phenotype).
}
  \item{datfilenames}{
vector of 1 or 2 (the number of loci involved in the \code{design} function) character strings giving the names of the Merlin format data files corresponding to the pedigree files.
}
  \item{freq.data}{
Either 
(1) a vector of 1 or 2 (the number of loci involved in the \code{design} function) character strings giving the names of the allele frequency files corresponding to the pedigree files. These files must be in Merlin Classic format.
or
(2) a list of length 1 or 2 (the number of loci involved in the \code{design} function), each element of which is a numeric vector of length 'number of SNPs in \code{datfilenames}' and specifies each SNP's minor allele.
}
  \item{ibdfilenames}{
vector of 1 or 2 (the number of loci involved in the \code{design} function) character strings giving the names of the Merlin format ibd files corresponding to the pedigree files. If NULL (the default), then we use the kinship coefficients multiplied by two, instead of the expectation of the IBD probabilities, in the computation of the score statistics. The kinship coefficients are obtained using the function \code{\link{kinship}} from the package \code{kinship2}.
}
  \item{snp.names.mat}{
matrix of one or two columns giving the names of the SNPs (if one column matrix) or pairs of SNPs (if two columns matrix) to be analyzed.  These SNPs represent all or part of the SNPs in the data files \code{datfilenames}.
}
  \item{ibd.loci}{
matrix of the same dimensions as \code{snp.names.mat}, giving the respective names of the markers (used to obtain the IBD results) closest to the corresponding SNPs. The marker names must be written exactly the same as in the ibd files \code{ibdfilenames} for extraction of IBD data. If the IBD data are specified by genetic positions instead of marker names, then for each SNP, specify the genetic position where IBD was inferred which is closest to the corresponding SNP.  If NULL (the default), then we use the kinship coefficients multiplied by two, instead of the expectation of the IBD probabilities, in the computation of the score statistics. The kinship coefficients are obtained using the function \code{\link{kinship}} from the package \code{kinship2}.
}
  \item{joint.tests}{
list of vectors of numbers between 1 and the total number of parameters in the \code{design} function. Each vector gives parameter indices to test the corresponding parameters jointly.  The default is not to perform any joint test.
}
  \item{contingency.file}{
if 'TRUE' (default is 'FALSE'), then a file called descriptive_statistics'date_and_time'.txt is created and contingency tables with the numbers of subjects per level are progressively added to this file.
}
  \item{design.constraint}{
function building the design matrices WITHIN each category, for constraints specific to each category.  It also returns the design matrices comprising only the loci main effects that are used for computing the covariances. An attribute \code{n.levels} must be added within the function, to the object it returns.
}
  \item{par.constrained}{ Optional matrix of dimensions (\code{n.levels}-1) x nc specifying the parameter in the linear predictor for each level involved in the nc constraints BETWEEN the logistic models for different levels of the response variable, one constraint per column. 
  This functionality is not yet implemented.
}
  \item{constraints}{
Optional matrix of dimensions (\code{n.levels}-1) x nc specifying the nc linear constraints BETWEEN the logistic models for different levels of the response variable, involving the parameters specified in \code{par.constrained}, one constraint per column.  A 0 means that the corresponding parameter is not involved in the constraint.
  This functionality is not yet implemented.
}
  \item{pairweights}{
  	function calculating the weights of the observation pair differences when conditioning on the first SNP in the test of the second SNP in a SNP pair. Default is calcule.poids.alphafixe, implementing the weighting function of equation (6) of reference [1]. An alternative is calcule.poids.Chen, implementing the weighting function of equation (7) of reference [1].
}
  \item{lc}{
  numerical identifier of the SNP (locus) on which to condition when testing model terms. Defaults to NULL, or no conditioning.
}
  \item{alpha}{
  vector of length \code{n.levels}-1 of the coefficients of the polytomous logistic model of association beween the phenotype and the conditionning SNP. Defaults to NULL. If \code{alpha = NULL} and \code{lc} is not NULL, an \code{alpha} is obtained by logistic regression (multinomial logistic regression if \code{n.levels>2}) of the phenotype on the genotype at locus \code{lc}.
}
}
\details{
All subjects included in the pedigree files must also be found in the IBD files.

All fields in the pedigree files must be numeric.  No letters allowed, even for family and subject ID's.

Families whose genotyped subjects are all in the same category (phenotype combination), are uninformative and will be excluded.

Conditioning on the first SNP in a SNP pair is implemented by weighting the observation pair differences according to a model of the polytomous outcome as a function of the first SNP genotypes. The function converting the coefficients of this regression model into weights is specified by the argument \code{pairweights}. The default function \code{calcule.poids.alphafixe} provided satisfactory power in simulations described in reference [1].

File "descriptive_statistics'date_and_time'.txt" (will be created if contingency.file='TRUE'): For each tested SNP, it shows contingency tables of the subjects in the 2 or 4 different categories, first for all families together and then for each individual family.

If one or both of the arguments ibd.loci and ibdfilenames are left unspecified (or NULL, their default), then we use the kinship coefficients multiplied by two, instead of the expectation of the IBD probabilities, in the computation of the score statistics. The kinship coefficients are obtained using the function \code{\link{kinship}} from the package \code{kinship2}.
}
\value{
returns a list of five objects:
  \item{scores.covs.all.SNPs}{
list of length 'nrow(\code{snp.names.mat})', each element of which contains the estimates of the scores and covariances of all the families.
}
 \item{p.values.scores}{
data frame of p-values for all the SNPs or SNP pairs in \code{snp.names.mat}, for the global test (all parameters tested jointly), the individual tests and other joint tests specified by the argument \code{joint.tests}.  The p-values are obtained from scores summed over all families.  These scores of individual tests are also included in this data frame.
}
 \item{MA.table}{
data frame giving the minor allele numbers of all the SNPs contained in the allele frequency files.
}
 \item{y1}{
affection name extracted from first line of the data file(s)
}
 \item{y2}{
affection name extracted from second line of the data file(s)
}
}
\references{
1. Bureau A., Croteau J., Chagnon, Y.C., Roy, M.-A. and Maziade, M. Extension of the Generalized Disequilibrium Test to polytomous phenotypes and two locus models. Frontiers in Genetics, 5: Article 258. 
2. http://www.sph.umich.edu/csg/abecasis/Merlin/tour/input_files.html
}
\author{
Alexandre Bureau and Jordie Croteau
}
\seealso{
\code{\link{fat2Lpoly.withinR}}
}
\examples{
path.data=paste(.libPaths()[which(unlist(lapply(.libPaths(),
function(x) length(grep("fat2Lpoly",dir(x)))))>0)],"/fat2Lpoly/extdata/",sep="")
if(length(path.data)>1) path.data=path.data[length(path.data)]

snps.anal=c("snp3.loc2","snp4.loc2")
microsat.names.loc2=c("2_3_mrk:","2_4_mrk:")

############ design.endo2disease with conditioning on locus 1 ################
\dontrun{
joint.tests=list(c(2,5))
snp.names.mat=cbind(rep("snp4.loc1",length(snps.anal)),snps.anal)
microsat.names.mat=cbind(rep("1_4_mrk:",length(snps.anal)),microsat.names.loc2)
test=fat2Lpoly(pedfilenames=paste(path.data,c("loc1.ped","loc2.ped"),sep=""),
               datfilenames=paste(path.data,c("loc1.dat","loc2.dat"),sep=""),
			   freq.data=paste(path.data,c("loc1.freq","loc2.freq"),sep=""),
               ibdfilenames=paste(path.data,c("loc1.ibd","loc2.ibd"),sep=""),
		       snp.names.mat=snp.names.mat,ibd.loci=microsat.names.mat,
		       joint.tests=joint.tests,contingency.file=TRUE,
		       design.constraint=design.endo2disease,lc=1)

test$p.values.scores
}		   
###############################################################################

################### design.endo2disease without conditioning ##################
joint.tests=list(c(2,5))
snp.names.mat=cbind(rep("snp4.loc1",length(snps.anal)),snps.anal)
microsat.names.mat=cbind(rep("1_4_mrk:",length(snps.anal)),microsat.names.loc2)
test=fat2Lpoly(pedfilenames=paste(path.data,c("loc1.ped","loc2.ped"),sep=""),
               datfilenames=paste(path.data,c("loc1.dat","loc2.dat"),sep=""),
			   freq.data=paste(path.data,c("loc1.freq","loc2.freq"),sep=""),
               ibdfilenames=paste(path.data,c("loc1.ibd","loc2.ibd"),sep=""),
		       snp.names.mat=snp.names.mat,ibd.loci=microsat.names.mat,
		       joint.tests=joint.tests,contingency.file=TRUE,
		       design.constraint=design.endo2disease)

test$p.values.scores   
###############################################################################

################# design.full with conditioning on locus 1 ##################
\dontrun{
joint.tests=list(c(2,3),c(5,6),c(8,9),c(2,3,5,6,8,9))
snp.names.mat=cbind(rep("snp4.loc1",length(snps.anal)),snps.anal)
microsat.names.mat=cbind(rep("1_4_mrk:",length(snps.anal)),microsat.names.loc2)
test=fat2Lpoly(pedfilenames=paste(path.data,c("loc1.ped","loc2.ped"),sep=""),
               datfilenames=paste(path.data,c("loc1.dat","loc2.dat"),sep=""),
			   freq.data=paste(path.data,c("loc1.freq","loc2.freq"),sep=""),
               ibdfilenames=paste(path.data,c("loc1.ibd","loc2.ibd"),sep=""),
		       snp.names.mat=snp.names.mat,ibd.loci=microsat.names.mat,
		       joint.tests=joint.tests,
               design.constraint=design.full,lc=1)

test$p.values.scores
}
##############################################################################

############################# design.1locus #################################
snp.names.mat=as.matrix(snps.anal)
microsat.names.mat=as.matrix(microsat.names.loc2)
test=fat2Lpoly(pedfilenames=paste(path.data,"loc2.ped",sep=""),
               datfilenames=paste(path.data,"loc2.dat",sep=""),
               freq.data=paste(path.data,"loc2.freq",sep=""),
			   ibdfilenames=paste(path.data,"loc2.ibd",sep=""),
		       snp.names.mat=snp.names.mat,ibd.loci=microsat.names.mat,
			   design.constraint=design.1locus)

test$p.values.scores			   
##############################################################################

############# design.dichotomous with conditioning on locus 1 ##############
\dontrun{
joint.tests=list(c(2,3))
snp.names.mat=cbind(rep("snp4.loc1",length(snps.anal)),snps.anal)
microsat.names.mat=cbind(rep("1_4_mrk:",length(snps.anal)),microsat.names.loc2)
test=fat2Lpoly(pedfilenames=paste(path.data,c("loc1.ped","loc2.ped"),sep=""),
               datfilenames=paste(path.data,c("loc1.dat","loc2.dat"),sep=""),
			   freq.data=paste(path.data,c("loc1.freq","loc2.freq"),sep=""),
               ibdfilenames=paste(path.data,c("loc1.ibd","loc2.ibd"),sep=""),
		       snp.names.mat=snp.names.mat,ibd.loci=microsat.names.mat,
		       joint.tests=joint.tests,
               design.constraint=design.dichotomous,lc=1)

test$p.values.scores
}			   
##############################################################################
}




