\name{fat2Lpoly}
\alias{fat2Lpoly}

\title{
Two-locus Family-based Association Test with Polytomous Outcome
}
\description{
For a list of SNP pairs, performs family-based association tests with a polytomous outcome under 2-locus and 1-locus models as described in [1]. Various functions \code{design.constraint} to create design matrices are provided in this package.  The tested SNP is the second one of each pair, while the first one is considered the conditioning SNP. The function may also perform one-locus tests if individual SNPs are specified instead of SNP pairs.
}
\usage{
fat2Lpoly(pedfilenames, datfilenames, freq.data, ibdfilenames = NULL, 
          snp.names.mat, ibd.loci = NULL, joint.tests = NULL, 
		  contingency.file = FALSE, design.constraint,
		  par.constrained, constraints)
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
}
\details{
All subjects included in the pedigree files must also be found in the IBD files.

All fields in the pedigree files must be numeric.  No letters allowed, even for family and subject ID's.

Families whose genotyped subjects are all in the same category (phenotype combination), are uninformative and will be excluded.

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
1. Bureau A., Croteau J., Chagnon, Y.C., Roy, M.-A. and Maziade, M. Extension of the Generalized Disequilibrium Test to polytomous phenotypes and two locus models.

2. http://www.sph.umich.edu/csg/abecasis/Merlin/tour/input_files.html
}
\author{
Alexandre Bureau and Jordie Croteau
}
\seealso{
\code{\link{fat2Lpoly.withinR}}
}
\examples{
path.data=paste(.libPaths()[which(unlist(lapply(.libPaths(),function(x) length(grep("fat2Lpoly",dir(x)))))>0)],"/fat2Lpoly/extdata/",sep="")
if(length(path.data)>1) path.data=path.data[length(path.data)]

snps.anal=c("snp3.loc2","snp4.loc2")
microsat.names.loc2=c("2_3_mrk:","2_4_mrk:")

############################## design.endo2disease ###########################################
joint.tests=list(c(2,5))
snp.names.mat=cbind(rep("snp4.loc1",length(snps.anal)),snps.anal)
microsat.names.mat=cbind(rep("1_4_mrk:",length(snps.anal)),microsat.names.loc2)
test=fat2Lpoly(pedfilenames=paste(path.data,c("loc1.ped","loc2.ped"),sep=""),
               datfilenames=paste(path.data,c("loc1.dat","loc2.dat"),sep=""),
			   freq.data=paste(path.data,c("loc1.freq","loc2.freq"),sep=""),
               ibdfilenames=paste(path.data,c("loc1.ibd","loc2.ibd"),sep=""),
			   snp.names.mat=snp.names.mat,ibd.loci=microsat.names.mat,joint.tests=joint.tests,
               contingency.file=TRUE,design.constraint=design.endo2disease)

test$p.values.scores			   
# > test$p.values.scores
   # snp.cond  snp.test global_p params.joint_2-5_p param_1_score param_2_score 
# 1 snp4.loc1 snp2.loc2 5.80e-03           7.12e-01         0.449         0.333 
# 2 snp4.loc1 snp4.loc2 2.14e-07           1.24e-05         0.449         0.890 
# 3 snp4.loc1 snp5.loc2 1.14e-03           1.44e-01         0.449         0.776 
# 4 snp4.loc1 snp6.loc2 5.59e-04           3.84e-02         0.449        -0.082 
# 5 snp4.loc1 snp8.loc2 1.15e-03           1.55e-01         0.449         0.869 
# param_3_score param_4_score param_5_score param_1_p param_2_p param_3_p param_4_p
#        -1.427         3.638         0.733     0.653     0.739     0.154  0.000275
# 	     -1.427         3.638         4.612     0.653     0.373     0.154  0.000275
#        -1.427         3.638         1.785     0.653     0.438     0.154  0.000275
 #       -1.427         3.638         2.553     0.653     0.934     0.154  0.000275
#        -1.427         3.638         1.695     0.653     0.385     0.154  0.000275	   
  # param_5_p
# 1  0.464000
# 2  0.000004
# 3  0.074200
# 4  0.010700
# 5  0.090100
# > test$MA.table
   # snp.names.dat MA
# 1      snp1.loc1  1
# 2      snp2.loc1  1
# 3      snp3.loc1  1
# 4      snp4.loc1  2
# 5      snp5.loc1  1
# 6      snp6.loc1  1
# 7      snp7.loc1  1
# 8      snp8.loc1  1
# 9      snp1.loc2  1
# 10     snp2.loc2  1
# 11     snp3.loc2  1
# 12     snp4.loc2  2
# 13     snp5.loc2  1
# 14     snp6.loc2  1
# 15     snp7.loc2  1
# 16     snp8.loc2  1
####################################################################################

##################################### design.full ##############################################
joint.tests=list(c(2,3),c(5,6),c(8,9),c(2,3,5,6,8,9))
snp.names.mat=cbind(rep("snp4.loc1",length(snps.anal)),snps.anal)
microsat.names.mat=cbind(rep("1_4_mrk:",length(snps.anal)),microsat.names.loc2)
test=fat2Lpoly(pedfilenames=paste(path.data,c("loc1.ped","loc2.ped"),sep=""),
               datfilenames=paste(path.data,c("loc1.dat","loc2.dat"),sep=""),
			   freq.data=paste(path.data,c("loc1.freq","loc2.freq"),sep=""),
               ibdfilenames=paste(path.data,c("loc1.ibd","loc2.ibd"),sep=""),
			   snp.names.mat=snp.names.mat,ibd.loci=microsat.names.mat,joint.tests=joint.tests,
               design.constraint=design.full)

test$p.values.scores
# > test$p.values.scores
 #   snp.cond  snp.test global_p params.joint_2-3_p params.joint_5-6_p params.joint_8-9_p 
# 1 snp4.loc1 snp2.loc2 4.73e-02              0.896             0.4760           7.64e-01 
# 2 snp4.loc1 snp4.loc2 2.74e-08              0.196             0.3070           1.97e-08 
# 3 snp4.loc1 snp5.loc2 1.60e-02              0.627             0.5960           1.99e-01 
# 4 snp4.loc1 snp6.loc2 1.38e-03              0.391             0.4060           1.68e-02 
# 5 snp4.loc1 snp8.loc2 1.16e-03              0.405             0.0131           2.35e-01 
# params.joint_2-3-5-6-8-9_p param_1_score param_2_score param_3_score
#                   9.09e-01         0.449        -0.448         0.133
#                   1.04e-06         0.449        -1.463        -1.058				  
#                   5.01e-01         0.449        -0.165        -0.953
#                   6.08e-02         0.449        -0.819         1.098
#                   5.15e-02         0.449        -0.552        -1.226				  
  # param_4_score param_5_score param_6_score param_7_score param_8_score param_9_score 
# 1        -1.427        -0.050        -1.218         3.638        -0.023         0.733 
# 2        -1.427        -0.467        -1.465         3.638         3.770         4.612 
# 3        -1.427        -0.641         0.791         3.638        -0.195         1.785 
# 4        -1.427         0.833        -1.054         3.638         1.286         2.553 
# 5        -1.427        -0.821        -2.828         3.638         0.158         1.695 
# param_1_p param_2_p param_3_p param_4_p param_5_p param_6_p param_7_p param_8_p
#     0.653     0.654     0.894     0.154     0.960   0.22300  0.000275  0.982000
#     0.653     0.143     0.290     0.154     0.640   0.14300  0.000275  0.000163	
#     0.653     0.869     0.341     0.154     0.522   0.42900  0.000275  0.845000
#     0.653     0.413     0.272     0.154     0.405   0.29200  0.000275  0.198000
#     0.653     0.581     0.220     0.154     0.412   0.00468  0.000275  0.875000	
  # param_9_p
# 1  0.464000
# 2  0.000004
# 3  0.074200
# 4  0.010700
# 5  0.090100
###############################################################################################

########################################## design.1locus ######################################
snp.names.mat=as.matrix(snps.anal)
microsat.names.mat=as.matrix(microsat.names.loc2)
test=fat2Lpoly(pedfilenames=paste(path.data,"loc2.ped",sep=""),
               datfilenames=paste(path.data,"loc2.dat",sep=""),
               freq.data=paste(path.data,"loc2.freq",sep=""),
			   ibdfilenames=paste(path.data,"loc2.ibd",sep=""),
		       snp.names.mat=snp.names.mat,ibd.loci=microsat.names.mat,
			   design.constraint=design.1locus)

test$p.values.scores			   
# > test$p.values.scores
   # snp.test global_p param_1_score param_2_score param_3_score 
# 1 snp2.loc2  0.97500        -0.448        -0.050        -0.023 
# 2 snp4.loc2  0.00121        -1.463        -0.467         3.770 
# 3 snp5.loc2  0.91000        -0.165        -0.641        -0.195 
# 4 snp6.loc2  0.36000        -0.819         0.833         1.286 
# 5 snp8.loc2  0.78400        -0.552        -0.821         0.158 
# param_1_p param_2_p param_3_p
#     0.654     0.960  0.982000
#     0.143     0.640  0.000163
#	  0.869     0.522  0.845000
#     0.413     0.405  0.198000
#     0.581     0.412  0.875000			
################################################################################################

################################## design.dichotomique #########################################
joint.tests=list(c(2,3))
snp.names.mat=cbind(rep("snp4.loc1",length(snps.anal)),snps.anal)
microsat.names.mat=cbind(rep("1_4_mrk:",length(snps.anal)),microsat.names.loc2)
test=fat2Lpoly(pedfilenames=paste(path.data,c("loc1.ped","loc2.ped"),sep=""),
               datfilenames=paste(path.data,c("loc1.dat","loc2.dat"),sep=""),
			   freq.data=paste(path.data,c("loc1.freq","loc2.freq"),sep=""),
               ibdfilenames=paste(path.data,c("loc1.ibd","loc2.ibd"),sep=""),
			   snp.names.mat=snp.names.mat,ibd.loci=microsat.names.mat,joint.tests=joint.tests,
               design.constraint=design.dichotomique)

test$p.values.scores			   
# > test$p.values.scores
   # snp.cond  snp.test global_p params.joint_2-3_p param_1_score param_2_score  
# 1 snp4.loc1 snp2.loc2 6.94e-02           9.95e-01         2.659        -0.050  
# 2 snp4.loc1 snp4.loc2 8.00e-07           5.97e-06         2.659         3.013 
# 3 snp4.loc1 snp5.loc2 8.06e-03           9.33e-02         2.659        -0.566 
# 4 snp4.loc1 snp6.loc2 3.03e-03           3.26e-02         2.659         1.680 
# 5 snp4.loc1 snp8.loc2 6.57e-02           9.35e-01         2.659        -0.345
# param_3_score param_1_p param_2_p param_3_p
#         0.085   0.00784   0.96000  0.932000
#         3.870   0.00784   0.00259  0.000109
#         2.103	  0.00784   0.57100  0.035500
#         2.006	  0.00784   0.09300  0.044900
#         0.124   0.00784   0.73000  0.902000		
###############################################################################################
}




