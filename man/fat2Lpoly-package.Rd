\name{fat2Lpoly-package}
\alias{fat2Lpoly-package}
\docType{package}
\title{
Two-locus Family-based Association Test with Polytomous Outcome
}
\description{
Performs 2-locus and 1-locus family-based association tests defined by some design constraints (see [1]).  
}
\details{
\tabular{ll}{
Package: \tab fat2Lpoly\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-08-31\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}

}
\author{
Alexandre BUREAU <alexandre.bureau@msp.ulaval.ca> and Jordie Croteau <jordie.croteau@crulrg.ulaval.ca> \cr
Maintainer: Alexandre BUREAU <alexandre.bureau@msp.ulaval.ca>
}
\references{
1. paper reference
}
\keyword{ package }

\examples{
path.data=paste(.libPaths()[which(unlist(lapply(.libPaths(),function(x) length(grep("fat2Lpoly",dir(x)))))>0)],"/fat2Lpoly/extdata/",sep="")
if(length(path.data)>1) path.data=path.data[length(path.data)]

snps.anal=c("snp3.loc2","snp4.loc2")
microsat.names.loc2=c("2_3_mrk:","2_4_mrk:")

##################################### design.endo2disease #########################################
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
#        -1.427         3.638         4.612     0.653     0.373     0.154  0.000275
#        -1.427         3.638         1.785     0.653     0.438     0.154  0.000275	   
#        -1.427         3.638         2.553     0.653     0.934     0.154  0.000275
#        -1.427         3.638         1.695     0.653     0.385     0.154  0.000275	   
  # param_5_p
# 1  0.464000
# 2  0.000004
# 3  0.074200
# 4  0.010700
# 5  0.090100
}
