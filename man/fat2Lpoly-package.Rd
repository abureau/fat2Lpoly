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
path.data=paste(.libPaths()[which(unlist(lapply(.libPaths(),
function(x) length(grep("fat2Lpoly",dir(x)))))>0)],"/fat2Lpoly/extdata/",sep="")
if(length(path.data)>1) path.data=path.data[length(path.data)]

snps.anal=c("snp3.loc2","snp4.loc2")
microsat.names.loc2=c("2_3_mrk:","2_4_mrk:")

################# design.endo2disease with conditioning on locus 1 #############
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
################################################################################  
}
