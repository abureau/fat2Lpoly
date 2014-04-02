\name{read.merlin.files}
\alias{read.merlin.files}

\title{
function to read input files in Merlin format
}
\description{
Reads the pedigree, data and allele frequency input files.  The data read is reformatted to be used by the function \code{fat2Lpoly.withinR}.
}
\usage{
read.merlin.files(pedfilenames, datfilenames, freq.data, ibdfilenames = NULL)
}
\arguments{
   \item{pedfilenames}{
vector of 1 or 2 (the number of loci involved in the \code{design} function) character strings giving the names of the pedigree files in Merlin format (see Merlin website [1]). Put the full path of the files if they are not in the current working directory. If the phenotype is polytomous with 4 levels created by all combinations of two dichotomous phenotypic variables \eqn{Y_{1}}{Y[1]} and \eqn{Y_{2}}{Y[2]}, then the sixth and seventh columns of each file contain respectively \eqn{Y_{1}}{Y[1]} (e.g. the endophenotype) and \eqn{Y_{2}}{Y[2]} (e.g. the disease phenotype). If the phenotype is dichotomous, then the sixth column of each file contains it.}
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
vector of 1 or 2 (the number of loci involved in the \code{design} function) character strings giving the names of the Merlin format ibd files corresponding to the pedigree files. If NULL (the default), the reading of the IBD files is skipped.
}
}
\details{
All subjects included in the pedigree files must also be found in the IBD files.

All fields in the pedigree files must be numeric.  No letters allowed, even for family and subject ID's.
}
\value{
returns a list of six objects:

  \item{ped}{
data frame with columns fam.id, subject.ids, endophenotype and phenotype (in the given order)
}
 \item{x.all}{
data frame of SNP genotypes in the format "(number of minor alleles)/2", for all SNPs listed in the file(s) in \code{datfilenames}.  It contains only the SNP data and it has as column names the SNP names in \code{datfilenames}.  The lines come in the same order as in \code{ped}.
}
 \item{MA.table}{
data frame giving the minor allele numbers of all the SNPs.  The first column consists of \code{x.all}'s column names and the second column the minor allele numbers.
}
 \item{ibd.dat.list}{
list of one or two data frames containing the columns of the IBD data file(s) in \code{ibdfilenames}.  
}
 \item{y1.name}{
affection name extracted from first line of the data file(s)
}
 \item{y2.name}{
affection name extracted from second line of the data file(s)
}
  \item{ibdfilenames}{
(same object as provided as argument) vector of 1 or 2 (the number of loci involved in the \code{design} function) character strings giving the names of the Merlin format ibd files corresponding to the pedigree files.
}
}
\references{
1. http://www.sph.umich.edu/csg/abecasis/Merlin/tour/input_files.html
}
\author{
Alexandre Bureau and Jordie Croteau
}
\seealso{
\code{\link{fat2Lpoly.withinR}}
}
\examples{
path.data=paste(.libPaths()[which(unlist(lapply(.libPaths(),
function(x) length(grep("fat2Lpoly",dir(x)))))>0)],
"/fat2Lpoly/extdata/",sep="")
if(length(path.data)>1) path.data=path.data[length(path.data)]

input.data=read.merlin.files(pedfilenames=
                  paste(path.data,c("loc1.ped","loc2.ped"),sep=""),
                  datfilenames=
				paste(path.data,c("loc1.dat","loc2.dat"),sep=""),
			      freq.data=
			      paste(path.data,c("loc1.freq","loc2.freq"),sep=""),
                       ibdfilenames=
				paste(path.data,c("loc1.ibd","loc2.ibd"),sep=""))

input.data2=read.merlin.files(pedfilenames=
					paste(path.data,"loc2.ped",sep=""),
                            datfilenames=
                            paste(path.data,"loc2.dat",sep=""),
                            freq.data=
                            paste(path.data,"loc2.freq",sep=""),
					ibdfilenames=
					paste(path.data,"loc2.ibd",sep=""))
}
