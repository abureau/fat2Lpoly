\name{ped.x.all}
\alias{ped.x.all}
\docType{data}
\title{
Example dataset returned by the function \code{read.merlin.files}
}
\description{
This list is an example of output from the function \code{read.merlin.files}. It is provided to test the function \code{fat2Lpoly.withinR} by executing the example code in the \code{fat2Lpoly.withinR} documentation.
}
\usage{data(ped.x.all)}
\format{
A list of six objects:
  \describe{
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
}
\examples{
data(ped.x.all)
}
\keyword{datasets}
