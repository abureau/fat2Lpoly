
# Jordie Croteau
# 17 août 2012

# correspond au fichier read_merlin_files_v3.R dans le dossier "programmes"

# modifié le 21 août pour avoir la possibilité de prendre comme argument des vecteurs d'allèles mineures au lieu de fichiers de fréquences d'allèles.

# modifié le 23 août pour intégrer la lecture des données d'IBD.

# 4 avril 2013: valeur NULL par défaut ajoutée pour l'argument ibdfilenames (dans ce cas, aucune lecture des données d'IBD)

read.merlin.files=function(pedfilenames,datfilenames,freq.data,ibdfilenames=NULL)
{
###################### Définition des arguments #####################################################################################
# pedfilenames : vecteur des noms de fichiers ped (un fichier par locus).  Les sujets inclus peuvent être un sous-ensemble de ceux 
#                inclus dans les fichiers d'IDB.
# datfilenames : vecteur des noms de fichiers dat (un fichier par locus). 
# freqfilenames: vecteur des noms de fichiers freq (un fichier par locus). 
# ibdfilenames: 
# Tous ces fichiers doivent être en format Merlin (voir http://www.sph.umich.edu/csg/abecasis/Merlin/tour/input_files.html pour la description détaillée de ce format)
#####################################################################################################################################

# extraction des fam.id, subject.ids, y1 et y2, à partir du premier fichier ped.
ped1.tmp=read.table(pedfilenames[1],header=FALSE,as.is=TRUE)
if(any(apply(ped1.tmp,2,is.character))) stop(paste(pedfilenames[1],"contains letters in some fields.  All fields must be numeric."))
fam.id=ped1.tmp[,1]
subject.ids=ped1.tmp[,2]
y1=ped1.tmp[,6]
y2=ped1.tmp[,7]
y1[is.na(y1)]=0
y2[is.na(y2)]=0

################ extraction des génotypes de tous les locus #####################################
ped=data.frame(fam.id,subject.ids,ped1.tmp[,3:4],y1,y2,ped1.tmp[,8:ncol(ped1.tmp)])
n.loc=length(pedfilenames)

if(n.loc>1)
 {
  for(loc.num in 2:n.loc)
   {
    ped.tmp=read.table(pedfilenames[loc.num],header=FALSE,as.is=TRUE)
	if(any(apply(ped.tmp,2,is.character))) stop(paste(pedfilenames[loc.num],"contains letters in some fields.  All fields must be numeric."))
    ped.tmp=ped.tmp[,c(1,2,8:ncol(ped.tmp))]
	colnames(ped.tmp)[1:2]=c("fam.id","subject.ids")
	ped=merge(ped,ped.tmp,by=c("fam.id","subject.ids"),all.x=FALSE,all.y=FALSE,sort=FALSE)
   }
 }
pere=ped[,3]
mere=ped[,4]
ped=ped[,-(3:4)]
if(nrow(ped)<nrow(ped1.tmp)) warning(paste("Subjects from 1 or more ped files differ from those of other ped files. Only subjects found in the",n.loc,"ped files at the same time are kept for analysis."))
##################################################################################################
 
############### lecture des fichiers datfile pour obtenir les noms de SNPs #######################
dat1=read.table(datfilenames[1],as.is=TRUE)
if(dat1[2,1]!="A"|dat1[1,1]!="A") stop(paste("The first and second lines of",datfilenames[1],"should be starting by 'A'."))
if(!all(dat1[3:nrow(dat1),1]=="M")) stop(paste("Lines number 3 to",nrow(dat1),"of",datfilenames[1],"should all be starting by 'M'."))

pheno.name=dat1[2,2]
endo.name=dat1[1,2]
cat("\n")
cat("Y1 data extracted from input files:   ",endo.name,"\n")
cat("Y2 data extracted from input files:   ",pheno.name,"\n")
cat("\n")

snp.names.dat=dat1[dat1[,1]=="M",2]

if(n.loc>1)
 {
  for(loc.num in 2:n.loc)
   {
    dat.tmp=read.table(datfilenames[loc.num],as.is=TRUE)
	if(dat.tmp[2,1]!="A"|dat.tmp[1,1]!="A") stop(paste("The first and second lines of",datfilenames[loc.num],"should be starting by 'A'."))
    if(!all(dat.tmp[3:nrow(dat.tmp),1]=="M")) stop(paste("Lines number 3 to",nrow(dat.tmp),"of",datfilenames[loc.num],"should all be starting by 'M'."))
    snp.names.dat=c(snp.names.dat,dat.tmp[dat.tmp[,1]=="M",2])
   }
 }
#################################################################################################

################################################ obtention des allèles mineures ##########################################################
if(is.list(freq.data)) MA=unlist(freq.data)
else
{ 
######################### lecture des fichiers freq pour obtenir les allèles mineures  ###########################################

# il y a 2 formats de fichiers freq possibles dans merlin (avec fréquences d'allèles une sous l'autre ou une à côté de l'autre).
# Donc, pour chaque fichier freq, il faut d'abord détecter quel est le format.
freq=readLines(freq.data[1])
n.lines=length(freq)
greps.M=grep("M",freq)
format.freq=as.numeric(length(greps.M)==(n.lines/2))

name.lines=freq[greps.M]
write(name.lines,"name.lines_asdhskjhfdak.txt")
snp.names.freq=as.character(read.table("name.lines_asdhskjhfdak.txt",as.is=TRUE)[,2])
unlink("name.lines_asdhskjhfdak.txt")

if(format.freq==0)
 {
  freq.col=read.table(freq.data[1],as.is=TRUE)[,2]
  freqs=vector(mode="list",length=length(greps.M)) 
  greps.M.mod=c(greps.M,length(freq.col)+1)
  for(j in 1:(length(greps.M.mod)-1)) freqs[[j]]=as.numeric(freq.col[(greps.M.mod[j]+1):(greps.M.mod[j+1]-1)])
 }
 
if(format.freq==1) freqs=lapply(strsplit(freq[-greps.M],split=" "),function(x){ x=as.numeric(x); return(x[!is.na(x)])})

# allèles mineures
MA=unlist(lapply(freqs,function(x) {y=x[x>0]; ifelse(length(y)==2,which(x==min(y)),NA)}))
# allèles majeures
MjA=unlist(lapply(freqs,which.max))

if(n.loc>1)
 {
  for(loc.num in 2:n.loc)
   {
    freq=readLines(freq.data[loc.num])
    n.lines=length(freq)
    greps.M=grep("M",freq)
    format.freq=as.numeric(length(greps.M)==(n.lines/2))

    name.lines=freq[greps.M]
    write(name.lines,"name.lines_asdhskjhfdak.txt")
    snp.names.freq.tmp=as.character(read.table("name.lines_asdhskjhfdak.txt",as.is=TRUE)[,2])
    snp.names.freq=c(snp.names.freq,snp.names.freq.tmp)

    unlink("name.lines_asdhskjhfdak.txt")

    if(format.freq==0)
     {
      freq.col=read.table(freq.data[loc.num],as.is=TRUE)[,2]
      freqs=vector(mode="list",length=length(greps.M)) 
      greps.M.mod=c(greps.M,length(freq.col)+1)
      for(j in 1:(length(greps.M.mod)-1)) freqs[[j]]=as.numeric(freq.col[(greps.M.mod[j]+1):(greps.M.mod[j+1]-1)])
     }
 
    if(format.freq==1) freqs=lapply(strsplit(freq[-greps.M],split=" "),function(x){ x=as.numeric(x); return(x[!is.na(x)])})
	
	# allèles mineures
    MA=c(MA,unlist(lapply(freqs,function(x) {y=x[x>0]; ifelse(length(y)==2,which(x==min(y)),NA)})))
    # allèles majeures
    MjA=c(MjA,unlist(lapply(freqs,which.max)))
   }
 }
if(any(is.na(MA))) stop("one or more SNP(s) with unknown minor allele, probably because of nul minor allele frequency(ies)")
 
if(!all(snp.names.freq==snp.names.dat)) stop("SNP names or order in freq files not the same as in the data files")

# vérifier que toutes les colonnes de génotypes ne contiennent pas d'autres allèles que les allèles mineures et majeures obtenues 
# à partir des fichiers freq.
for(j in 1:length(MA))
 {
  genos.tmp=as.vector(ped[,c(3+2*j,4+2*j)])
  if(!all(is.na(genos.tmp)|genos.tmp==0))
   {
    genos.tmp=genos.tmp[!(is.na(genos.tmp)|genos.tmp==0)]
    if(sum(!(genos.tmp%in%c(MA[j],MjA[j])))!=0) stop(paste("one or more alleles differ from minor and major alleles expected from freq files, for SNP",snp.names.dat[j]))
   }
 }
################################################################################################################################
}

# conversion des génotypes en valeurs 0,0.5,1 (mode="allelic")
x.all=alleles2sums(ped[,5:ncol(ped)],MA.vec=MA,snp.names=snp.names.dat,mode="allelic")

MA.table=data.frame(snp.names.dat,MA)

# lecture des fichiers d'IBD
ibd.dat.list=vector(mode="list",length=n.loc)
if(!is.null(ibdfilenames)) for(loc.num in 1:n.loc) ibd.dat.list[[loc.num]]=as.data.frame(read.table(ibdfilenames[loc.num],header=TRUE,as.is=TRUE))

list(ped=ped[,1:4],x.all=x.all,MA.table=MA.table,ibd.dat.list=ibd.dat.list,y1.name=endo.name,y2.name=pheno.name,ibdfilenames=ibdfilenames,pere=pere,mere=mere)
}                      