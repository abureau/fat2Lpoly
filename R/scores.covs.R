# Fonction globale qui calcule les scores et les covariances pour toutes les familles et tous les locus.
# ATTENTION, MÊME SI CERTAINES PARTIES DE CE PROGRAMME FONCTIONNENT POUR UN NOMBRE ARBITRAIRE DE LOCUS,
# D'AUTRES PARTIES NE SONT VALIDES QUE POUR UN OU 2 LOCUS.  PAR EXEMPLE, LES LIGNES 255, 258 ET 283 DANS LA VERSION 12, 
# POUR LE MOMENT (EN DATE DU 18 MAI 2012) NE FONCTIONNENT BIEN QUE POUR LES CAS DE 1 OU 2 LOCUS.

# correspond au fichier fonction_scores_covs_v2.R dans le dossier "programmes"


# Jordie Croteau
# 11 mars 2011

# Modifiée par Alexandre Bureau
# 24 mars 2011

# Ajout de l'argument design.constraint pour créer les matrices de design à l'INTÉRIEUR chaque catégorie 
# pour des contraintes spécifiques à chaque catégorie, et incluant possiblement des termes d'interaction.
# Retrait de l'argument x.loc. Les listes de locus impliqués dans chaque terme sont maintenant produites
# par la fonction design.constraint.

# Re-modifiée par Jordie Croteau
# 25 mars 2011

# Re-modifiée par Jordie Croteau (seulement les commentaires)
# 31 mars 2011

# Re-modifiée par Jordie Croteau
# 1er avril 2011

# Ajout de la gestion des analyses quand certains sujets ont un phénotype inconnu.
# Ajout aussi d'un traitement séparé des familles dont tous les sujets sont dans 
# la même catégorie phénotypique. Ces familles ne contribuent ni au score, ni à sa variance.
# Cependant, pour éviter les divisions par 0, une petite variance de 1e-10 est comptabilisée
# (Il faudra vérifier quand ces divisions par 0 se produisent)

# Re-modifié par Alexandre Bureau
# 4 avril 2011

# Correction d'une erreur qui faisait que les calculs de variance pour toutes les catégories se faisaient  
# tous avec les termes de variance par locus pour la catégorie 1
# Ajout de la gestion des contraintes inter-catégories spécifiées par par.constrained et constraints.

# Re-modifié par Alexandre Bureau
# 29 avril 2011

# Correction d'une erreur dans dans le recodage des y1 et y2 en catégories 1 à 4

# Re-modifié par Alexandre Bureau
# 9 juin 2011

# Correction d'une erreur dans le calcul de la variance des termes d'interaction: il faut calculer la somme des variances du produit 
# des scores de chaque locus par paires et non le produit des variances de la somme des scores de chaque locus.
# Implantation du calcul des covariances entre scores pour le même locus dans des fonctions de régression distincte

# Re-modifié par Alexandre Bureau
# 20 janvier 2012
# Correction des longueurs des vecteurs ind.par et ind.cat 

# Re-modifié par Jordie les 17 et 18 mai 2012 pour corriger quelques "degenerate cases".

# Re-modifié par Jordie le 21 mai 2012 pour permettre d’avoir des fichiers de pedigree avec des ensembles de sujets différents 
# (ou seulement dans un ordre différent) d’un fichier de pedigree à l’autre.  Un avertissement est toutefois affiché si 
# certains ensembles de sujets diffèrent.  Dans ce cas, seulement les sujets en commun à tous les fichiers de pedigree
# sont conservés.  De plus, on exclut les sujets dont le génotype d'un des SNPs à analyser est manquant.  On exclut 
# aussi les sujets dont le phénotype ou l'endophénotype est manquant (valeur 0).

# Re-modifié par Alexandre Bureau
# 23 mai 2012
# ajout de conditions "if (n.loc>1)" pour éviter des opérations qui ne s'appliquent pas avec 1 seul locus
# S'il y a des termes d'IBD qui sont NA parce qu'une catégorie n'est pas représentée, une correction a été apportée 
# pour que la contribution à la covariance soit 0, pas NA

# Re-modifié par Alexandre Bureau
# 25 mai 2012
# Changement à l'appel de cov.score.interfunction pour passer de la version 2 à la version 3 de cette fonction

# Re-modifié par Alexandre Bureau
# 5 juin 2012

# Implantation d'une nouvelle façon d'estimer la variance des termes de produit

# Re-modifié le 25 juin 2012 par Jordie:
# partie du début (lecture des données) reléguée à une fonction séparée: read.merlin.files,
# où on lit les données de format merlin.

# petites modifs par Jordie le 17 août 2012.

# modifié par Jordie le 22 août pour enlever des arguments superflus.

scores.covs=function(subject.ids,fam.id,y,n.levels,ibd.dat,n.loc,xp,xp.loc,xl,il,xibd.loc,ind.par,rep.par,ind.catl,ind.cat,contingency.file,descrip.file,calculpoids=calcule.poids.alphafixe,lc=NULL,alpha.vec=rep(0,n.levels-1))
{
###################### Définition des arguments #####################################################################################
# subject.ids: vecteur des id de chaque sujet
# fam.id: vecteur des id de la famille de chaque sujet
# y: vecteur de type "factor" donnant le statut, qui est la combinaison des 6e colonne (endophénotype) et 7e colonne (phénotype) des fichiers ped
# n.levels : nombre de catégories du statut (y)
# ibd.dat : data frame des données d'IBD combinées
# n.loc : nombre de locus
# xp : matrice de design pour le calcul du score
# xp.loc : Liste du vecteur de locus impliqués dans chaque paramètre
# xl : matrice de design pour le calcul des covariances (contient seulement les effets principaux des locus)
# il : indices après conversion des termes de produits en indices des variables dans le produit
# ind.par : donne les indices des locus pour la catégorie à laquelle chaque terme appartient
# lc : locus sur lequel on conditionne le test du score
# alpha.vec : vecteur de log rapports de cote entre phénotype et compte d'allèle au locus précisé par lc, un élément pour chaque catégorie k vs. les autres catégories
# autres paramètres à commenter..............
#####################################################################################################################################


######################################## données d'IBD ##########################################################
# Prendre le produit des IBD des paires de locus pour le calcul de la variance des termes d'interaction
# s'il y a plus d'un locus et il y a présence de termes d'interaction.
# Attention! Pour l'instant, soit tous les termes d'interaction sont présents, soit il y en a aucun
l1=ibd.dat$ID1
l2=ibd.dat$ID2
fam.id.ibd=ibd.dat$FAMILY
pim=as.matrix(ibd.dat[,4:ncol(ibd.dat)])

if (n.loc>1 & max(xibd.loc) > n.loc)
{
pim2 <- apply(pim,1,produits.paires)

# Combiner les termes originaux et les produits en une seule matrice
# Jordie, 17 mai 2012: ajout de transposes pour avoir la matrice dans le bon sens lorsqu'il y a plus de 2 locus.
pim <- t(rbind(t(pim),pim2))
}
#################################################################################################################

#################### calculs des scores et covariances pour toutes les familles #####################################################
fam.id.u=unique(fam.id)
nb.fam=length(fam.id.u)

ibd.terms.mat=sigma2i.mat=array(NA,c(nb.fam,ncol(pim),n.levels-1,n.levels-1))
sigma2.mat=array(NA,c(nb.fam,dim(xl)[2],n.levels-1,n.levels-1))
scores.mat=array(NA,c(nb.fam,dim(xp)[2]))

for(i in 1:nb.fam)
 {
  indices.y=fam.id==fam.id.u[i]
  if(contingency.file)
   {
    cat("family",fam.id.u[i],"\n",file=descrip.file,append=TRUE)
    cat(table(y[indices.y]),"\n",file=descrip.file,append=TRUE)
    cat("\n",file=descrip.file,append=TRUE)
   }
  # Si plus d'une catégorie de y est représentée dans la famille
  if(length(unique(y[indices.y]))>1)
   {
    # nombre de sujets dans la famille
    ni = sum(indices.y)
 
    # dans le cas où l'analyse est faite avec un sous-ensemble des sujets de la famille,
    # limiter les données d'IBD à ce sous-ensemble de sujets.
    sub.tmp=subject.ids[indices.y]
    indices.ibd=fam.id.ibd==fam.id.u[i] & l1%in%sub.tmp & l2%in%sub.tmp
  
	  sigma2.mat[i,,,]=cov.score.poly(array(xl[indices.y,,],c(ni,dim(xl)[2],dim(xl)[3])),y[indices.y],subject.ids[indices.y],l1[indices.ibd],l2[indices.ibd],array(pim[indices.ibd,],c(sum(indices.ibd),ncol(pim))),x.loc=xibd.loc)
	  # Les covariances entre scores pour le même locus dans différente fonctions de régression impliquent des estimations de la variance des scores
	  # On les copie de sigma2.mat dans les éléments appropriés de sigma2i.mat
    # Appel de la version 4
    sigma2i.mat[i,,,]=cov.score.interfunction(xibd.loc,ind.catl,array(sigma2.mat[i,,,],dim(sigma2.mat)[2:4]))
    
    # Si on ne conditionne pas sur un locus  
    if (is.null(lc))
      {
      ibd.terms.mat[i,,,]=ibd.terms(y[indices.y],subject.ids[indices.y],l1[indices.ibd],l2[indices.ibd],array(pim[indices.ibd,],c(sum(indices.ibd),ncol(pim))))
      scores.mat[i,]=score.poly(array(xp[indices.y,,],c(ni,dim(xp)[2],dim(xp)[3])),y[indices.y])
      }
    # Sinon on conditionne sur un locus en pondérant
    else
      {
      if (!(lc %in% 1:dim(xl)[2])) stop("The index lc does not correspond to a valid locus index in",1:dim(xl)[2])
      if (length(alpha.vec)!=dim(xl)[3]) stop("The number of alpha coefficients for the computation of weights does not equal the 3rd dimension of xl (",dim(xl)[3],")")
      # Calcul des poids pour la famille
      # On suppose que la matrice de design pour le niveau 1 contient le locus lc
      w = calculpoids(array(xl[indices.y,,],c(ni,dim(xl)[2],dim(xl)[3])),y[indices.y],ind.par,rep.par,alpha.vec,lc,klc=1)      
      ibd.terms.mat[i,,,]=ibd.terms.w(y[indices.y],subject.ids[indices.y],l1[indices.ibd],l2[indices.ibd],array(pim[indices.ibd,],c(sum(indices.ibd),ncol(pim))),w)
      scores.mat[i,]=score.poly.w(array(xp[indices.y,,],c(ni,dim(xp)[2],dim(xp)[3])),y[indices.y],w)      
      }
    }
  # s'il y a une seule catégorie de y représentée dans la famille, les programmes ne s'appliquent pas 
  # et on donne une valeur de 0 au score et une valeur très petite à ibd.terms.mat et sigma2.mat.
  else 
   {
	ibd.terms.mat[i,,,]=1e-10
    scores.mat[i,]=0
   }
 }

# la matrice (array) des covariances des effets principaux est la somme sur les deux dimensions de catégories
# des produits des termes d'ibd et des sigma2.
#petite modif Jordie, 24 mars (il y avait une erreur si une ou des dimensions sont de longueur 1):
#ibd.terms.xl.mat=array(ibd.terms.mat[,xl.loc,,],c(dim(ibd.terms.mat)[1],length(xl.loc),dim(ibd.terms.mat)[3],dim(ibd.terms.mat)[4]))
ibd.terms.xl.mat=array(ibd.terms.mat[,xibd.loc,,],c(dim(ibd.terms.mat)[1],length(xibd.loc),dim(ibd.terms.mat)[3],dim(ibd.terms.mat)[4]))
cov.mat.l=apply(ibd.terms.xl.mat*sigma2.mat,1:2,sum,na.rm=TRUE)

# Calcul de covariances entre les différentes catégories pour les mêmes locus
#     Hypothèse: les matrices de design sont les mêmes pour différentes catégories (ce qu'on obtient 
#     s'il n'y a pas de contrainte)
# Il faudra vérifier si on a besoin de forcer un array s'il y a un seul locus
cov.inter.mat <- ibd.terms.mat*sigma2i.mat
# S'il y a des termes d'IBD qui sont NA parce qu'une catégorie n'est pas représentée, la contribution 
# à la covariance doit être 0, pas NA
cov.inter.mat[is.na(cov.inter.mat)] = 0
# Vérification que tous les termes de produits sont présents dans sigm2i.mat
if (dim(sigma2i.mat)[2] > n.loc & dim(sigma2i.mat)[2] != n.loc*(n.loc+1)/2 ) stop ("sigma2i.mat does not contain variances for all product terms.")
####################################################################################################################################################

# On obtient la matrice des covariances de tous les effets en prenant la combinaison des termes appropriés
cov.mat <- array(0,c(nb.fam,dim(xp)[2],dim(xp)[2]))
# Boucle sur les paramètres
for (j in 1:length(xp.loc))
  {
  # Boucle sur les termes associés à chaque paramètre (ordinairement un seul)
  for (h in 1:length(xp.loc[[j]]))
    {
	# Extraction des locus impliqués dans le terme h
    il <- as.numeric(unlist(strsplit(xp.loc[[j]][h],split="")))
    # On va chercher la variance du produit des locus impliqués dans l'élément h du terme xp.loc[[j]]
	# Les indices des locus doivent s'interpréter dans la portion de la matrice de covariance
	# qui s'applique aux paramètres du terme xp.loc[[j]]
	# Pour l'instant, on suppose qu'il y a juste 2 locus
    if (length(il) > 1) 
	{   
    # AB: Voici ma solution pour le cas d'un seul terme d'interaction qui est considéré dans le programme.
    # indice où se trouve le terme d'IBD
    # calculé selon la formule n.loc + somme_i=1^(il[1]-1) (n.loc - i) + il[2] - il[1]
    # (on additionne n.loc pour sauter les termes d'IBD pour les effets principaux des locus)
    ii = n.loc + ifelse(il[1]>1,(il[1] - 1)*(n.loc - il[1]/2),0) + il[2] - il[1]
    cov.tmp <- cov.mat.l[,ind.par[[j]][ii]]
	}
    else cov.tmp <- cov.mat.l[,ind.par[[j]][il]]
	# On additionne le terme courant à la valeur de la covariance
    cov.mat[,j,j] <- cov.mat[,j,j] + cov.tmp
    }
  # s'il y a plus d'un paramètre
  if (length(xp.loc) > 1 & j > 1)
	{
	for (jj in 1:(j-1))
	  {
	  # On trouve l'intersection entre les termes des 2 listes
	  inter <- intersect(xp.loc[[j]],xp.loc[[jj]])
	  
	  if (length(inter) > 0) 
		{
        for (h in 1:length(inter))
          {
		  # On laisse tomber complètement pour l'instant la covariance entre termes d'une même fonction de régression
		  		  # Extraction des locus impliqués dans le terme h
		  il <- as.numeric(unlist(strsplit(inter[h],split="")))

	      # Si les indices de paramètres réfèrent à des catégories différentes 
          if (ind.cat[j] != ind.cat[jj])
            {
	        if (length(il) > 1) 
	        {
          # AB: Voici ma solution pour le cas d'un seul terme d'interaction qui est considéré dans le programme.
          # indice où se trouve le terme d'IBD
          # calculé selon la formule n.loc + somme_i=1^(il[1]-1) (n.loc - i) + il[2] - il[1]
          # (on additionne n.loc pour sauter les termes d'IBD pour les effets principaux des locus)
          ii = n.loc + ifelse(il[1]>1,(il[1] - 1)*(n.loc - il[1]/2),0) + il[2] - il[1]
          
          cov.tmp <- cov.inter.mat[,ii,ind.cat[j],ind.cat[jj]]
          }
	        # Sinon il y a un seul locus dans l'intersection
	        else cov.tmp <- cov.inter.mat[,il,ind.cat[j],ind.cat[jj]]
	        # On additionne le terme courant à la valeur de la covariance
            cov.mat[,j,jj] <- cov.mat[,j,jj] + cov.tmp
		    }		  
	      }
		}
	  cov.mat[,jj,j] <- cov.mat[,j,jj]
	  }
	}
  }
list(scores.mat=scores.mat,cov.mat=cov.mat)
}
