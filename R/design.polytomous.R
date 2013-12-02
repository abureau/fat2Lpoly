# Fonction pour créer matrice de design 

# correspond au fichier design_polytomique_v7.r dans le dossier "programmes"


# Version 2: 
# - la matrice de données x est mainteant spécifique à chaque dimension (et peut donc varier). 
# - la fonction traite maintenant la liste des locus impliqués dans chaque effet

# Version 3: 
# Cette nouvelle version résout un problème avec l’option « constraints » qui faisait que la vérification 
# qu'un paramètre n'est pas impliqué dans plus d'une contrainte était incorrecte. 
# Le problème n’affectait pas les cas sans contraintes inter-catégories.

# Version 4:
# Ajout de la liste ind.vec qui donne pour chaque paramètre les indices des paramètres (locus) de la catégorie à laquelle il appartient.
# Attention! avec cette version, les indices ne sont pas dans l'ordre des paramètres. À corriger.

# Version 5:
# Ajout de l'argument rep.par qui permet de passer un nombre de paramètres différent du nombre de
# colonnes de x pour le nombre de répétition des indices de paramètres du vecteur ind.par. 
# Ceci est requis quand on crée la matrice de design avec seulement les effets principaux
# mais qu'on veut garder un vecteur ind.par de longueur différente au nombre de colonnes
# dans  la matrice de design

# par Alexandre Bureau
# janvier 2012

# Version 6:
# Jordie, 23 mai 2012: modifications pour correction d'un bug potentiel aux lignes 125 et 177.

# Version 7:
# Changement dans la déclaration de ind.par pour avoir toujours la bonne longueur


# par Alexandre Bureau
# juin 2012

design.polytomous <- function(x,K,x.loc.in,par.constrained,constraints,rep.par=NULL)
{
# x: liste de matrices de design pour un modèle logistique entre chaque catégorie k de 1 à K-1  et la catégorie K.
#     Il s'agit du modèle le plus général.
# K: nombre de catégories de la variable réponse
# x.loc.in : liste des vecteurs énumérant les locus impliqués dans chaque colonne de x
# par.constrained : matrice des indices du paramètre impliqué dans chaque contrainte (nc colonnes) pour chaque catégorie (K-1 lignes)
# constraints: matrice (K-1) x nc spécifiant des contraintes entre les paramètres des modèles logistiques
#              pour différents niveaux de la variable réponse, une contrainte par colonne. La valeur 0 veut dire que le paramètre n'est pas impliqué.
# rep.par: vecteur du nombre de paramètres associé à chaque matrice x  (pas nécessairement égal au nombre de colonnes de x)
if (K<2) stop ("K must be >= 2.")

if (length(unique(unlist(lapply(x,nrow)))) > 1) stop("All matrices in x must have the same number of rows.")

# np est le vecteur du nombre de paramètres dans chaque matrice x
np <- sapply(x,ncol)

if (missing(constraints))
  {
  # Cas particulier sans contrainte: on copie simplement la matrice x sur les K - 1 tranches de xp,
  # décallée à chaque tranche.
  # La liste des indices où sont situés les paramètres pour la catégorie 1 est 1 à np[1]
  if (is.null(rep.par))
    {
    ind.par <- vector("list",length=sum(np))
    ind.par[1:np[1]] <- list(1:np[1])
	}
  else
    { 
    ind.par <- vector("list",length=length(rep.par))
    ind.par[1:rep.par[1]] <- list(1:np[1])
    }
  
  xp <- array(0,c(nrow(x[[1]]),sum(np),K-1))
  xp[,1:np[1],1] <- x[[1]]
  if (K > 2)
  {
  for (k in 2:(K-1))
  {
    xp[,(sum(np[1:k-1])+1):sum(np[1:k]),k] <- x[[k]]
    if (is.null(rep.par))	ind.par[(sum(np[1:k-1])+1):sum(np[1:k])] <- list((sum(np[1:k-1])+1):sum(np[1:k]))
    else  ind.par[(sum(rep.par[1:k-1])+1):sum(rep.par[1:k])] <- list((sum(np[1:k-1])+1):sum(np[1:k]))
	}
  }
  # On garde la même liste de vecteurs de locus en entrée et en sortie
  x.loc.out <- x.loc.in
  }
else
  {
  # Vérification des dimensions de la matrice de contraintes et du vecteur de paramètres contraints
  if (nrow(constraints) != K-1 | nrow(par.constrained) != K-1) stop("The constraint matrix does not have K-1 = ",K-1," rows.")
  if (ncol(par.constrained) != ncol(constraints)) stop ("The number of columns of par.constrained (",ncol(par.constrained),") is not equal to the number of constraints (",ncol(constraints),").")
  # On s'assure que chaque paramètre est impliqué dans au plus une contrainte
  # On boucle sur chaque vecteur de paramètre impliqué dans au moins une contrainte
  # Attention! La vérification est faite uniquement pour les indices de paramètres s'appliquant à la catégorie 1
  # Il faudra s'assurer qu'on n'a pas besoin de le faire pour d'autres
  for (i in unique(par.constrained))
    {
	# Identification des indices dans lequel un élément du vecteur de paramètres est impliqué
	constr.indices <- apply(par.constrained,1,function (vec,i) which(vec==i),i=i)
	# Calcul du nombre de contraintes dans lequel un élément du vecteur de paramètres est impliqué
	if (is.list(constr.indices))
	{
	for (k in 1:(K-1))
	  {
	  if (length(constr.indices[[k]]) > 1)
	    {
	    nconstr.element <- apply(constraints[,constr.indices[[k]]],1,function(vec) sum(vec!=0))
		#print(nconstr.element)
	    if (any(nconstr.element > 1)) stop("Element(s) ",which(nconstr.element > 1)," of parameter vector ",k," is involved in more than one constraint.")
	    }
	  }
	}
	}
  # Calcul du nombre de colonnes de xp, la matrice de devis polytomique contraint
  nc <- ncol(constraints)
  # Le nombre de paramètres contraint est égal à un de moins que les paramètres impliqués dans la contrainte.
  npar.contraints <- apply(constraints,2,function(vec) sum(vec!=0)-1)
  # Nombre total de paramètres libres
  ncolxp <- sum(np) - sum(npar.contraints)
  xp <- array(0,c(nrow(x[[1]]),ncolxp,K-1))
  # Remplissage de xp
  # Pour la catégorie 1, on prend la matrice x, mais on multiplie s'il y a lieu par le coefficient de la contrainte
  xp[,1:np[1],1] <- x[[1]]
  # On identifie les coefficients impliqués dans des contraintes avec les catégories 2 à K-1
  constr.ind.avenir <- which(constraints[1,] != 0)
  for (i in constr.ind.avenir) xp[,par.constrained[1,i],1] = xp[,par.constrained[1,i],1] * constraints[1,i]
  
  # La liste de locus impliqués dans les colonnes de xp pour la catégorie 1 est la liste de locus en entrée
  x.loc.out <- vector("list",length=ncolxp)
  x.loc.out[1:np[1]] <- x.loc.in[1:np[1]]
  
  # La liste des indices où sont situés les paramètres pour la catégorie 1 est 1 à np[1]
  ind.par <- vector("list",length=ncolxp)
  ind.par[1:np[1]] <- list(1:np[1])
  
  if (K > 2)
  {
  # Nombre cumulatif de paramètres libres après chaque catégorie
  nparcumul <- c(np[1],rep(NA,K-2))
  # Pour les catégories de 2 à K-1
  for (k in 2:(K-1))
    {
#	cat("k = ",k)
	# On identifie les coefficients impliqués dans des contraintes avec les catégories 1 à k-1. Note: correction par Jordie le 23 mai 2012 (remplacement du "as.matrix" par un "array")
    constr.ind <- which(constraints[k,] != 0 & apply(array(constraints[1:(k-1),],c(k-1,ncol(constraints))),2,function (vec) any(vec!= 0)))
	if (length(constr.ind)>0)
	  {
	  npk <- np[k] - length(constr.ind)
	  # On alloue un vecteur de longueur égale au nombre de paramètres de la catégorie
	  # Pour l'instant, on le fait seulement s'il y a au moins un nouveau paramètre.
	  # Il faudra réfléchir à ce qu'on fait s'il n'y a pas de nouveau paramètres
	  if(npk>0)	  
	    {
		for (h in 1:npk) ind.par[[nparcumul[k-1]+h]] <- numeric(np[k])
		}
	  # On copie les valeurs de x pour les coefficients impliqués dans des contraintes avec les catégories 1 à k-1
      for (i in constr.ind) 
	    {
		# Première catégorie impliquée dans la contrainte avec le paramètre courant
		cat.constr <- which (constraints[,i] != 0)[1]
		# Nombre de paramètres qu'il faut sauter pour arriver à ceux de la catégorie impliquée dans la contrainte
		decalage <- ifelse(cat.constr>1,nparcumul[cat.constr-1],0)
        xp[,decalage + par.constrained[k,i],k] = x[[k]][,par.constrained[k,i]] * constraints[k,i]
		# On ajoute les locus associés aux paramètres contraints s'ils ne sont pas déjà dans la liste (unique() élimine alors les redondances)
		x.loc.out[[decalage + par.constrained[k,i]]] = unique(c(x.loc.out[[decalage + par.constrained[k,i]]],x.loc.in[(sum(np[1:k-1])+1):sum(np[1:k])][[par.constrained[k,i]]]))
        # On ajoute les indices des paramètres contraints
	    if(npk>0)	  
	      {
		  for (h in 1:npk) ind.par[[nparcumul[k-1]+h]][which(constr.ind==i)] <- decalage + par.constrained[k,i]
		  }
#		  print (decalage + par.constrained[k,i])
#		  print (ind.par)
	    }
	  # On ajoute les valeurs de x restantes
      if(npk>0)
      {
      	  xp[,(nparcumul[k-1]+1):(nparcumul[k-1]+npk),k] = x[[k]][,(1:np[k])[-par.constrained[k,constr.ind]]]
	    # On ajoute à x.loc.out les listes de locus pour les nouveaux paramètres
	    x.loc.out[(nparcumul[k-1]+1):(nparcumul[k-1]+npk)] <- x.loc.in[(sum(np[1:k-1])+1):sum(np[1:k])][-par.constrained[k,constr.ind]]
	    # On ajoute à ind.par les indices des nouveaux paramètres 
	    for (h in 1:npk) ind.par[[nparcumul[k-1]+h]][(length(constr.ind)+1):np[k]] <- (nparcumul[k-1]+1):(nparcumul[k-1]+npk)
	    }
	  }
	else
	# Aucun coefficient n'est impliqué dans des contraintes avec les catégories 1 à k-1
	# On ajoute les nouveaux coefficients à la suite
	  {
	  npk <- np[k]
	  xp[,(nparcumul[k-1]+1):(nparcumul[k-1]+npk),k] = x[[k]]
	  x.loc.out[(nparcumul[k-1]+1):(nparcumul[k-1]+npk)] <- x.loc.in[(sum(np[1:k-1])+1):sum(np[1:k])]
	  for (j in (nparcumul[k-1]+1):(nparcumul[k-1]+npk))
	    ind.par[[j]] <- (nparcumul[k-1]+1):(nparcumul[k-1]+npk)
	  }
	if (k < K-1)
	  {
	  # On identifie les coefficients impliqués dans des contraintes avec les catégories k+1 à K-1. Note: correction par Jordie le 23 mai 2012 (remplacement du "as.matrix" par un "array")
      constr.ind.avenir <- which(constraints[k,] != 0 & apply(array(constraints[1:(k-1),],c(k-1,ncol(constraints))),2,function (vec) all(vec== 0)))
	  
      # On multiplie s'il y a lieu par le coefficient de la contrainte
	  for (i in constr.ind.avenir)
	    {
	    xp[,nparcumul[k-1] + par.constrained[k,i],k] = xp[,nparcumul[k-1]  + par.constrained[k,i],k] * constraints[k,i]
	    }
      }
	nparcumul[k] <- nparcumul[k-1] + npk
    }
	}
  }
  list(xp=xp,x.loc.out=x.loc.out,ind.par=ind.par)
}