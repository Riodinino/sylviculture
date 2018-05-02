rm(list = ls()); invisible(gc())

path <- "C:/Users/nino.page/Desktop/TROLL project/data/traits"
path %>% setwd

library(tidyverse)

dat <- read.csv(file = "BRIDGE_data_ind_2013_10_29.csv", header = T, sep = ";", dec = ".")
dat <- dat %>% mutate(species = paste(Genus, species, sep = "_"))
dat <- dat %>% mutate(LMA = 10000*dry_mass/traits_surf_area)
dat <- dat %>% mutate(wsg = sapwood_dens)

library(mice)

### Fonction qui prend en compte les deux tableaux de traits Bridge et d'histoire de vie (Hmax,WD et seed mass)
## Le dernier argument InvSp est le tableau qui associe Famille, Genre, espèce pour Paracou 
#(et donc potentiellement des espèces qui ne sont pas dans les bases de données)

Traits_filling<-function(Traits,Seltraits<-c("LMA","N","P","wsg","dry_mass","C","thoughness","thickness","Chlorophyll_content")){ #invSp
  

  #Elimination des lignes vides 
  Traits<-Traits[which(apply(Traits,1,function(li){return(any(!is.na(li)))})),]
  
  #Selecion des traits qu'on va garder, pour les indices de diversité
  #"S_mass")
  traits<-traits[,c("Family","Genus","name","bar_code",Seltraits)]
  # maybe thickness, thoughness
  
  ### espèces dont seulement certains traits manquent
  
  # On prend ligne par ligne (fonction apply) et on renvoie les numéro de lignes où au moins un élément est différent de NA
  traitsPartial<-traits[which(apply(traits[,Seltraits],1,function(li){return(any(!is.na(li)))})),]
  # On vérifie juste qu'on a bien aucune ligne où il n'y a que des NA... des fois que!
  traitsPartial<-traitsPartial[which(apply(traitsPartial[,Seltraits],1,function(li){return(any(is.na(li)))})),]
  
  # lignes où tous les traits sont renseignés. On prend ligne par ligne (fonction apply) et 
  #on renvoie les numéros de colonnes aucun élément n'est NA
  
  traitsComp<-traits[which(apply(traits[,Seltraits],1,function(li){return(!any(is.na(li)))})),]
  
  # espèces où il n'y a aucun élément, ligne par ligne, on prend celles ou aucun élément n'est différent de NA
  # En pratique, il n'y en a pas, mais on ne sait jamais!
  traitsEmpty<-traits[which(apply(traits[,Seltraits],1,function(li){return(!any(!is.na(li)))})),]
  
  # on récupère les espèce squi ne sont pas dans les bases de données mais inventoriées à Paracou, 
  # on les ajoute aux espèces sans aucun trait

  ##invSp<-InvSp[which(InvSp[,"name"]%in%setdiff(InvSp[,"name"],unique(traits[,"name"]))),]
  ##invSp<-cbind(invSp,matrix(NA,nrow=nrow(invSp),ncol=ncol(traits)-ncol(invSp)))
  ##colnames(invSp)<-colnames(traits)
  ##traitsEmpty<-rbind(traitsEmpty,invSp)
  
  # Pour les espèces aux traits partiellement renseignés, on découpe en une liste par genre
  traitsmix <- rbind(traitsPartial,traitsComp)
  traitsPartial_gen<-lapply(unique(traitsmix[,"Genus"]),function(gen){
    return(traitsPartial[which(traitsmix[,"Genus"]==gen),])})
  
  # On vérifie que pour chaque genre on a suffisemment d'individus pour faire le gap filling, 
  # on ne garde que les genre où il y a au moins 20 échantillons (purement arbitraire)
  traitsPartial_fam<-do.call(rbind,traitsPartial_gen[which(unlist(lapply(traitsPartial_gen,nrow))<=20)])
  
  # Si on a moins de 20 individus, on estimera avec les individus de la même famille
  traitsPartial_gen<-traitsPartial_gen[which(unlist(lapply(traitsPartial_gen,nrow))>20)]
  
  # gap filling à proporement parler, à chaque élément de la liste de sgenre, on applique la fonction mice
  traitsPartial_gen<-lapply(traitsPartial_gen,function(sub){
    ret<-complete(mice(sub[,Seltraits],printFlag=F))
    return(cbind(sub[,c("Family","Genus","name","bar_code")],ret))})
  
  # J'avais de soucis avec les NA persistants pour la Wood Density, 
  #si le problème persiste je ramène ces espèces à la famille et je refais l'estimation
  ##traitsPartial_fam<-rbind(traitsPartial_fam,
  ##                         do.call(rbind,traitsPartial_gen[which(unlist(lapply(traitsPartial_gen,function(sub){anyNA(sub[,"WD"])})))]))
  ##traitsPartial_gen<-do.call(rbind,
  ##                           traitsPartial_gen[which(!unlist(lapply(traitsPartial_gen,function(sub){anyNA(sub[,"WD"])})))])
  
  # Idem, soucis avec les Arecaceae, pas de données
  traitsPartial_fam<-lapply(unique(traitsPartial_fam[,"Family"]),function(fam){
    return(traitsPartial_fam[which(traitsPartial_fam[,"Family"]==fam & traitsPartial_fam[,"Family"]!="Arecaceae"),])})
  traitsPartial_all<-do.call(rbind,traitsPartial_fam[which(unlist(lapply(traitsPartial_fam,nrow))<=20)])
  traitsPartial_all<-rbind(traitsPartial_all,traitsPartial[which(traitsPartial[,"Family"]=="Arecaceae"),])
  traitsPartial_fam<-traitsPartial_fam[which(unlist(lapply(traitsPartial_fam,nrow))>20)]
  
  # Gap filling des groupes par famille
  traitsPartial_fam<-lapply(traitsPartial_fam,function(sub){
    ret<-complete(mice(sub[,Seltraits],printFlag=F))
    return(cbind(sub[,c("Family","Genus","name","bar_code")],ret))})
  
  # Pour les espèce squi n'avaient toujours pas assez d'individus dans leur famille pour permettre l'estimation
  idScarce<-as.character(traitsPartial_all[,"bar_code"],
                         do.call(rbind,traitsPartial_fam[which(unlist(
                           lapply(traitsPartial_fam,function(sub){anyNA(sub[,"WD"])})))])[,"bar_code"])
  # Pour ces individus on utilise l'ensemble du jeu de données pour estimer les valeurs manquantes
  comScarce<-cbind(traits[,c("Family","Genus","name","bar_code")],
                   complete(mice(traits[,Seltraits],printFlag=F)))
  comScarce<-comScarce[which(comScarce[,"bar_code"]%in%idScarce),]
  ##traitsPartial_fam<-do.call(rbind,
  ##                           traitsPartial_fam[which(!unlist(lapply(traitsPartial_fam,function(sub){anyNA(sub[,"WD"])})))])
  
  # On rassemble tout
  traitsPartial_filled<-rbind(traitsPartial_gen,traitsPartial_fam,comScarce)
  
  # Pour les espèces dont on a aucune valeur de traits, 
  # on tire un ensemble de traits complet correspondant à l'une des espcèces de l'inventaire.
  # on tire au hasard une espèce dans le même genre, ou dans la même famille, 
  # ou dans tout le jeu de données si l'espèce est seule dans sa famille
  ##if(nrow(traitsEmpty)!=0){
  ##  traitsEmpty_filled<-do.call(rbind,lapply(1:nrow(traitsEmpty),function(li){
  ##    sp<-traitsEmpty[li,]
  ##    Tosample<-traitsComp[which(traitsComp[,"name"]==as.character(sp["name"])),Seltraits]
  ##    if(nrow(Tosample)!=0){sp[Seltraits]<-Tosample[sample(rownames(Tosample),1),]}
  ##    
  ##    if(nrow(Tosample)==0){
  ##      
  ##      Tosample<-traitsComp[which(traitsComp[,"Genus"]==as.character(sp["Genus"])),Seltraits]
  ##      if(nrow(Tosample)!=0){sp[Seltraits]<-Tosample[sample(rownames(Tosample),1),]}
  ##      
  ##      if(nrow(Tosample)==0){
  ##        Tosample<-traitsComp[which(traitsComp[,"Family"]==as.character(sp["Family"])),Seltraits]
  ##        
  ##        if(nrow(Tosample)!=0){sp[Seltraits]<-Tosample[sample(rownames(Tosample),1),]}
  ##        if(nrow(Tosample)==0){sp[Seltraits]<-traitsComp[sample(rownames(traitsComp),1),Seltraits]}
  ##      }}
  ##    return(sp)
  ##  }))
  ##  traits_filled<-rbind(traitsEmpty_filled,traitsPartial_filled,traitsComp)
  ##}
  
  ##if(nrow(traitsEmpty)==0){
    ##traits_filled<-rbind(traitsPartial_filled,traitsComp)
  ##}
  traits_filled <- -traitsPartial_filled
  
  # Quelques derniers ajustements, on a en sortie un tableau individus x traits complet
  traits_filled<-traits_filled[order(sort(traits_filled[,"name"])),which(colnames(traits_filled)!="moisture")]
  ##traits_filled<-traits_filled[which(!is.na(traits_filled[,"Hmax"])),]
  ##traits_filled<-traits_filled[which(!traits_filled[,"Genus"]=="Astrocaryum"),]
  ##traits_filled<-traits_filled[which(traits_filled[,"SLA"]>=10),]
  return(traits_filled)
}
