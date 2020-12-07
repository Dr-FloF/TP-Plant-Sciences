
####################################################################
########## Fonction pour simulation individuelle ###################
###################################################################

Archi <- function(dmin,dmax,igravi,ibd,Dicot,Enzyme,P.Surf,P.Depth,P.Availability,Graphic,Graph.Treshold,Legend)
  
  {
  
  
#############################################
# Sol homogène et non limitant pour commencer 

SOLUNC <- data.frame(Croiss=rep(1.0,times=100),Ramif=rep(1.0,times=100),ICMeca=rep(0.02,times=100),DirectionConstraint=rep(0,times=100)) #création d'un tableau avec pour longueur 100, indice de ramification = 1 ie homogène et parallèle, direction de contrainte mécanique: le sol joue nn seulement sur la vitesse de croissance et sur sa direction (contrainte angle) 
write.table(x=SOLUNC,file="./sol.txt",sep = "\t",row.names = FALSE,col.names = TRUE)


###################################
# Allocation Biomasse non imitante

biom00 <- rep(10,times=300) 

write(x=biom00,file="./biomrac.txt",sep="\n")  

###################################
if (Dicot==TRUE)
{
  PAR<- data.frame(
    duration=60,
    emissratesem=0.9,
    propDiamSem=1,
    nbmaxsem=1,
    ageadvemis=100,
    distmaxadv=80,
    emissrateadv=0.1,
    propDiamadv=1.0,
    nbmaxadv=5,
    dmin=dmin,
    dmax=dmax,
    elongvsdiam=22,
    tgravi=1,
    igravi=igravi,
    durdp=5.0,
    ibd=ibd,
    dldm=0.30,
    vard=0.20,
    tmd=0.10,
    growthdur=600,
    lifeexp=2500,
    radgrowth=0.30 )
} else {
  PAR<- data.frame(
    duration=50,
    emissratesem=0.9,
    propDiamSem=0.65,
    nbmaxsem=3,
    ageadvemis=3,
    distmaxadv=4,
    emissrateadv=5,
    propDiamadv=1.0,
    nbmaxadv=35,
    dmin=dmin,
    dmax=dmax,
    elongvsdiam=22,
    tgravi=1,
    igravi=igravi,
    durdp=5.0,
    ibd=ibd,
    dldm=0.10,
    vard=0.20,
    tmd=0.10,
    growthdur=250,
    lifeexp=2500,
    radgrowth=0 )
}



# Ecriture du fichier des paramètres
write.table(x=PAR,file="paramarch9.txt",sep = "\n",row.names = FALSE,col.names = FALSE)

# Simulation
system("archisimp9.exe 0.0 0.0 5.0") # position de la semence (x, y, z) enterre un peu à 10 mm

# Lecture du fichier de la structure
SR <- read.table(file="seg.txt",header=T)

# Calcul des paramètres racines 
SR$depth <- cut(0.5*(SR$Z1+SR$Z2),breaks=seq(from=0,to=1000,by=50)) # on fait des classes
PROF <- data.frame(depth=seq(from=0,to=950,by=50))
PROF$length <- tapply(sqrt(((SR$X1-SR$X2)^2)+((SR$Y1-SR$Y2)^2)+((SR$Z1-SR$Z2)^2)),SR$depth,sum) 
PROF$surface <- tapply(pi*SR$Diam*sqrt(((SR$X1-SR$X2)^2)+((SR$Y1-SR$Y2)^2)+((SR$Z1-SR$Z2)^2)),SR$depth,sum) # mm²
PROF$Prop.length <- (PROF$length/sum(PROF$length,na.rm=TRUE))*100
PROF$Prop.surface <- (PROF$surface/sum(PROF$surface,na.rm=TRUE))*100

#### Module sol

SOL <- data.frame(depth=seq(from=0,to=950,by=50))

Vol<-NULL

X<-cbind(tapply(SR$X1,SR$depth,min),tapply(SR$X2,SR$depth,min),tapply(SR$X1,SR$depth,max),tapply(SR$X2,SR$depth,max))
Y<-cbind(tapply(SR$Y1,SR$depth,min),tapply(SR$Y2,SR$depth,min),tapply(SR$Y1,SR$depth,max),tapply(SR$Y2,SR$depth,max))
for (i in 1:nrow(X))
{
Vol[i]<-((max(X[i,3:4])-min(X[i,1:2]))*(max(Y[i,3:4])-min(Y[i,1:2]))*50)/1000
}
SOL$Vol.Explo<-Vol # volume exploré en cm3
SOL$Masse.Explo<-(SOL$Vol.Explo*1.3)/1000

SOL$C.P.Total<-ifelse (SOL$depth<200,P.Surf,P.Depth) #  mg.kg

SOL$P.Total.Explo<-SOL$C.P.Total*SOL$Masse.Explo  # mg accessible 

if (Enzyme==TRUE)
{SOL$P.Available<-ifelse((SOL$P.Total.Explo*P.Availability*1.3)<SOL$P.Total.Explo,SOL$P.Total.Explo*P.Availability*1.3,SOL$P.Total.Explo)
}else{
SOL$P.Available<-SOL$P.Total.Explo*P.Availability}

SOL$P.Export<-ifelse(SOL$P.Available*(PROF$surface/100)*0.0005<SOL$P.Available,SOL$P.Available*(PROF$surface/100)*0.0005,SOL$P.Available) # fraction de P accessible autour de la racine

### Summary

Summary<-data.frame(matrix(nrow=1,ncol=11))
names(Summary)<-c("Length","Surface","Shallow.Length","Shallow.Surface","P.Export","dmin","dmax","igravi","ibd","Dicot","Enzyme")
Summary$Length<-sum(PROF$length,na.rm=TRUE)
Summary$Surface<-sum(PROF$surface,na.rm=TRUE)
Summary$Shallow.Length<-sum(PROF$Prop.length[1:4],na.rm=TRUE)
Summary$Shallow.Surface<-sum(PROF$Prop.surface[1:4],na.rm=TRUE)
Summary$P.Export<-sum(SOL$P.Export,na.rm=TRUE)
Summary$igravi<-igravi
Summary$dmin<-dmin
Summary$dmax<-dmax
Summary$igravi<-igravi
Summary$ibd<-ibd
Summary$Dicot<-Dicot
Summary$Enzyme<-Enzyme

### Sorties


ifelse(Graphic==FALSE,
result<-list(Summary=Summary, Plant=PROF , SOL=SOL),
result<-list(XZ=DessineSRXZ(SR[SR$Diam>Graph.Treshold,]),
             if(Legend==TRUE) 
                {leg=legend("bottomleft",paste(
               "dmin=",dmin," ",
               "dmax=",dmax," ",
               "igravi=",igravi," ",
               "ibd=",ibd," ",
               "Dicot=",Dicot),
               cex=0.8)},
             Summary=Summary, 
             Layers=PROF, 
             SOL=SOL))

return(result)

}




########################################################
########## Fonction pour sensibilité ###################
########################################################


Sensitivity.Archi<-function(Var,dmin,dmax,igravi,ibd,Dicot,Enzyme,P.Surf,P.Depth,P.Availability)
{
  Values<-NULL
  sum1<-NULL
  
## sensibilité dmin
  if (Var=="dmin")
  {Values<-dmin
  for (dmin in Values)
{
  simul<-Archi(
    #### Crop characteristics
    dmin=dmin,
    dmax=dmax,
    igravi=igravi,
    ibd=ibd,
    Dicot=Dicot,
    Enzyme=Enzyme,
    
    ### Soil characteristics
    P.Surf = P.Surf,
    P.Depth = P.Depth,
    P.Availability = P.Availability, 
    
    
    ### Type of analysis 
    Graphic=FALSE)
  
  sum1<-rbind(sum1,simul$Summary)
  }
  }
  
## sensibilité dmax
  if (Var=="dmax")
  {Values<-dmax
  for (dmax in Values)
  {
    simul<-Archi(
      #### Crop characteristics
      dmin=dmin,
      dmax=dmax,
      igravi=igravi,
      ibd=ibd,
      Dicot=Dicot,
      Enzyme=Enzyme,
      
      ### Soil characteristics
      P.Surf = P.Surf,
      P.Depth = P.Depth,
      P.Availability = P.Availability, 
 
      ### Type of analysis 
      Graphic=FALSE)
    
    sum1<-rbind(sum1,simul$Summary)
  }
  }
  
## sensibilité igravi
  if (Var=="igravi")
  {Values<-igravi
  for (igravi in Values)
  {
    simul<-Archi(
      #### Crop characteristics
      dmin=dmin,
      dmax=dmax,
      igravi=igravi,
      ibd=ibd,
      Dicot=Dicot,
      Enzyme=Enzyme,
      
      ### Soil characteristics
      P.Surf = P.Surf,
      P.Depth= P.Depth,
      P.Availability = P.Availability, 
      
      ### Type of analysis 
      Graphic=FALSE)
    
    sum1<-rbind(sum1,simul$Summary)
  }
  }

  
## sensibilité ibd
  if (Var=="ibd")
  {Values<-ibd
  for (ibd in Values)
  {
    simul<-Archi(
      #### Crop characteristics
      dmin=dmin,
      dmax=dmax,
      igravi=igravi,
      ibd=ibd,
      Dicot=Dicot,
      Enzyme=Enzyme,
      
      ### Soil characteristics
      P.Surf = P.Surf,
      P.Depth= P.Depth,
      P.Availability = P.Availability, 
      
      ### Type of analysis 
      Graphic=FALSE)
    
    sum1<-rbind(sum1,simul$Summary)
  }
  }
  
  for (k in 1:9)
  {
    sum1["cv",k]<-round(sd(sum1[,k],na.rm=TRUE)/mean(sum1[,k],na.rm=TRUE),2)
  }
  return(Summary.Sensitivity=sum1)
}



