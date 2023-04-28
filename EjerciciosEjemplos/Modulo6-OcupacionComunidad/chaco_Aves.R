##########################################################
##### CURSO Modelado y estimaci?n de ocupaci?n para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########       Ejemplo Ocupaci?n en JAGS      #########
###########              Comunidades             #########
# ##########################################################
# Macchi, L., Decarre, J., Goijman, A. P., Mastrangelo, M., 
# Blendinger, P. G., Gavier‐Pizarro, G. I., ... & Kuemmerle,
#T. (2020). Trade‐offs between biodiversity and agriculture
#are moving targets in dynamic landscapes. Journal of 
#Applied Ecology, 57(10), 2054-2063.
############################################################

#################################
rm(list=ls(all=TRUE))

###set working directory
setwd("C:/Users/andrea/Documents/GitHub/Curso-Ocupacion23/EjerciciosEjemplos/Modulo6-OcupacionComunidad")  #Andrea

library(jagsUI)    #paquete JAGS
library(reshape)
library(plyr)

data<-read.csv("species_sitiosAves.csv",header = T)   #datatest.csv

###Reshape the data using the R package "reshape"
####The detection/non-detection data is reshaped into a three dimensional 
#array "y" where the first dimension, j, is the point; the second dimension, k, 
#is the rep; the third dimension, i, is the sp 

data.melt=melt(data,id.var=c("species", "site", "point"), measure.var="n.ind")
y1=cast(data.melt, site ~ point ~ species)
dim(y1)

y1[is.na(y1)] <- 0 

### All observations >1 <-1, otherwise 0 (this changes numbers to only 0 or 1)
y<-aaply(y1,1,function(x) {x[x>1]<-1;x})

#Remove "NA" sp (hay sitios sin registros entonces hay NA)
y<-y[,,-198]  

str(y)
dimnames(y)
dim(y)

# Names of each dimension
species.list<-as.vector(names(y[1,1,]))
site.list<-as.vector(names(y[,1,1]))
rep.list<-as.vector(names(y[1,,1]))

# Dimensions of each element
nspec<-length(species.list)
nsite<-length(site.list)
nrep<-length(rep.list)


##########################################################################
# Covariate data #########################################################
##########################################################################

###Read in the covariate data
covas <- read.table("covas_sitiosAves.csv", header=TRUE,sep=",",na.strings=c("NA"))   #covas.csv

names(covas)

#MODEL 1
#Detectabilidad = source (mmdet+ierdet) + habitat.type (htype) + julian.date

###Standardize covariates and create new column
# covas ocupacion
attach(covas)
m.annual.rain<-mean(annual.rain); sd.annual.rain<-sd(annual.rain)
m.aridity<-mean(aridity); sd.aridity<-sd(aridity)
#m.forest_6km<-mean(forest_6km, na.rm=T); sd.forest_6km<-sd(forest_6km, na.rm=T)
m.forest_10km<-mean(forest_10km, na.rm=T); sd.forest_10km<-sd(forest_10km, na.rm=T)
#m.yieldP<-mean(yieldP); sd.yieldP<-sd(yieldP)
#m.yieldE<-mean(yieldE); sd.yieldE<-sd(yieldE)
m.yieldM<-mean(yieldM); sd.yieldM<-sd(yieldM)

arainz<-covas$arainz <- (annual.rain - m.annual.rain)/sd.annual.rain
aridityz<-covas$aridityz <- (aridity - m.aridity)/sd.aridity
forest10z<-covas$forest10z <- (forest_10km - m.forest_10km)/sd.forest_10km
yieldmz<-covas$yieldmz <- (yieldM - m.yieldM)/sd.yieldM

#covas de deteccion
m.julian.date<-mean(julian.date); sd.julian.date<-sd(julian.date)
jdatez<-covas$jdatez <- (julian.date - m.julian.date)/sd.julian.date
detmm<-covas$mmdet
detier<-covas$ierdet
htype<-covas$htype

detach(covas)

### Vector con numero de puntos por sitio
puntos<-read.table("puntos_por_sitioAves.csv", header=TRUE,sep=",",na.strings=c("NA"))   #vector.puntos
J<-puntos$puntos
length(J)

##########################################################################
##### Model - Random effects community occupancy model####################
######Allows for estimation over less frequent species####################
##########################################################################

# Augment data set (DA part)
nz <- 215                # Number of potential species in superpopulation
nz <- nz-nspec          # Size of augmented data set ('superpopulation')
yaug <- array(0, dim=c(nsite,nrep,nz+nspec)) # Add all zero histories
yaug[,,1:nspec]<-y
dim(yaug)

dim(y)

#site, rep, sp
#sp, site, rep
yaug2<-aperm(yaug,c(3,1,2))
yaug<-yaug2

# Armar paquete de datos y resumirlos
str(chaco.data <- list(y=yaug, forest=forest10z, yield=yieldmz, aridity=aridityz, 
                       detmm=detmm, htype=htype, detier=detier, nsite=nsite, nspec=nspec, 
                       M=nspec+nz, nz=nz, J=J))

### JAGS code
sink("chacomb.jags")  
cat("
    model{
  
    # Priors 
    omega~dunif(0,1)
    
    # Hyperpriors for species effects          
    mu.lpsi ~ dnorm(0, 1/2.25^2)  
    tau.lpsi <- pow(sd.lpsi,-2)   
    sd.lpsi ~ dunif(0,8)              
    mu.betalpsi1~ dnorm(0, 1/2.25^2)
    tau.betalpsi1 <- pow(sd.betalpsi1,-2)
    sd.betalpsi1 ~ dunif(0,8)
    mu.betalpsi2~ dnorm(0, 1/2.25^2)
    tau.betalpsi2 <- pow(sd.betalpsi2,-2)
    sd.betalpsi2 ~ dunif(0,8)
    tau.betalpsi3 <- pow(sd.betalpsi3,-2)
    sd.betalpsi3 ~ dunif(0,8)
    mu.betalpsi4~ dnorm(0, 1/2.25^2)
    tau.betalpsi4 <- pow(sd.betalpsi4,-2)
    sd.betalpsi4 ~ dunif(0,8)
    mu.lp ~ dnorm(0, 1/2.25^2)
    tau.lp <- pow(sd.lp,-2)
    sd.lp ~ dunif(0,5)
    mu.betalp1 ~ dnorm(0, 1/2.25^2)
    tau.betalp1 <- pow(sd.betalp1,-2)
    sd.betalp1 ~ dunif(0,5)
    mu.betalp2 ~ dnorm(0, 0.01)
    tau.betalp2 <- pow(sd.betalp2,-2)
    sd.betalp2 ~ dunif(0,5)
    mu.betalp3 ~ dnorm(0, 1/2.25^2)
    tau.betalp3 <- pow(sd.betalp3,-2)
    sd.betalp3 ~ dunif(0,5)
    
    # priors for species specific effects in occupancy and detection 
    # sp represents per species
    for (k in 1:M){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    lp[k] ~ dnorm(mu.lp, tau.lp)
    betalpsi1[k] ~ dnorm (mu.betalpsi1, tau.betalpsi1)
    betalpsi2[k] ~ dnorm (mu.betalpsi2, tau.betalpsi2)
    betalpsi3[k] ~ dnorm (mu.betalpsi3, tau.betalpsi3)
    betalpsi4[k] ~ dnorm (mu.betalpsi4, tau.betalpsi4)
    betalp1[k] ~ dnorm(mu.betalp1,tau.betalp1)
    betalp2[k] ~ dnorm(mu.betalp2,tau.betalp2)
    betalp3[k] ~ dnorm(mu.betalp3,tau.betalp3)
    
    # Superpopulation process
    w[k]~dbern(omega)
    
    # Ecological model, process model (true occurrence at site i) 
    for (i in 1:nsite) {                                           #loop sobre sitios 
    logit(psi[k,i]) <- lpsi[k] + betalpsi1[k]*forest[i]  
    betalpsi2[k]*yield[i]+ betalpsi3[k]*aridity[i] 
    mu.psi[k,i]<-w[k]*psi[k,i]
    z[i,k] ~ dbern(mu.psi[k,i])     
    
    
    # Observation model for site i, replicate nrep=j, nspec=k
    for (j in 1:J[i]) {                                   #loop sobre unbalanced reps J   
    logit(p[k,i,j]) <-  lp[k] + betalp1[k]* htype[i]+ betalp2[k]* detmm[i]
    + betalp3[k]* detier[i]
    mu.p[k,i,j] <- p[k,i,j]*z[i,k]
    y[k,i,j] ~ dbern(mu.p[k,i,j])
    } #nrep
    } #nsite
    } #nspec
    
    # Derived quantities
    for (i in 1:nsite) {                                           
    Nsite [i] <- sum (z[i,])     # number of sp per site
    }
    for (k in 1:nspec){                                         
    Nocc.fs [k] <- sum (z[,k])   #number of occupied sites
    }
    n0<-sum(w[(nspec+1):(nspec+nz)]) #Number of unseen species
    Ntotal<- sum(w[])                #total metacommunity size
    
    lpsiS[1:(nspec+1)]<- lpsi[1:(nspec+1)]
    betalpsi1S[1:(nspec+1)]<- betalpsi1[1:(nspec+1)]
    betalpsi2S[1:(nspec+1)]<- betalpsi2[1:(nspec+1)]
    betalpsi3S[1:(nspec+1)]<- betalpsi3[1:(nspec+1)]
    betalpsi4S[1:(nspec+1)]<- betalpsi4[1:(nspec+1)]
    lpS[1:(nspec+1)]<- lp[1:(nspec+1)]
    betalp1S[1:(nspec+1)]<- betalp1[1:(nspec+1)]
    betalp2S[1:(nspec+1)]<- betalp2[1:(nspec+1)]
    betalp3S[1:(nspec+1)]<- betalp3[1:(nspec+1)]
    
    } #model
    ",fill=TRUE)
sink()

# valores iniciales
wst <- rep(1,nspec+nz)
zst <- array(1,dim=c(nsite, nspec+nz))

inits <- function() list(z=zst,w=wst,lpsi=rnorm(nspec+nz,0,2.25), betalpsi1=rnorm(nspec+nz,0,2.25),
                         betalpsi2=rnorm(nspec+nz,0,2.25),betalpsi3=rnorm(nspec+nz,0,2.25),betalpsi4=rnorm(nspec+nz,0,2.25),
                         lp=rnorm(nspec+nz,0,2.25),betalp1=rnorm(nspec+nz,0,2.25),betalp2=rnorm(nspec+nz,0,2.25),
                         betalp3=rnorm(nspec+nz,0,2.25))

# parametros a monitorear
## parameters for model selection
par.ms <- c("lpsiS","betalpsi1S","betalpsi2S","betalpsi3S","betalpsi4S","lpS", "betalp1S", "betalp2S", 
            "betalp3S","mu.betalpsi1","mu.betalpsi2", "mu.betalpsi3","mu.betalpsi4","mu.lpsi","Nsite")

# ajustes de MCMC
ni <- 500 #30000 
nt <- 10    
nb <- 100 #1000  
nc <- 3
na <- 100 #10000 


outms = jags(chaco.data, inits, par.ms, "chacomb.jags", n.chains=nc, n.iter=ni, n.burnin=nb, 
             n.thin=nt, n.adapt=na, parallel=TRUE)

# corrida completa 9hs
#save results
#save(outms, file='modelaves.rda') 


# -----------------------------------------------
# Actividad 
# encontrar el error en el modelo
# revisar los graficos y graficar las covariables que falten
# -----------------------------------------------
  
#####################################################################################################################
#####################################################################################################################
# Nueva Sesion, correr la parte inicial, luego load mcmc.objects
load('modelaves.rda') 
str(outms)

# summary of the posteriors
outms

# tabla de salidas
#write.csv(outms$summary, "model_outms.csv")

#### Chequeo convergencia de las cadenas ####
par(mfrow=c(3,3))
traceplot(outms, parameters = NULL)


### Density plots for community effects
d1<-density(outms$sims.list$mu.betalpsi1)
d2<-density(outms$sims.list$mu.betalpsi2)
d3<-density(outms$sims.list$mu.betalpsi3)
d4<-density(outms$sims.list$mu.betalpsi4)

par(mfrow=c(2,2))
plot(d1, main = "Forest 10km on community")
abline(v=0)
plot(d2, main = "Meat yield on community")
abline(v=0)
plot(d3, main = "Aridity on community")
abline(v=0)
plot(d4, main = "Forest-Yield interaction")
abline(v=0)

# Rename
out1<-outms

##COVARIATES
#forest10z
#yieldmz
#aridityz
#detmm
#htype
#detier

## Litado de especies para graficar
sp<-species.list

# OCUPACION
lpsi <- apply(out1$sims.list$lpsi,c(2),mean)
cri.lpsi <- apply(out1$sims.list$lpsi,c(2),function(x) quantile(x,prob=c(0.025,0.975)))
betalpsi1 <- apply(out1$sims.list$betalpsi1,c(2),mean)
cri.betalpsi1 <- apply(out1$sims.list$betalpsi1,c(2),function(x) quantile(x,prob=c(0.025,0.975)))

# DATOS ORDENADOS
# Forest 10km
betalpsi1_all<-cbind(betalpsi1[1:196],cri.betalpsi1[1,1:196],cri.betalpsi1[2,1:196],seq(1:196))
betalpsi1_all<-betalpsi1_all[order(-betalpsi1_all[,1]),]
betalpsi1_all

#######################################################
##### Plot covariate effects on species occupancy #####
#####                   (Betas)                   #####         
#######################################################

par(mfrow=c(1,1), cex=1, cex.lab=1, cex.axis=0.9, cex.main=1)
par(mar=c(2.8,2.5,1.2,1))
par(mgp = c(1.5,0.5,0))

## Yield in 196 species (ORDENADO)
plot(betalpsi2_all[1:196,1],1:196, xlim=c(-2.5,2),xlab='Yield meat', ylab='Species',
     cex=0.8,pch=16,tck=-0.02, main='Covariate effects on species (Beta)')
abline(v=0, lwd=2, col='black')
segments(betalpsi2_all[1:196,2],1:196,betalpsi2_all[1:196,3],1:196,col='grey', lwd=1)
sig2<- betalpsi2_all[,2]*betalpsi2_all[,3]>0    #mayor a cero son los dos + o dos - == efecto signif   
segments(betalpsi2_all[,2][sig2==1],(1:196)[sig2==1],betalpsi2_all[,3][sig2==1],(1:196)[sig2==1],lwd=1, col='red')
abline(v=out1$summary[1784,1],lwd=1,col='red')
abline(v=out1$summary[1784,c(3,7)],lwd=1,col='red',lty=2)

# -----------------------------------------------
# Actividad 
# 
# Graficar efectos del resto de los parametros
# -----------------------------------------------

#################################################
#### DETECTION
#################################################
# posterior means for code convenience
pm<-out1$mean

par(mfrow=c(1,1))
# Trees
detp<-array (NA,dim = c(nspec))
for(i in 1:nspec){
  detp[i]<-plogis(pm$lpS[i]+pm$betalp1S[i])}
barplot(detp, ylim = c(0,1))

# estimar detectabilidad c arboles
detp<-array (NA,dim = c(nspec))
for(i in 1:nspec){
  detp[i]<-plogis(pm$lpS[i]+pm$betalp1S[i])}

# ordenando parametros
det_order<-cbind(detp,seq(1:195))
det_order<-det_order[order(det_order[,1]),]

barplot(det_order[,1], ylim = c(0,0.7),xlab = 'Species',
        ylab = 'Detection probabilities', names.arg = seq(1:195), 
        axisnames = TRUE, main = 'Detection with trees')


###########################################
#####     Community psi & Covariates  #####
########################################### 
## Graficos de relaciones c covas
# Armar vectores de covariables para graficar predicciones
# Forest 10km (0 - 100) 
# YieldM (0,582)
# Aridity (0.56,0.8)

o.forest <- seq(2,75,,500)   # con los valores umbrales
pred.forest <- (o.forest - m.forest_10km)/sd.forest_10km
o.aridity <- seq(0.56,0.80,,500)
pred.aridity <- (o.aridity - m.aridity)/sd.aridity
o.yieldm <- seq(0,582,,500)
pred.yieldm <- (o.yieldm - m.yieldM)/sd.yieldM



##### Community and and interactions
## Simulation results
tmp<-out1$sims.list
nsamp<-out1$mcmc.info$n.samples

## Create array for predicted simulations
pred<- array(NA, dim=c(500,nsamp,3))

### Yield M y Forest (Aridity to the mean = 0)
### Community
for(i in 1:nsamp){
  pred[,i,1]<-plogis(tmp$mu.lpsi[i]+tmp$betalpsi1[i]*(-2.29)+tmp$mu.betalpsi2[i]*pred.yieldm+
                       tmp$mu.betalpsi4[i]*pred.yieldm*(-2.29))                                  #fors min (-2.29)
  pred[,i,2]<-plogis(tmp$mu.lpsi[i]+tmp$betalpsi1[i]*(-0.88)+tmp$mu.betalpsi2[i]*pred.yieldm+
                       tmp$mu.betalpsi4[i]*pred.yieldm*(-0.88))                                 #forsc med (-0.88)
  pred[,i,3]<-plogis(tmp$mu.lpsi[i]+tmp$betalpsi1[i]*(0.52)+tmp$mu.betalpsi2[i]*pred.yieldm+
                       tmp$mu.betalpsi4[i]*pred.yieldm*(0.52))                                   #fors max (0.52)
}

cri<-function(x) quantile(x,prob=c(0.025,0.975))
cri1<-apply(pred[,,1],1,cri)
cri3<-apply(pred[,,3],1,cri)


par(mfrow=c(1,3))
selection<-sample(1:nsamp,20) #500 
matplot(pred.yieldm, pred[,selection,1], ylab="Occupancy mean", xlab="Meat Yield",
        type="l",lty=1,lwd=1,col="grey",ylim=c(0,1),frame=F,main="Community - Min. forest (2%)")
lines(pred.yieldm, apply(pred[,,1],1,mean),lwd=3,col="black")
matplot(pred.yieldm, pred[,selection,2], ylab="Occupancy mean", xlab="Meat Yield",
        type="l",lty=1,lwd=1,col="grey",ylim=c(0,1),frame=F,main="Community - Med. forest")
lines(pred.yieldm, apply(pred[,,2],1,mean),lwd=3,col="black")
matplot(pred.yieldm, pred[,selection,3], ylab="Occupancy mean", xlab="Meat Yield",
        type="l",lty=1,lwd=1,col="grey",ylim=c(0,1),frame=F,main="Community - Max. forest (75%)")
lines(pred.yieldm, apply(pred[,,3],1,mean),lwd=3,col="black")


############################
##### PLOTS DE RIQUEZA ##### 
############################

## Sample of posterior distributions of local species richness (alpha diversity) at all sites                                                                                              
par(mfrow=c(1,1), cex=1, cex.lab=1, cex.axis=0.9, cex.main=1)
plot(table(out1$sims.list$Nsite[,]),ylab="",ylim=c(0,60000),xlab='Number of species',pch=16,
     main="Posterior distribution of total number of species")

## Relationship between Community size (alpha diversity - Derived parameter) and covariates
# Percent Forest
plot(covas$forest_10km,out1$mean$Nsite, xlab="Percent Forest 10km", ylab="Community size estimate", pch=16)


##############################################################################
######## GRAFICOS PARA GRUPOS DE ESPECIES ####################################
##############################################################################

## PLOTS VARIAS SP, NO ERROR
# Forest
occp<-array (NA,dim = c(196,500))
for(i in 1:196){
  occp[i,]<-plogis(pm$lpsiS[i]+pm$betalpsi1S[i]*pred.forest)}

par(mfrow=c(1,1), cex=1, cex.lab=1, cex.axis=0.9, cex.main=1)
par(mar=c(2.8,2.5,1.2,1))
par(mgp = c(1.5,0.5,0))
# Plot Occupancy probabilities with Forest
plot(o.forest,occp[194,], type='l', lty=1, ylab="Occupancy probabilities", main="Edge foliage",
     xlab="Forest",ylim=c(0,1), col=1, lwd=3)
lines(o.forest,occp[91,],col=2, lwd=3)
lines(o.forest,occp[176,],col=3, lwd=3)
lines(o.forest,occp[155,],col=4, lwd=3)
lines(o.forest,occp[107,],col=5, lwd=3)
lines(o.forest,occp[139,],col=6, lwd=3)
lines(o.forest,occp[173,],col=7, lwd=3)
lines(o.forest,occp[189,],col=8, lwd=3)
lines(o.forest,occp[6,],col=9, lwd=3)
lines(o.forest,occp[163,],col=10, lwd=3)
lines(o.forest,occp[47,],col=11, lwd=3)
lines(o.forest,occp[196,],col=12, lwd=3)
lines(o.forest,occp[162,],col=13, lwd=3)
lines(o.forest,occp[144,],col=14, lwd=3)
lines(o.forest,occp[157,],col=15, lwd=3)
lines(o.forest,occp[185,],col=16, lwd=3)
lines(o.forest,occp[96,],col=17, lwd=3)
lines(o.forest,occp[49,],col=18, lwd=3)
lines(o.forest,occp[28,],col=19, lwd=3)


