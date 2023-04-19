#comn los dos tiempos en el array de datos
#variables escaladas 
#escalas 150 200 300 400 500 600  - con % de tierra cultivada
#SIN ambientes, solo el %tC

############################################################
rm(list=ls(all=TRUE))

setwd("/home/vanesa/Documentos/Vanesa UNRC/Publicaciones/Paper complejidad/analisis/Analisis revision landscape 2")

library(jagsUI)    #paquete JAGS
library(boot)
library(matrixStats)

data<-read.csv("datos7.csv",header = T)
str(data)
head(data)

attach(data)

#para poder ordenar los datos en la matriz es importante que est?n todos los sitios 
#por mas que no hayan capturado nada y todas las potenciales sp por sitios
#total despues se sacan del an?lisis, es para hacer la matrix 4d mas facil

#tienen que estar acomodados por sesion y despues por sp
data.season.1<-data[1:450,1:7]
data.season.2<-data[451:900,1:7]

species.list <- list("Aa","Ad","Cl","Cm","Cv","Md","Of","Or","Tp")   # alphabetic list
season.list<- list("1", "2")
site.list<- list("1", "2","3","4","5","6","7","8","9","10",
                 "11","12","13","14","15","16","17","18","19","20",
                 "21","22","23","24","25","26","27","28","29","30",
                 "31","32","33","34","35","36","37","38","39","40",
                 "41","42","43","44","45","46","47","48","49","50")
rep.list<- list("1", "2","3","4")

counts.season.1 <- cbind(data.season.1$X1, data.season.1$X2, data.season.1$X3,data.season.1$X4)  # Counts 2014
DET.season.1 <- counts.season.1
DET.season.1[DET.season.1 > 1] <- 1       # now turned into detection/nondetection data

counts.season.2 <- cbind(data.season.2$X1, data.season.2$X2, data.season.2$X3,data.season.2$X4)  # Counts 2014
DET.season.2 <- counts.season.2
DET.season.2[DET.season.2 > 1] <- 1       # now turned into detection/nondetection data

nsite <- 50                    # numero de sitios
nrep <- 4       # number of replicate surveys per season  - cada una de las 4 noches
count<-1
nspec <- length(species.list)   # 9 species en data.1
nseason<-2

Y.season.1 <- array(NA, dim = c(nsite, nrep, nspec))

dimnames(Y.season.1) <- list(site=site.list,rep= rep.list, especie=species.list)
dim(Y.season.1)
Y.season.1

#es un loop para cada conjunto de datos
#lo eque hace es que le digo que me rellene para caada una de las especies con las filas que 
#correspondan 
#para probar probar con i=1

for(i in 1:nspec){
  Y.season.1[,,i] <- DET.season.1[((i-1)*nsite+1):(i*nsite),]     #i es de la especie asi te rellena para cada sp con 50 sitios
}

# for(i in 1:nspec){
#   J1[,,i] <- j1[((i-1)*nsite+1):(i*nsite)]     #i es de la especie asi te rellena para cada sp con 50 sitios
# }

Y.season.2 <- array(NA, dim = c(nsite, nrep, nspec))
dimnames(Y.season.2) <- list(site=site.list,rep= rep.list, especie=species.list)
dim(Y.season.2)
Y.season.2

# J2<-array(NA, dim = c(nsite,count, nspec))
# dimnames(J2) <- list(site=site.list,count="count", especie=species.list)

for(i in 1:nspec){
  Y.season.2[,,i] <- DET.season.2[((i-1)*nsite+1):(i*nsite),]     #i es de la especie asi te rellena para cada sp con 50 sitios
}

Y.season.2

# for(i in 1:nspec){
#   J2[,,i] <- j2[((i-1)*nsite+1):(i*nsite)]     #i es de la especie asi te rellena para cada sp con 50 sitios
# }
# 
# J2
# 
library(abind) #necesaria para combinar las matrices

y<-abind(Y.season.1,Y.season.2,along = 4)  #es 4 porque es en una dimensi?n mas de las que ya tienen (3)

dimnames(y) <- list(site=site.list,rep= rep.list, especie=species.list,season=season.list)  #nombres de las dimensiones
#dim(y[,,5,1])

#Armar paquete de datos y resumirlos
#ysum<- apply(y, c(1,3,4),sum, na.rm=T)
y
# ysumj<- apply(y, c(1,4),length(unique(y)), na.rm=T)
# 


# j1<- data[1:450,8]
# j2<- data[451:900,8]
# Put detection data into 3D array: site x rep x species
# 
# j<-abind(J1,J2,along = 4)  #es 4 porque es en una dimensi?n mas de las que ya tienen (3)
# nrep<- j

# ##### cargo y transformo en matriz las covariables
is150 <- (scale(IS150))
is200 <- (scale(IS200))
is300 <- (scale(IS300))
is400 <- (scale(IS400))
is500 <- (scale(IS500))
is600 <- (scale(IS600))
pa150 <- (scale(pa150))
pa200 <- (scale(pa200))
pa300 <- (scale(pa300))
pa400 <- (scale(pa400))
pa500 <- (scale(pa500))
pa600 <- (scale(pa600))

######NOOO
# obtener en una tabla todos los valores estandarizados y los no estandarizados
# para todas las escalas de estidio para los 50 sitios
# a<-sort(IS150[1:50])
# a1<-round(sort(is150[1:50]),3)
# b<-sort(IS200[1:50])
# b1<-round(sort(is200[1:50]),3)
# c<-sort(IS300[1:50])
# c1<-round(sort(is300[1:50]),3)
# d<-sort(IS400[1:50])
# d1<-round(sort(is400[1:50]),3)
# e<-sort(IS500[1:50])
# e1<-round(sort(is500[1:50]),3)
# f<-sort(IS600[1:50])
# f1<-round(sort(is600[1:50]),3)
# 
# is <-data.frame(a,a1, b,b1,c,c1,d,d1,e,e1,f,f1 )
# colnames(is)<-c("150", "150 stand","200" ,"200 stand","300" ,"300 stand",
#                 "400" ,"400 stand","500" ,"500 stand","600" ,"600 stand")
# # ismean<-apply(is,2,mean)
# # boxplot(is)

#write.csv(is, "valores is.csv")


###Las covariables que quedan son al150, al300, al600, cultivos
### Build individual matrixes (site,year) for each covariate
is150<- array(is150, dim = c(nsite))
dimnames(is150) <- list(site=site.list)  #nombres de las dimensiones
is200<- array(is200, dim = c(nsite))
dimnames(is200) <- list(site=site.list)  #nombres de las dimensiones
is300<- array(is300, dim = c(nsite))
dimnames(is300) <- list(site=site.list)  #nombres de las dimensiones
is400<- array(is400, dim = c(nsite))
dimnames(is400) <- list(site=site.list)  #nombres de las dimensiones
is500<- array(is500, dim = c(nsite))
dimnames(is500) <- list(site=site.list)  #nombres de las dimensiones
is600<- array(is600, dim = c(nsite))
dimnames(is600) <- list(site=site.list)  #nombres de las dimensiones

pa150<- array(pa150, dim = c(nsite))
dimnames(pa150) <- list(site=site.list)  #nombres de las dimensiones
pa200<- array(pa200, dim = c(nsite))
dimnames(pa200) <- list(site=site.list)  #nombres de las dimensiones
pa300<- array(pa300, dim = c(nsite))
dimnames(pa300) <- list(site=site.list)  #nombres de las dimensiones
pa400<- array(pa400, dim = c(nsite))
dimnames(pa300) <- list(site=site.list)  #nombres de las dimensiones
pa500<- array(pa500, dim = c(nsite))
dimnames(pa500) <- list(site=site.list)  #nombres de las dimensiones
pa600<- array(pa600, dim = c(nsite))
dimnames(pa600) <- list(site=site.list)  #nombres de las dimensiones

J1<- data$J[1:nsite]
J2<- data$J[451:500]

J<-abind(matrix(J1),matrix(J2))

  str(roed.data.pa <- list(y=y,J=J, nsite=dim(y)[1],nspec=dim(y)[3],nrep=dim(y)[2],
                      nseason=dim(y)[4], pa150=pa150, pa200=pa200, pa300=pa300,pa400=pa400, 
                      pa500=pa500, pa600=pa600))
###Modelos solo PA
### JAGS code
sink("mod150_pa_r1.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~   dnorm(0, 0.1)
    mu.lp ~   dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)

    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    
    # priors for species-specific effects
    for (k in 1:nspec){
      lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
      b1[k] ~ dnorm (mu.b1, tau.b1)
    
    for (j in 1:nrep) {  
      lp[j,k] ~ dnorm (mu.lp, tau.lp)    
      
    } #nrep
  
  # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*pa150[i]
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    for (l in 1:nseason){                                         #loop sobre tiempo
    for (j in 1:nrep) {  
    logit(p[j,k,l]) <- lp[j,k]     #poner por fuera de i
    } #season
    }

    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (j in 1:nrep) {                                           #loop sobre noche 
    for (l in 1:nseason){                                         #loop sobre tiempo 
      mu.p[i,j,k,l] <- p[j,k,l]*z[i,k,l]
      y[i,j,k,l] ~ dbern(mu.p[i,j,k,l])
    
    } #nseason
    } #nrep
    } #nsite
    } #nspec
}
      ",fill=TRUE)

sink()
zst <- array(1,c(nsite, nspec,2))

# valores iniciales
#inits <- function()     list(z=zst  )

inits <- function() list(z=zst,lpsi=rnorm(roed.data.pa$nspec), b1=rnorm(roed.data.pa$nspec),
                         b2=rnorm(roed.data.pa$nspec),#lp=c(rnorm(roed.data.pa$nrep),rnorm(roed.data.pa$nspec)),
                         mu.lpsi=rnorm(1),mu.b1=rnorm(1),mu.b2=rnorm(1),
                         mu.lp=rnorm(1))
# ,sd.lpsi=runif(roed.data$nspec,0.1,5),
# sd.lp=runif(roed.data$nspec,0.1,5), sd.b1=runif(roed.data$nspec,0.1,5))


params1 <- c("lpsi", "b1", "lp","mu.lpsi","mu.lp","mu.b1")
params2 <- c("Nsite", "Nocc.fs")

# ajustes de MCMC
ni <- 50000
nt <- 10
nb <- 25000
nc <- 3

# llamar JAGS desde R 
out150.a= jags(roed.data.pa,inits, params1, "mod150_pa_r1.jags", n.chains=nc, 
                  n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out150.b= jags(roed.data.pa,inits, params2, "mod150_pa_r1.jags", n.chains=nc, 
                  n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out150.a, file='out150_pa1_r1.rda') 
save(out150.b, file='out150_pa2_r1.rda') 


###################################

### JAGS code 200 m 
sink("mod200_pa_r1.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~   dnorm(0, 0.1)
    mu.lp ~   dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    
    # priors for species-specific effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    
    for (j in 1:nrep) {  
    lp[j,k] ~ dnorm (mu.lp, tau.lp)    
    
    } #nrep
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*pa200[i]
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    for (l in 1:nseason){                                         #loop sobre tiempo
    for (j in 1:nrep) {  
    logit(p[j,k,l]) <- lp[j,k]     #poner por fuera de i
    } #season
    }
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (j in 1:nrep) {                                           #loop sobre noche 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    mu.p[i,j,k,l] <- p[j,k,l]*z[i,k,l]
    y[i,j,k,l] ~ dbern(mu.p[i,j,k,l])
    
    } #nseason
    } #nrep
    } #nsite
    } #nspec
    }
    ",fill=TRUE)

sink()

# llamar JAGS desde R 
out200.a= jags(roed.data.pa,inits, params1, "mod200_pa_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out200.b= jags(roed.data.pa,inits, params2, "mod200_pa_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out200.a, file='out200_pa1_r1.rda') 
save(out200.b, file='out200_pa2_r1.rda') 


####################
### JAGS code 300 m 
sink("mod300_pa_r1.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~   dnorm(0, 0.1)
    mu.lp ~   dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    
    # priors for species-specific effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    
    for (j in 1:nrep) {  
    lp[j,k] ~ dnorm (mu.lp, tau.lp)    
    
    } #nrep
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*pa300[i]
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    for (l in 1:nseason){                                         #loop sobre tiempo
    for (j in 1:nrep) {  
    logit(p[j,k,l]) <- lp[j,k]     #poner por fuera de i
    } #season
    }
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (j in 1:nrep) {                                           #loop sobre noche 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    mu.p[i,j,k,l] <- p[j,k,l]*z[i,k,l]
    y[i,j,k,l] ~ dbern(mu.p[i,j,k,l])
    
    } #nseason
    } #nrep
    } #nsite
    } #nspec
    }
    ",fill=TRUE)

sink()

# llamar JAGS desde R 
out300.a= jags(roed.data.pa,inits, params1, "mod300_pa_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out300.b= jags(roed.data.pa,inits, params2, "mod300_pa_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out300.a, file='out300_pa1_r1.rda') 
save(out300.b, file='out300_pa2_r1.rda') 


####################
### JAGS code 400 m 
sink("mod400_pa_r1.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~   dnorm(0, 0.1)
    mu.lp ~   dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    
    # priors for species-specific effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    
    for (j in 1:nrep) {  
    lp[j,k] ~ dnorm (mu.lp, tau.lp)    
    
    } #nrep
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*pa400[i]
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    for (l in 1:nseason){                                         #loop sobre tiempo
    for (j in 1:nrep) {  
    logit(p[j,k,l]) <- lp[j,k]     #poner por fuera de i
    } #season
    }
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (j in 1:nrep) {                                           #loop sobre noche 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    mu.p[i,j,k,l] <- p[j,k,l]*z[i,k,l]
    y[i,j,k,l] ~ dbern(mu.p[i,j,k,l])
    
    } #nseason
    } #nrep
    } #nsite
    } #nspec
    }
    ",fill=TRUE)

sink()

# llamar JAGS desde R 
out400.a= jags(roed.data.pa,inits, params1, "mod400_pa_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out400.b= jags(roed.data.pa,inits, params2, "mod400_pa_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out400.a, file='out400_pa1_r1.rda') 
save(out400.b, file='out400_pa2_r1.rda') 
##################

####################
### JAGS code 500 m 
sink("mod500_pa_r1.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~   dnorm(0, 0.1)
    mu.lp ~   dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    
    # priors for species-specific effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    
    for (j in 1:nrep) {  
    lp[j,k] ~ dnorm (mu.lp, tau.lp)    
    
    } #nrep
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*pa500[i]
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    for (l in 1:nseason){                                         #loop sobre tiempo
    for (j in 1:nrep) {  
    logit(p[j,k,l]) <- lp[j,k]     #poner por fuera de i
    } #season
    }
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (j in 1:nrep) {                                           #loop sobre noche 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    mu.p[i,j,k,l] <- p[j,k,l]*z[i,k,l]
    y[i,j,k,l] ~ dbern(mu.p[i,j,k,l])
    
    } #nseason
    } #nrep
    } #nsite
    } #nspec
    }
    ",fill=TRUE)

sink()

# llamar JAGS desde R 
out500.a= jags(roed.data.pa,inits, params1, "mod500_pa_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out500.b= jags(roed.data.pa,inits, params2, "mod500_pa_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out500.a, file='out500_pa1_r1.rda') 
save(out500.b, file='out500_pa2_r1.rda') 
##################


####################
### JAGS code 600 m 
sink("mod600_pa_r1.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~   dnorm(0, 0.1)
    mu.lp ~   dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    
    # priors for species-specific effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    
    for (j in 1:nrep) {  
    lp[j,k] ~ dnorm (mu.lp, tau.lp)    
    
    } #nrep
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*pa600[i]
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    for (l in 1:nseason){                                         #loop sobre tiempo
    for (j in 1:nrep) {  
    logit(p[j,k,l]) <- lp[j,k]     #poner por fuera de i
    } #season
    }
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (j in 1:nrep) {                                           #loop sobre noche 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    mu.p[i,j,k,l] <- p[j,k,l]*z[i,k,l]
    y[i,j,k,l] ~ dbern(mu.p[i,j,k,l])
    
    } #nseason
    } #nrep
    } #nsite
    } #nspec
    }
    ",fill=TRUE)

sink()

# llamar JAGS desde R 
out600.a= jags(roed.data.pa,inits, params1, "mod600_pa_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out600.b= jags(roed.data.pa,inits, params2, "mod600_pa_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out600.a, file='out600_pa1_r1.rda') 
save(out600.b, file='out600_pa2_r1.rda') 
##################


########################################## 
### IS
##############################

str(roed.data.is <- list(y=y,J=J, nsite=dim(y)[1],nspec=dim(y)[3],nrep=dim(y)[2],
                         nseason=dim(y)[4], is150=is150, is200=is200, is300=is300,is400=is400, 
                         is500=is500, is600=is600))
###Modelos solo PA
### JAGS code
sink("mod150_is_r1.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~   dnorm(0, 0.1)
    mu.lp ~   dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    
    # priors for species-specific effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    
    for (j in 1:nrep) {  
    lp[j,k] ~ dnorm (mu.lp, tau.lp)    
    
    } #nrep
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*is150[i]
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    for (l in 1:nseason){                                         #loop sobre tiempo
    for (j in 1:nrep) {  
    logit(p[j,k,l]) <- lp[j,k]     #poner por fuera de i
    } #season
    }
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (j in 1:nrep) {                                           #loop sobre noche 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    mu.p[i,j,k,l] <- p[j,k,l]*z[i,k,l]
    y[i,j,k,l] ~ dbern(mu.p[i,j,k,l])
    
    } #nseason
    } #nrep
    } #nsite
    } #nspec
    }
    ",fill=TRUE)

sink()
zst <- array(1,c(nsite, nspec,2))

# valores iniciales
#inits <- function()     list(z=zst  )

inits <- function() list(z=zst,lpsi=rnorm(roed.data.pa$nspec), b1=rnorm(roed.data.pa$nspec),
                         b2=rnorm(roed.data.pa$nspec),#lp=c(rnorm(roed.data.pa$nrep),rnorm(roed.data.pa$nspec)),
                         mu.lpsi=rnorm(1),mu.b1=rnorm(1),mu.b2=rnorm(1),
                         mu.lp=rnorm(1))
# ,sd.lpsi=runif(roed.data$nspec,0.1,5),
# sd.lp=runif(roed.data$nspec,0.1,5), sd.b1=runif(roed.data$nspec,0.1,5))


params1 <- c("lpsi", "b1", "lp","mu.lpsi","mu.lp","mu.b1")
params2 <- c("Nsite", "Nocc.fs")

# ajustes de MCMC
ni <- 100000
nt <- 2
nb <- 50000
nc <- 3              #ojo!!! de pa me corrio 25 t 25

# llamar JAGS desde R 
out150.is.a= jags(roed.data.is,inits, params1, "mod150_is_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out150.is.b= jags(roed.data.is,inits, params2, "mod150_is_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out150.is.a, file='out150_is1_r1.rda') 
save(out150.is.b, file='out150_is2_r1.rda') 

###################################
### JAGS code 200 m 
sink("mod200_is_r1.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~   dnorm(0, 0.1)
    mu.lp ~   dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    
    # priors for species-specific effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    
    for (j in 1:nrep) {  
    lp[j,k] ~ dnorm (mu.lp, tau.lp)    
    
    } #nrep
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*is200[i]
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    for (l in 1:nseason){                                         #loop sobre tiempo
    for (j in 1:nrep) {  
    logit(p[j,k,l]) <- lp[j,k]     #poner por fuera de i
    } #season
    }
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (j in 1:nrep) {                                           #loop sobre noche 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    mu.p[i,j,k,l] <- p[j,k,l]*z[i,k,l]
    y[i,j,k,l] ~ dbern(mu.p[i,j,k,l])
    
    } #nseason
    } #nrep
    } #nsite
    } #nspec
    }
    ",fill=TRUE)

sink()

# llamar JAGS desde R 
out200.is.a= jags(roed.data.is,inits, params1, "mod200_is_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out200.is.b= jags(roed.data.is,inits, params2, "mod200_is_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out200.is.a, file='out200_is1_r1.rda') 
save(out200.is.b, file='out200_is2_r1.rda') 

####################
### JAGS code 300 m 
sink("mod300_is_r1.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~   dnorm(0, 0.1)
    mu.lp ~   dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    
    # priors for species-specific effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    
    for (j in 1:nrep) {  
    lp[j,k] ~ dnorm (mu.lp, tau.lp)    
    
    } #nrep
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*is300[i]
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    for (l in 1:nseason){                                         #loop sobre tiempo
    for (j in 1:nrep) {  
    logit(p[j,k,l]) <- lp[j,k]     #poner por fuera de i
    } #season
    }
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (j in 1:nrep) {                                           #loop sobre noche 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    mu.p[i,j,k,l] <- p[j,k,l]*z[i,k,l]
    y[i,j,k,l] ~ dbern(mu.p[i,j,k,l])
    
    } #nseason
    } #nrep
    } #nsite
    } #nspec
    }
    ",fill=TRUE)

sink()

# llamar JAGS desde R 
out300.is.a= jags(roed.data.is,inits, params1, "mod300_is_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out300.is.b= jags(roed.data.is,inits, params2, "mod300_is_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out300.is.a, file='out300_is1_r1.rda') 
save(out300.is.b, file='out300_is2_r1.rda') 


####################
### JAGS code 400 m 
sink("mod400_is_r1.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~   dnorm(0, 0.1)
    mu.lp ~   dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    
    # priors for species-specific effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    
    for (j in 1:nrep) {  
    lp[j,k] ~ dnorm (mu.lp, tau.lp)    
    
    } #nrep
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*is400[i]
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    for (l in 1:nseason){                                         #loop sobre tiempo
    for (j in 1:nrep) {  
    logit(p[j,k,l]) <- lp[j,k]     #poner por fuera de i
    } #season
    }
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (j in 1:nrep) {                                           #loop sobre noche 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    mu.p[i,j,k,l] <- p[j,k,l]*z[i,k,l]
    y[i,j,k,l] ~ dbern(mu.p[i,j,k,l])
    
    } #nseason
    } #nrep
    } #nsite
    } #nspec
    }
    ",fill=TRUE)

sink()

# llamar JAGS desde R 
out400.is.a= jags(roed.data.is,inits, params1, "mod400_is_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out400.is.b= jags(roed.data.is,inits, params2, "mod400_is_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out400.is.a, file='out400_is1_r1.rda') 
save(out400.is.b, file='out400_is2_r1.rda') 
##################

####################
### JAGS code 500 m 
sink("mod500_is_r1.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~   dnorm(0, 0.1)
    mu.lp ~   dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    
    # priors for species-specific effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    
    for (j in 1:nrep) {  
    lp[j,k] ~ dnorm (mu.lp, tau.lp)    
    
    } #nrep
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*is500[i]
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    for (l in 1:nseason){                                         #loop sobre tiempo
    for (j in 1:nrep) {  
    logit(p[j,k,l]) <- lp[j,k]     #poner por fuera de i
    } #season
    }
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (j in 1:nrep) {                                           #loop sobre noche 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    mu.p[i,j,k,l] <- p[j,k,l]*z[i,k,l]
    y[i,j,k,l] ~ dbern(mu.p[i,j,k,l])
    
    } #nseason
    } #nrep
    } #nsite
    } #nspec
    }
    ",fill=TRUE)

sink()

# llamar JAGS desde R 
out500.is.a= jags(roed.data.is,inits, params1, "mod500_is_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out500.is.b= jags(roed.data.is,inits, params2, "mod500_is_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out500.is.a, file='out500_is1_r1.rda') 
save(out500.is.b, file='out500_is2_r1.rda') 

##################
## JAGS code 600 m 
sink("mod600_is_r1.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~   dnorm(0, 0.1)
    mu.lp ~   dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    
    # priors for species-specific effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    
    for (j in 1:nrep) {  
    lp[j,k] ~ dnorm (mu.lp, tau.lp)    
    
    } #nrep
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*is600[i]
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    for (l in 1:nseason){                                         #loop sobre tiempo
    for (j in 1:nrep) {  
    logit(p[j,k,l]) <- lp[j,k]     #poner por fuera de i
    } #season
    }
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (j in 1:nrep) {                                           #loop sobre noche 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    mu.p[i,j,k,l] <- p[j,k,l]*z[i,k,l]
    y[i,j,k,l] ~ dbern(mu.p[i,j,k,l])
    
    } #nseason
    } #nrep
    } #nsite
    } #nspec
    }
    ",fill=TRUE)

sink()

# llamar JAGS desde R 
out600.is.a= jags(roed.data.is,inits, params1, "mod600_is_r1.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out600.is.b= jags(roed.data.is,inits, params2, "mod600_is_r1.jags", n.chains=nc, 
              n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out600.is.a, file='out600_is1_r1.rda') 
save(out600.is.b, file='out600_is2_r1.rda') 
##################















































###############Modelos con PA e IS todo junto
### JAGS code
sink("mod150-ispa.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~   dnorm(0, 0.1)
    mu.b2 ~   dnorm(0, 0.1)
    mu.lp ~   dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    tau.b2 <- pow(sd.b2,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    sd.b2 ~ dunif(0,1)
    
    # priors for species-specific effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    lp[k] ~ dnorm (mu.lp, tau.lp)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    b2[k] ~ dnorm (mu.b2, tau.b2)
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*is150[i] + b2[k]*pa150[i]
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    for (l in 1:nseason){                                         #loop sobre tiempo
    logit(p[k,l]) <- lp[k]    #poner por fuera de i
    } #season
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    mu.p[i,k,l] <- p[k,l]*z[i,k,l]
    ysum[i,k,l] ~ dbin(mu.p[i,k,l],J[i,l])
    
    } #nseason
    } #nsite
    } #nspec
    
    #Derived quantities
    for (i in 1:nsite) {
    for(l in 1:nseason){
    Nsite [i,l] <- sum (z[i,,l])
    }
    }
    for (k in 1:nspec){
    for(l in 1:nseason){
    Nocc.fs [k,l] <- sum (z[,k,l])
    }   
    }    
    } #model
    
    ",fill=TRUE)

sink()
zst <- apply(y,c(1,3,4),max)
zst[is.na(zst)]<-1
# valores iniciales
#inits <- function()     list(z=zst  )

inits <- function() list(z=zst,lpsi=rnorm(roed.data$nspec), b1=rnorm(roed.data$nspec),
                         b2=rnorm(roed.data$nspec),lp=rnorm(roed.data$nspec),mu.lpsi=rnorm(1),mu.b1=rnorm(1),mu.b2=rnorm(1),
                         mu.lp=rnorm(1))
# ,sd.lpsi=runif(roed.data$nspec,0.1,5),
# sd.lp=runif(roed.data$nspec,0.1,5), sd.b1=runif(roed.data$nspec,0.1,5))


params1 <- c("lpsi", "b1","b2", "lp","mu.lpsi","mu.lp","mu.b1")
params2 <- c("Nsite", "Nocc.fs")

# ajustes de MCMC
ni <- 100000
nt <- 2
nb <- 50000
nc <- 3

# llamar JAGS desde R 
out150.v1a2= jags(roed.data,inits, params1, "mod150-ispa.jags", n.chains=nc, 
                  n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out150.v1b2= jags(roed.data,inits, params2, "mod150-ispa.jags", n.chains=nc, 
                  n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out150.v1a2, file='out150-ispa1.rda') 
save(out150.v1b2, file='out150.ispa2.rda') 


###################################
##############
# mod 2 - 200m
sink("mod200-ispa.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1   ~ dnorm(0, 0.1)
    mu.b2   ~ dnorm(0, 0.1)
    mu.lp   ~ dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp   <- pow(sd.lp,-2)   
    tau.b1   <- pow(sd.b1,-2)
    tau.b2   <- pow(sd.b2,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp   ~ dunif(0,1)
    sd.b1   ~ dunif(0,1)
    sd.b2   ~ dunif(0,1)
    
    # priors for species-specific time effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    lp[k] ~ dnorm (mu.lp, tau.lp)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    b2[k] ~ dnorm (mu.b2, tau.b2)
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*is200[i] + b2[k]*pa200[i]
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    
    for (l in 1:nseason){                                         #loop sobre tiempo
    logit(p[k,l]) <- lp[k]    #poner por fuera de i
    } #season
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    
    mu.p[i,k,l] <- p[k,l]*z[i,k,l]
    ysum[i,k,l] ~ dbin(mu.p[i,k,l],J[i,l])
    
    } #season
    } #nsite
    } #nspec
    
    # Derived quantities
    for (i in 1:nsite) {
    for(l in 1:nseason){
    Nsite [i,l] <- sum (z[i,,l])
    }
    }
    for (k in 1:nspec){
    for(l in 1:nseason){
    Nocc.fs [k,l] <- sum (z[,k,l])
    }
    }
    } #model
    
    
    ",fill=TRUE)


sink()
zst <- apply(y,c(1,3,4),max)
zst[is.na(zst)]<-1
# valores iniciales
#inits <- function()     list(z =zst)
inits <- function() list(z=zst,lpsi=rnorm(roed.data$nspec), b1=rnorm(roed.data$nspec),
                         b2=rnorm(roed.data$nspec),lp=rnorm(roed.data$nspec),mu.lpsi=rnorm(1),mu.b1=rnorm(1),
                         mu.lp=rnorm(1),mu.b2=rnorm(1))
params1 <- c("lpsi", "b1", "b2","lp","mu.lpsi","mu.b1","mu.lp")
params2 <- c("Nsite", "Nocc.fs")

#ajustes de MCMC
ni <- 100000
nt <- 2
nb <- 50000
nc <- 3

# llamar JAGS desde R 
out200a2= jags(roed.data,inits, params1, "mod200-ispa.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out200b2= jags(roed.data,inits, params2, "mod200-ispa.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out200a2, file='out200-ispa1.rda') 
save(out200b2, file='out200-ispa2.rda') 

##############
# mod 3 - 300m
sink("mod300-ispa.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~ dnorm(0, 0.1)
    mu.b2 ~ dnorm(0, 0.1)
    mu.lp ~ dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    tau.b2 <- pow(sd.b2,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    sd.b2 ~ dunif(0,1)
    
    # priors for species-specific time effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    lp[k] ~ dnorm (mu.lp, tau.lp)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    b2[k] ~ dnorm (mu.b2, tau.b2)
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*is300[i] + b2[k]*pa300[i]
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    
    for (l in 1:nseason){                                         #loop sobre tiempo
    logit(p[k,l]) <- lp[k]    #poner por fuera de i
    } #season
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    
    mu.p[i,k,l] <- p[k,l]*z[i,k,l]
    ysum[i,k,l] ~ dbin(mu.p[i,k,l],J[i,l])
    
    } #season
    } #nsite
    } #nspec
    
    # Derived quantities
    for (i in 1:nsite) {
    for(l in 1:nseason){
    Nsite [i,l] <- sum (z[i,,l])
    }
    }
    for (k in 1:nspec){
    for(l in 1:nseason){
    Nocc.fs [k,l] <- sum (z[,k,l])
    }
    }
    } #model
    
    ",fill=TRUE)

sink()
zst <- apply(y,c(1,3,4),max)
zst[is.na(zst)]<-1
# valores iniciales
#inits <- function()     list(z =zst)
inits <- function() list(z=zst,lpsi=rnorm(roed.data$nspec), b1=rnorm(roed.data$nspec),
                         b2=rnorm(roed.data$nspec),lp=rnorm(roed.data$nspec),mu.lpsi=rnorm(1),mu.b1=rnorm(1),mu.b2=rnorm(1),
                         mu.lp=rnorm(1))

params1 <- c("lpsi", "b1","b2","lp","mu.lpsi","mu.b1","mu.lp")
params2 <- c("Nsite", "Nocc.fs")

#ajustes de MCMC
ni <- 100000
nt <- 2
nb <- 50000
nc <- 3

# llamar JAGS desde R 
out300a2= jags(roed.data,inits, params1, "mod300-ispa.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out300b2= jags(roed.data,inits, params2, "mod300-ispa.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out300a2, file='out300-ispa1.rda') 
save(out300b2, file='out300-ispa2.rda') 

####### mod 4 - 400m
sink("mod400-ispa.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~ dnorm(0, 0.1)
    mu.b2 ~ dnorm(0, 0.1)
    mu.lp ~ dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    tau.b2 <- pow(sd.b2,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    sd.b2 ~ dunif(0,1)
    
    # priors for species-specific time effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    lp[k] ~ dnorm (mu.lp, tau.lp)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    b2[k] ~ dnorm (mu.b2, tau.b2)
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*is400[i] + b2[k]*pa400[i] 
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    
    for (l in 1:nseason){                                         #loop sobre tiempo
    logit(p[k,l]) <- lp[k]    #poner por fuera de i
    } #season
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    
    mu.p[i,k,l] <- p[k,l]*z[i,k,l]
    ysum[i,k,l] ~ dbin(mu.p[i,k,l],J[i,l])
    
    } #season
    } #nsite
    } #nspec
    
    # Derived quantities
    for (i in 1:nsite) {
    for(l in 1:nseason){
    Nsite [i,l] <- sum (z[i,,l])
    }
    }   
    for (k in 1:nspec){
    for(l in 1:nseason){
    Nocc.fs [k,l] <- sum (z[,k,l])
    }
    }
    } #model
    
    ",fill=TRUE)


sink()
zst <- apply(y,c(1,3,4),max)
zst[is.na(zst)]<-1
# valores iniciales
# inits <- function()     list(z =zst)
inits <- function() list(z=zst,lpsi=rnorm(roed.data$nspec), b1=rnorm(roed.data$nspec),
                         b2=rnorm(roed.data$nspec),lp=rnorm(roed.data$nspec),mu.lpsi=rnorm(1),mu.b1=rnorm(1),mu.b2=rnorm(1),
                         mu.lp=rnorm(1))
params1 <- c("lpsi", "b1", "b2","lp","mu.lpsi","mu.b1","mu.lp")
params2 <- c("Nsite", "Nocc.fs")

#ajustes de MCMC
ni <- 100000
nt <- 2
nb <- 50000
nc <- 3

# llamar JAGS desde R 
out400a2= jags(roed.data,inits, params1, "mod400-ispa.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out400b2= jags(roed.data,inits, params2, "mod400-ispa.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out400a2, file='out400-ispa1.rda') 
save(out400b2, file='out400-ispa2.rda') 

########################################
##############
# mod 5 - 500m
sink("mod500-ispa.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~ dnorm(0, 0.1)
    mu.b2 ~ dnorm(0, 0.1)
    mu.lp ~ dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    tau.b2 <- pow(sd.b2,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    sd.b2 ~ dunif(0,1)
    
    # priors for species-specific time effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    lp[k] ~ dnorm (mu.lp, tau.lp)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    b2[k] ~ dnorm (mu.b2, tau.b2)
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*is500[i] + b2[k]*pa500[i] 
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    
    for (l in 1:nseason){                                         #loop sobre tiempo
    logit(p[k,l]) <- lp[k]    #poner por fuera de i
    } #season
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    
    mu.p[i,k,l] <- p[k,l]*z[i,k,l]
    ysum[i,k,l] ~ dbin(mu.p[i,k,l],J[i,l])
    
    } #season
    } #nsite
    } #nspec
    
    # Derived quantities
    for (i in 1:nsite) {
    for(l in 1:nseason){
    Nsite [i,l] <- sum (z[i,,l])
    }
    }
    for (k in 1:nspec){
    for(l in 1:nseason){
    Nocc.fs [k,l] <- sum (z[,k,l])
    }
    }
    } #model
    
    ",fill=TRUE)

sink()
zst <- apply(y,c(1,3,4),max)
zst[is.na(zst)]<-1
# valores iniciales
#inits <- function()     list(z =zst)
inits <- function() list(z=zst,lpsi=rnorm(roed.data$nspec), b1=rnorm(roed.data$nspec),
                         b2=rnorm(roed.data$nspec),lp=rnorm(roed.data$nspec),mu.lpsi=rnorm(1),mu.b1=rnorm(1),
                         mu.b2=rnorm(1),             mu.lp=rnorm(1))

params1 <- c("lpsi", "b1", "b2","lp","mu.lpsi","mu.b1","mu.lp")
params2 <- c("Nsite", "Nocc.fs")

#ajustes de MCMC
ni <- 100000
nt <- 2
nb <- 50000
nc <- 3

# llamar JAGS desde R 
out500a2= jags(roed.data,inits, params1, "mod500-ispa.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out500b2= jags(roed.data,inits, params2, "mod500-ispa.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out500a2, file='out500-ispa1.rda') 
save(out500b2, file='out500-ispa2.rda') 

##################################################
# mod 6 - 600m
sink("mod600-ispa.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~ dnorm(0, 0.1)
    mu.b2 ~ dnorm(0, 0.1)
    mu.lp ~ dnorm(0, 0.1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)
    tau.b2 <- pow(sd.b2,-2)
    
    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)
    sd.b2 ~ dunif(0,1)
    
    # priors for species-specific time effects
    for (k in 1:nspec){
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    lp[k] ~ dnorm (mu.lp, tau.lp)
    b1[k] ~ dnorm (mu.b1, tau.b1)
    b2[k] ~ dnorm (mu.b2, tau.b2)
    
    # Ecological model, process model (true occurrence at site i) 
    ## time as random effect
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    logit(psi[i,k,l]) <- lpsi[k] + b1[k]*is600[i] + b2[k]*pa600[i] 
    
    z[i,k,l] ~ dbern(psi[i,k,l])
    } #season
    } #site
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    
    for (l in 1:nseason){                                         #loop sobre tiempo
    logit(p[k,l]) <- lp[k]    #poner por fuera de i
    } #season
    
    for (i in 1:nsite) {                                           #loop sobre sitios 
    for (l in 1:nseason){                                         #loop sobre tiempo 
    
    mu.p[i,k,l] <- p[k,l]*z[i,k,l]
    ysum[i,k,l] ~ dbin(mu.p[i,k,l],J[i,l])
    
    } #season
    } #nsite
    } #nspec
    
    # Derived quantities
    for (i in 1:nsite) {
    for(l in 1:nseason){
    Nsite [i,l] <- sum (z[i,,l])
    }
    }
    for (k in 1:nspec){
    for(l in 1:nseason){
    Nocc.fs [k,l] <- sum (z[,k,l])
    }
    }
    } #model
    
    ",fill=TRUE)

sink()
zst <- apply(y,c(1,3,4),max)
zst[is.na(zst)]<-1
# valores iniciales
#inits <- function()     list(z =zst)
inits <- function() list(z=zst,lpsi=rnorm(roed.data$nspec), b1=rnorm(roed.data$nspec),
                         b2=rnorm(roed.data$nspec),
                         lp=rnorm(roed.data$nspec),mu.lpsi=rnorm(1),mu.b1=rnorm(1),
                         mu.b2=rnorm(1),mu.lp=rnorm(1))

params1 <- c("lpsi", "b1","b2", "lp","mu.lpsi","mu.b1","mu.lp")
params2 <- c("Nsite", "Nocc.fs")

#ajustes de MCMC
ni <- 100000
nt <- 2
nb <- 50000
nc <- 3

# llamar JAGS desde R 
out600a2= jags(roed.data,inits, params1, "mod600-ispa.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
out600b2= jags(roed.data,inits, params2, "mod600-ispa.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

save(out600a2, file='out600-ispa1.rda') 
save(out600b2, file='out600-ispa2.rda') 
# 
#######################################

load(file='out150_is1_r1.rda')  #out150.is.a
load(file='out150_is2_r1.rda')  #out150.is.b
load(file='out200_is1_r1.rda')  #out200.is.a
load(file='out200_is2_r1.rda')  #out200.is.b
load(file='out300_is1_r1.rda')  #out300.is.a
load(file='out300_is2_r1.rda')  #out300.is.b
load(file='out400_is1_r1.rda')  #out400.is.a
load(file='out400_is2_r1.rda')  #out400.is.b
load(file='out500_is1_r1.rda')  #out500.is.a
load(file='out500_is2_r1.rda')  #out500.is.b
load(file='out600_is1_r1.rda')  #out600.is.a
load(file='out600_is2_r1.rda')  #out600.is.b




###############################
# Miodel Regression para Riqueza e IS
is150 <- (scale(IS150[1:50]))
is200 <- (scale(IS200[1:50]))
is300 <- (scale(IS300[1:50]))
is400 <- (scale(IS400[1:50]))
is500 <- (scale(IS500[1:50]))
is600 <- (scale(IS600[1:50]))


#promedio de riqueza y intervalos para todas las escalas
N.pm.150 <-(apply(out150.is.b$mean$Nsite,1,mean))
N.psd.150 <- apply(out150.is.b$sd$Nsite,1,mean)
N.cri.150<- apply(out150.is.b$sims.list$Nsite, 2,function(x) quantile(x,prob=c(0.125,0.875)))

N.pm.200 <-(apply(out200.is.b$mean$Nsite,1,mean))
N.psd.200 <- apply(out200.is.b$sd$Nsite,1,mean)
N.cri.200<- apply(out200.is.b$sims.list$Nsite, 2,function(x) quantile(x,prob=c(0.125,0.875)))

N.pm.300 <-(apply(out300.is.b$mean$Nsite,1,mean))
N.psd.300 <- apply(out300.is.b$sd$Nsite,1,mean)
N.cri.300<- apply(out300.is.b$sims.list$Nsite, 2,function(x) quantile(x,prob=c(0.125,0.875)))

N.pm.400 <-(apply(out400.is.b$mean$Nsite,1,mean))
N.psd.400 <- apply(out400.is.b$sd$Nsite,1,mean)
N.cri.400<- apply(out400.is.b$sims.list$Nsite, 2,function(x) quantile(x,prob=c(0.125,0.875)))

N.pm.500 <-(apply(out500.is.b$mean$Nsite,1,mean))
N.psd.500 <- apply(out500.is.b$sd$Nsite,1,mean)
N.cri.500<- apply(out500.is.b$sims.list$Nsite, 2,function(x) quantile(x,prob=c(0.125,0.875)))

N.pm.600 <-(apply(out600.is.b$mean$Nsite,1,mean))
N.psd.600 <- apply(out600.is.b$sd$Nsite,1,mean)
N.cri.600<- apply(out600.is.b$sims.list$Nsite, 2,function(x) quantile(x,prob=c(0.125,0.875)))


is.pred<- seq(-1.436,2.6,length=500)

str(win.data<- list(is=is150,N=N.pm.150,psd=N.psd.150, n= length(N.pm.150), pred.is=is.pred, npred=length(is.pred)))

sink("meta.analysis.150_revision.txt")
cat("
model{
  
  #Priors
  for (v in 1:3) {     #previas for intercep and plonomial coefficients
    beta[v] ~dnorm (0,0.0001)
  }
  tau.site<-pow(sd.site, -2)
  sd.site~dunif(0,10)
  
  # Likelihood
  for (i in 1:n) {
    N[i]~dnorm(muN[i], tau.psd[i])
    tau.psd[i]<-pow(psd[i],-2)
    muN[i]<- beta[1]+ beta[2]*is[i] +beta[3]*pow(is[i],2) + eps.site[i]
    eps.site [i] ~ dnorm (0, tau.site)
  }
  
  #get predictions for plot
  for (i in 1:npred){
    Npred[i] <- beta[1] + beta[2]*pred.is[i] + beta[3]*pow(pred.is[i],2) 
    
  }
}
", fill= TRUE)
sink()

inits<- function() list(beta=rnorm(3))
params<- c("beta", "sd.site","Npred")
ni<- 12000 ; nt<- 10 ; nb <- 2000 ; nc<- 3

out.150<- jags(win.data, params,inits = inits, "meta.analysis.150_revision.txt", n.chains = nc,
           n.thin = nt, n.iter = ni, n.burnin = nb)


str(win.data<- list(is=is200,N=N.pm.200,psd=N.psd.200, n= length(N.pm.200), pred.is=is.pred, npred=length(is.pred)))

sink("meta.analysis.200_revision.txt")
cat("
    model{
    
    #Priors
    for (v in 1:3) {     #previas for intercep and plonomial coefficients
    beta[v] ~dnorm (0,0.0001)
    }
    tau.site<-pow(sd.site, -2)
    sd.site~dunif(0,10)
    
    # Likelihood
    for (i in 1:n) {
    N[i]~dnorm(muN[i], tau.psd[i])
    tau.psd[i]<-pow(psd[i],-2)
    muN[i]<- beta[1]+ beta[2]*is[i] +beta[3]*pow(is[i],2) + eps.site[i]
    eps.site [i] ~ dnorm (0, tau.site)
    }
    
    #get predictions for plot
    for (i in 1:npred){
    Npred[i] <- beta[1] + beta[2]*pred.is[i] + beta[3]*pow(pred.is[i],2) 
    
    }
    }
    ", fill= TRUE)
sink()

inits<- function() list(beta=rnorm(3))
params<- c("beta", "sd.site","Npred")
ni<- 12000 ; nt<- 10 ; nb <- 2000 ; nc<- 3

out.200<- jags(win.data, params,inits = inits, "meta.analysis.200_revision.txt", n.chains = nc,
           n.thin = nt, n.iter = ni, n.burnin = nb)


str(win.data<- list(is=is300,N=N.pm.300,psd=N.psd.300, n= length(N.pm.300), pred.is=is.pred, npred=length(is.pred)))

sink("meta.analysis.300_revision.txt")
cat("
    model{
    
    #Priors
    for (v in 1:3) {     #previas for intercep and plonomial coefficients
    beta[v] ~dnorm (0,0.0001)
    }
    tau.site<-pow(sd.site, -2)
    sd.site~dunif(0,10)
    
    # Likelihood
    for (i in 1:n) {
    N[i]~dnorm(muN[i], tau.psd[i])
    tau.psd[i]<-pow(psd[i],-2)
    muN[i]<- beta[1]+ beta[2]*is[i] +beta[3]*pow(is[i],2) + eps.site[i]
    eps.site [i] ~ dnorm (0, tau.site)
    }
    
    #get predictions for plot
    for (i in 1:npred){
    Npred[i] <- beta[1] + beta[2]*pred.is[i] + beta[3]*pow(pred.is[i],2) 
    
    }
    }
    ", fill= TRUE)
sink()

inits<- function() list(beta=rnorm(3))
params<- c("beta", "sd.site","Npred")
ni<- 12000 ; nt<- 10 ; nb <- 2000 ; nc<- 3

out.300<- jags(win.data, params,inits = inits, "meta.analysis.300_revision.txt", n.chains = nc,
           n.thin = nt, n.iter = ni, n.burnin = nb)


str(win.data<- list(is=is400,N=N.pm.400,psd=N.psd.400, n= length(N.pm.400), pred.is=is.pred, npred=length(is.pred)))

sink("meta.analysis.400_revision.txt")
cat("
    model{
    
    #Priors
    for (v in 1:3) {     #previas for intercep and plonomial coefficients
    beta[v] ~dnorm (0,0.0001)
    }
    tau.site<-pow(sd.site, -2)
    sd.site~dunif(0,10)
    
    # Likelihood
    for (i in 1:n) {
    N[i]~dnorm(muN[i], tau.psd[i])
    tau.psd[i]<-pow(psd[i],-2)
    muN[i]<- beta[1]+ beta[2]*is[i] +beta[3]*pow(is[i],2) + eps.site[i]
    eps.site [i] ~ dnorm (0, tau.site)
    }
    
    #get predictions for plot
    for (i in 1:npred){
    Npred[i] <- beta[1] + beta[2]*pred.is[i] + beta[3]*pow(pred.is[i],2) 
    
    }
    }
    ", fill= TRUE)
sink()

inits<- function() list(beta=rnorm(3))
params<- c("beta", "sd.site","Npred")
ni<- 12000 ; nt<- 10 ; nb <- 2000 ; nc<- 3

out.400<- jags(win.data, params,inits = inits, "meta.analysis.400_revision.txt", n.chains = nc,
           n.thin = nt, n.iter = ni, n.burnin = nb)



str(win.data<- list(is=is500,N=N.pm.500,psd=N.psd.500, n= length(N.pm.500), pred.is=is.pred, npred=length(is.pred)))

sink("meta.analysis.500_revision.txt")
cat("
    model{
    
    #Priors
    for (v in 1:3) {     #previas for intercep and plonomial coefficients
    beta[v] ~dnorm (0,0.0001)
    }
    tau.site<-pow(sd.site, -2)
    sd.site~dunif(0,10)
    
    # Likelihood
    for (i in 1:n) {
    N[i]~dnorm(muN[i], tau.psd[i])
    tau.psd[i]<-pow(psd[i],-2)
    muN[i]<- beta[1]+ beta[2]*is[i] +beta[3]*pow(is[i],2) + eps.site[i]
    eps.site [i] ~ dnorm (0, tau.site)
    }
    
    #get predictions for plot
    for (i in 1:npred){
    Npred[i] <- beta[1] + beta[2]*pred.is[i] + beta[3]*pow(pred.is[i],2) 
    
    }
    }
    ", fill= TRUE)
sink()

inits<- function() list(beta=rnorm(3))
params<- c("beta", "sd.site","Npred")
ni<- 12000 ; nt<- 10 ; nb <- 2000 ; nc<- 3

out.500<- jags(win.data, params,inits = inits, "meta.analysis.500_revision.txt", n.chains = nc,
           n.thin = nt, n.iter = ni, n.burnin = nb)


str(win.data<- list(is=is600,N=N.pm.600,psd=N.psd.600, n= length(N.pm.600), pred.is=is.pred, npred=length(is.pred)))

sink("meta.analysis.600_revision.txt")
cat("
    model{
    
    #Priors
    for (v in 1:3) {     #previas for intercep and plonomial coefficients
    beta[v] ~dnorm (0,0.0001)
    }
    tau.site<-pow(sd.site, -2)
    sd.site~dunif(0,10)
    
    # Likelihood
    for (i in 1:n) {
    N[i]~dnorm(muN[i], tau.psd[i])
    tau.psd[i]<-pow(psd[i],-2)
    muN[i]<- beta[1]+ beta[2]*is[i] +beta[3]*pow(is[i],2) + eps.site[i]
    eps.site [i] ~ dnorm (0, tau.site)
    }
    
    #get predictions for plot
    for (i in 1:npred){
    Npred[i] <- beta[1] + beta[2]*pred.is[i] + beta[3]*pow(pred.is[i],2) 
    
    }
    }
    ", fill= TRUE)
sink()

inits<- function() list(beta=rnorm(3))
params<- c("beta", "sd.site","Npred")
ni<- 12000 ; nt<- 10 ; nb <- 2000 ; nc<- 3

out.600<- jags(win.data, params,inits = inits, "meta.analysis.600_revision.txt", n.chains = nc,
           n.thin = nt, n.iter = ni, n.burnin = nb)





plot(is.150.pred, out$mean$Npred, col="blue", ylim = c(0,6.5),type = "l")

points(is200, N.pm.200, col= "lightgoldenrod4", pch = 15,cex= 0.8)
points(is300, N.pm.300, col= "deeppink4", pch = 13,cex= 0.8)
points(is400, N.pm.400, col= "red2", pch = 17,cex= 0.8)
points(is500, N.pm.500, col= "darkcyan", pch = 18,cex= 0.8)
points(is600, N.pm.600, col= "darkorange", pch = 19,cex= 0.8)
lines(is.200.pred, out.200$mean$Npred, col="red")
lines(is.300.pred, out.300$mean$Npred, col="green")
lines(is.400.pred, out.400$mean$Npred, col="grey")
lines(is.500.pred, out.500$mean$Npred, col="brown")
lines(is.600.pred, out.600$mean$Npred, col="yellow")




tiff(file = "fig 5b_nuevo_estimada.tiff",                #Guardados como Tiff
     width = 84, height = 70,
     units = "mm", res =1200)

par(mai=c(0.35,0.4,0.1,0.3),xpd=F)
#plot estimares as a function of is
plot(is150, N.pm.150,  ylab= "Estimated species richness", 
     ylim = c(0,6.5), xlim=c(-1.8,2.7),pch= 8, bty="l",cex= 0.8,cex.lab=0.9,
     xaxt="n",yaxt="n", col="grey10")
axis(2,at=0:6,ylab= "Estimated species richness",labels = c(seq(0,6.5,1)), las=1,
     cex.axis=0.65,tck=-0.02, mgp=c(3,0.4,0))
axis(1,at=-2:2.8,labels = c(seq(-2,2.8)), las=1,
     cex.axis=0.6,tck=-0.02, mgp=c(3,-0.1,0))
mtext(1,line=0.55,text="Shannon Habitat Index",cex = 0.65)
mtext(2,line=1,text="Estimated species richness",cex = 0.65)

points(is200, N.pm.200, col= "lightgoldenrod4", pch = 15,cex= 0.8)
points(is300, N.pm.300, col= "deeppink4", pch = 13,cex= 0.8)
points(is400, N.pm.400, col= "red2", pch = 17,cex= 0.8)
points(is500, N.pm.500, col= "darkcyan", pch = 18,cex= 0.8)
points(is600, N.pm.600, col= "darkorange", pch = 19,cex= 0.8)

lines(is.pred, out.150$mean$Npred, col="grey10", bty="l", lwd=1.2)
lines(is.pred, out.200$mean$Npred, col="lightgoldenrod4", lty=2, lwd=1.2)
lines(is.pred, out.300$mean$Npred, col="deeppink4", lty=3, lwd=1.2)
lines(is.pred, out.400$mean$Npred, col="red2", lty=4, lwd=1.2)
lines(is.pred, out.500$mean$Npred, col="darkcyan", lty=5, lwd=1.2)
lines(is.pred, out.600$mean$Npred, col="darkorange", lty=6, lwd=1.2)
# 
# 
# mtext("150 m",4,las=1,cex = 0.3, at=3.75, line=0.35)
# mtext("200 m",4,las=1,cex = 0.3, at=3.6, line=0.35)
# mtext("600 m",4,las=1,cex = 0.3, at=3.45, line=0.35)
# mtext("300 m",4,las=1, cex = 0.3, at=3.3, line=0.35)
# mtext("500 m",4,las=1,cex = 0.3, at=3.15, line=0.35)
# mtext("400 m",4,las=1,cex = 0.3, at=3, line=0.35)
# # 
# par(mai=c(1,1,0.3,0.6),xpd=T)
# segments(2.88,3.45,3,3.75, lwd=0.5)  #150
# segments(2.88,3.35,3,3.6, lwd=0.5)  #200
# segments(2.88,3.28,3,3.45, lwd=0.5)  #600
# segments(2.88,3.22,3,3.3, lwd=0.5)  #300
# segments(2.88,3.06,3,3.15, lwd=0.5)   #500
# segments(2.88,3.01,3,3, lwd=0.5)   #400

legend(1.2,6.5, c("150 m", "200 m", "300 m", "400 m", "500 m", "600 m"),
       col = c("black", "lightgoldenrod4", "deeppink4", "red2", "darkcyan","darkorange"),
       lty = c(1,2,3,4,5,6), pch = c(8,15,13,17,18,19),
       pt.cex = c(0.6,0.6,0.6,0.6,0.7,0.6),
       ncol =2,y.intersp = 0.75,merge = F, cex = 0.4, box.lwd=0.7)

dev.off()


#####################



###############################
# Miodel Regression para Riqueza e IS
pa150 <- (scale(pa150[1:50]))
pa200 <- (scale(pa200[1:50]))
pa300 <- (scale(pa300[1:50]))
pa400 <- (scale(pa400[1:50]))
pa500 <- (scale(pa500[1:50]))
pa600 <- (scale(pa600[1:50]))


#promedio de riqueza y intervalos para todas las escalas
N.pm.150 <-(apply(out150.b$mean$Nsite,1,mean))
N.psd.150 <- apply(out150.is.b$sd$Nsite,1,mean)
N.cri.150<- apply(out150.is.b$sims.list$Nsite, 2,function(x) quantile(x,prob=c(0.125,0.875)))

N.pm.200 <-(apply(out200.is.b$mean$Nsite,1,mean))
N.psd.200 <- apply(out200.is.b$sd$Nsite,1,mean)
N.cri.200<- apply(out200.is.b$sims.list$Nsite, 2,function(x) quantile(x,prob=c(0.125,0.875)))

N.pm.300 <-(apply(out300.is.b$mean$Nsite,1,mean))
N.psd.300 <- apply(out300.is.b$sd$Nsite,1,mean)
N.cri.300<- apply(out300.is.b$sims.list$Nsite, 2,function(x) quantile(x,prob=c(0.125,0.875)))

N.pm.400 <-(apply(out400.is.b$mean$Nsite,1,mean))
N.psd.400 <- apply(out400.is.b$sd$Nsite,1,mean)
N.cri.400<- apply(out400.is.b$sims.list$Nsite, 2,function(x) quantile(x,prob=c(0.125,0.875)))

N.pm.500 <-(apply(out500.is.b$mean$Nsite,1,mean))
N.psd.500 <- apply(out500.is.b$sd$Nsite,1,mean)
N.cri.500<- apply(out500.is.b$sims.list$Nsite, 2,function(x) quantile(x,prob=c(0.125,0.875)))

N.pm.600 <-(apply(out600.is.b$mean$Nsite,1,mean))
N.psd.600 <- apply(out600.is.b$sd$Nsite,1,mean)
N.cri.600<- apply(out600.is.b$sims.list$Nsite, 2,function(x) quantile(x,prob=c(0.125,0.875)))


is.pred<- seq(-1.436,2.6,length=500)

str(win.data<- list(is=is150,N=N.pm.150,psd=N.psd.150, n= length(N.pm.150), pred.is=is.pred, npred=length(is.pred)))

sink("meta.analysis.150_revision.txt")
cat("
    model{
    
    #Priors
    for (v in 1:3) {     #previas for intercep and plonomial coefficients
    beta[v] ~dnorm (0,0.0001)
    }
    tau.site<-pow(sd.site, -2)
    sd.site~dunif(0,10)
    
    # Likelihood
    for (i in 1:n) {
    N[i]~dnorm(muN[i], tau.psd[i])
    tau.psd[i]<-pow(psd[i],-2)
    muN[i]<- beta[1]+ beta[2]*is[i] +beta[3]*pow(is[i],2) + eps.site[i]
    eps.site [i] ~ dnorm (0, tau.site)
    }
    
    #get predictions for plot
    for (i in 1:npred){
    Npred[i] <- beta[1] + beta[2]*pred.is[i] + beta[3]*pow(pred.is[i],2) 
    
    }
    }
    ", fill= TRUE)
sink()

inits<- function() list(beta=rnorm(3))
params<- c("beta", "sd.site","Npred")
ni<- 12000 ; nt<- 10 ; nb <- 2000 ; nc<- 3

out.150<- jags(win.data, params,inits = inits, "meta.analysis.150_revision.txt", n.chains = nc,
               n.thin = nt, n.iter = ni, n.burnin = nb)


str(win.data<- list(is=is200,N=N.pm.200,psd=N.psd.200, n= length(N.pm.200), pred.is=is.pred, npred=length(is.pred)))

sink("meta.analysis.200_revision.txt")
cat("
    model{
    
    #Priors
    for (v in 1:3) {     #previas for intercep and plonomial coefficients
    beta[v] ~dnorm (0,0.0001)
    }
    tau.site<-pow(sd.site, -2)
    sd.site~dunif(0,10)
    
    # Likelihood
    for (i in 1:n) {
    N[i]~dnorm(muN[i], tau.psd[i])
    tau.psd[i]<-pow(psd[i],-2)
    muN[i]<- beta[1]+ beta[2]*is[i] +beta[3]*pow(is[i],2) + eps.site[i]
    eps.site [i] ~ dnorm (0, tau.site)
    }
    
    #get predictions for plot
    for (i in 1:npred){
    Npred[i] <- beta[1] + beta[2]*pred.is[i] + beta[3]*pow(pred.is[i],2) 
    
    }
    }
    ", fill= TRUE)
sink()

inits<- function() list(beta=rnorm(3))
params<- c("beta", "sd.site","Npred")
ni<- 12000 ; nt<- 10 ; nb <- 2000 ; nc<- 3

out.200<- jags(win.data, params,inits = inits, "meta.analysis.200_revision.txt", n.chains = nc,
               n.thin = nt, n.iter = ni, n.burnin = nb)


str(win.data<- list(is=is300,N=N.pm.300,psd=N.psd.300, n= length(N.pm.300), pred.is=is.pred, npred=length(is.pred)))

sink("meta.analysis.300_revision.txt")
cat("
    model{
    
    #Priors
    for (v in 1:3) {     #previas for intercep and plonomial coefficients
    beta[v] ~dnorm (0,0.0001)
    }
    tau.site<-pow(sd.site, -2)
    sd.site~dunif(0,10)
    
    # Likelihood
    for (i in 1:n) {
    N[i]~dnorm(muN[i], tau.psd[i])
    tau.psd[i]<-pow(psd[i],-2)
    muN[i]<- beta[1]+ beta[2]*is[i] +beta[3]*pow(is[i],2) + eps.site[i]
    eps.site [i] ~ dnorm (0, tau.site)
    }
    
    #get predictions for plot
    for (i in 1:npred){
    Npred[i] <- beta[1] + beta[2]*pred.is[i] + beta[3]*pow(pred.is[i],2) 
    
    }
    }
    ", fill= TRUE)
sink()

inits<- function() list(beta=rnorm(3))
params<- c("beta", "sd.site","Npred")
ni<- 12000 ; nt<- 10 ; nb <- 2000 ; nc<- 3

out.300<- jags(win.data, params,inits = inits, "meta.analysis.300_revision.txt", n.chains = nc,
               n.thin = nt, n.iter = ni, n.burnin = nb)


str(win.data<- list(is=is400,N=N.pm.400,psd=N.psd.400, n= length(N.pm.400), pred.is=is.pred, npred=length(is.pred)))

sink("meta.analysis.400_revision.txt")
cat("
    model{
    
    #Priors
    for (v in 1:3) {     #previas for intercep and plonomial coefficients
    beta[v] ~dnorm (0,0.0001)
    }
    tau.site<-pow(sd.site, -2)
    sd.site~dunif(0,10)
    
    # Likelihood
    for (i in 1:n) {
    N[i]~dnorm(muN[i], tau.psd[i])
    tau.psd[i]<-pow(psd[i],-2)
    muN[i]<- beta[1]+ beta[2]*is[i] +beta[3]*pow(is[i],2) + eps.site[i]
    eps.site [i] ~ dnorm (0, tau.site)
    }
    
    #get predictions for plot
    for (i in 1:npred){
    Npred[i] <- beta[1] + beta[2]*pred.is[i] + beta[3]*pow(pred.is[i],2) 
    
    }
    }
    ", fill= TRUE)
sink()

inits<- function() list(beta=rnorm(3))
params<- c("beta", "sd.site","Npred")
ni<- 12000 ; nt<- 10 ; nb <- 2000 ; nc<- 3

out.400<- jags(win.data, params,inits = inits, "meta.analysis.400_revision.txt", n.chains = nc,
               n.thin = nt, n.iter = ni, n.burnin = nb)



str(win.data<- list(is=is500,N=N.pm.500,psd=N.psd.500, n= length(N.pm.500), pred.is=is.pred, npred=length(is.pred)))

sink("meta.analysis.500_revision.txt")
cat("
    model{
    
    #Priors
    for (v in 1:3) {     #previas for intercep and plonomial coefficients
    beta[v] ~dnorm (0,0.0001)
    }
    tau.site<-pow(sd.site, -2)
    sd.site~dunif(0,10)
    
    # Likelihood
    for (i in 1:n) {
    N[i]~dnorm(muN[i], tau.psd[i])
    tau.psd[i]<-pow(psd[i],-2)
    muN[i]<- beta[1]+ beta[2]*is[i] +beta[3]*pow(is[i],2) + eps.site[i]
    eps.site [i] ~ dnorm (0, tau.site)
    }
    
    #get predictions for plot
    for (i in 1:npred){
    Npred[i] <- beta[1] + beta[2]*pred.is[i] + beta[3]*pow(pred.is[i],2) 
    
    }
    }
    ", fill= TRUE)
sink()

inits<- function() list(beta=rnorm(3))
params<- c("beta", "sd.site","Npred")
ni<- 12000 ; nt<- 10 ; nb <- 2000 ; nc<- 3

out.500<- jags(win.data, params,inits = inits, "meta.analysis.500_revision.txt", n.chains = nc,
               n.thin = nt, n.iter = ni, n.burnin = nb)


str(win.data<- list(is=is600,N=N.pm.600,psd=N.psd.600, n= length(N.pm.600), pred.is=is.pred, npred=length(is.pred)))

sink("meta.analysis.600_revision.txt")
cat("
    model{
    
    #Priors
    for (v in 1:3) {     #previas for intercep and plonomial coefficients
    beta[v] ~dnorm (0,0.0001)
    }
    tau.site<-pow(sd.site, -2)
    sd.site~dunif(0,10)
    
    # Likelihood
    for (i in 1:n) {
    N[i]~dnorm(muN[i], tau.psd[i])
    tau.psd[i]<-pow(psd[i],-2)
    muN[i]<- beta[1]+ beta[2]*is[i] +beta[3]*pow(is[i],2) + eps.site[i]
    eps.site [i] ~ dnorm (0, tau.site)
    }
    
    #get predictions for plot
    for (i in 1:npred){
    Npred[i] <- beta[1] + beta[2]*pred.is[i] + beta[3]*pow(pred.is[i],2) 
    
    }
    }
    ", fill= TRUE)
sink()

inits<- function() list(beta=rnorm(3))
params<- c("beta", "sd.site","Npred")
ni<- 12000 ; nt<- 10 ; nb <- 2000 ; nc<- 3

out.600<- jags(win.data, params,inits = inits, "meta.analysis.600_revision.txt", n.chains = nc,
               n.thin = nt, n.iter = ni, n.burnin = nb)





plot(is.150.pred, out$mean$Npred, col="blue", ylim = c(0,6.5),type = "l")

points(is200, N.pm.200, col= "lightgoldenrod4", pch = 15,cex= 0.8)
points(is300, N.pm.300, col= "deeppink4", pch = 13,cex= 0.8)
points(is400, N.pm.400, col= "red2", pch = 17,cex= 0.8)
points(is500, N.pm.500, col= "darkcyan", pch = 18,cex= 0.8)
points(is600, N.pm.600, col= "darkorange", pch = 19,cex= 0.8)
lines(is.200.pred, out.200$mean$Npred, col="red")
lines(is.300.pred, out.300$mean$Npred, col="green")
lines(is.400.pred, out.400$mean$Npred, col="grey")
lines(is.500.pred, out.500$mean$Npred, col="brown")
lines(is.600.pred, out.600$mean$Npred, col="yellow")




tiff(file = "fig 5b_nuevo_estimada.tiff",                #Guardados como Tiff
     width = 84, height = 70,
     units = "mm", res =1200)

par(mai=c(0.35,0.4,0.1,0.3),xpd=F)
#plot estimares as a function of is
plot(is150, N.pm.150,  ylab= "Estimated species richness", 
     ylim = c(0,6.5), xlim=c(-1.8,2.7),pch= 8, bty="l",cex= 0.8,cex.lab=0.9,
     xaxt="n",yaxt="n", col="grey10")
axis(2,at=0:6,ylab= "Estimated species richness",labels = c(seq(0,6.5,1)), las=1,
     cex.axis=0.65,tck=-0.02, mgp=c(3,0.4,0))
axis(1,at=-2:2.8,labels = c(seq(-2,2.8)), las=1,
     cex.axis=0.6,tck=-0.02, mgp=c(3,-0.1,0))
mtext(1,line=0.55,text="Shannon Habitat Index",cex = 0.65)
mtext(2,line=1,text="Estimated species richness",cex = 0.65)

points(is200, N.pm.200, col= "lightgoldenrod4", pch = 15,cex= 0.8)
points(is300, N.pm.300, col= "deeppink4", pch = 13,cex= 0.8)
points(is400, N.pm.400, col= "red2", pch = 17,cex= 0.8)
points(is500, N.pm.500, col= "darkcyan", pch = 18,cex= 0.8)
points(is600, N.pm.600, col= "darkorange", pch = 19,cex= 0.8)

lines(is.pred, out.150$mean$Npred, col="grey10", bty="l", lwd=1.2)
lines(is.pred, out.200$mean$Npred, col="lightgoldenrod4", lty=2, lwd=1.2)
lines(is.pred, out.300$mean$Npred, col="deeppink4", lty=3, lwd=1.2)
lines(is.pred, out.400$mean$Npred, col="red2", lty=4, lwd=1.2)
lines(is.pred, out.500$mean$Npred, col="darkcyan", lty=5, lwd=1.2)
lines(is.pred, out.600$mean$Npred, col="darkorange", lty=6, lwd=1.2)
# 
# 
# mtext("150 m",4,las=1,cex = 0.3, at=3.75, line=0.35)
# mtext("200 m",4,las=1,cex = 0.3, at=3.6, line=0.35)
# mtext("600 m",4,las=1,cex = 0.3, at=3.45, line=0.35)
# mtext("300 m",4,las=1, cex = 0.3, at=3.3, line=0.35)
# mtext("500 m",4,las=1,cex = 0.3, at=3.15, line=0.35)
# mtext("400 m",4,las=1,cex = 0.3, at=3, line=0.35)
# # 
# par(mai=c(1,1,0.3,0.6),xpd=T)
# segments(2.88,3.45,3,3.75, lwd=0.5)  #150
# segments(2.88,3.35,3,3.6, lwd=0.5)  #200
# segments(2.88,3.28,3,3.45, lwd=0.5)  #600
# segments(2.88,3.22,3,3.3, lwd=0.5)  #300
# segments(2.88,3.06,3,3.15, lwd=0.5)   #500
# segments(2.88,3.01,3,3, lwd=0.5)   #400

legend(1.2,6.5, c("150 m", "200 m", "300 m", "400 m", "500 m", "600 m"),
       col = c("black", "lightgoldenrod4", "deeppink4", "red2", "darkcyan","darkorange"),
       lty = c(1,2,3,4,5,6), pch = c(8,15,13,17,18,19),
       pt.cex = c(0.6,0.6,0.6,0.6,0.7,0.6),
       ncol =2,y.intersp = 0.75,merge = F, cex = 0.4, box.lwd=0.7)

dev.off()

