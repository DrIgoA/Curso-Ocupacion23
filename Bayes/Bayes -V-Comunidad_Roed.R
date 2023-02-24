############################################################
##### CURSO METODOS CUANTITATIVOS DETECCION IMPERFECTA #####
###########        UNRC - JUNIO JULIO 2016       ###########
###########             A. P. GOIJMAN            ###########
############################################################
###########       Ejemplo Ocupacion en JAGS      ###########
###########              Comunidades             ###########
############################################################

######################################
### Borrar documentacion anterior de R
######################################
rm(list=ls(all=TRUE))

###################################### 
### Elegir directorio de trabajo
######################################
setwd("C:/Users/Andrea/Dropbox/TEACHING/APG_MetCuant_UNRC_Jun2016/Docentes/Ejercicios/Bayes")

######################################
### Cargar Paquetes 
######################################
library(jagsUI)    #paquete JAGS

######################################
### Cargar Datos Roedores 
######################################

#load("C:\\Users\\Andrea\\Dropbox\\TEACHING\\APG_MetCuant_UNRC_Jun2016\\Docentes\\Ejercicios\\Bayes\\Roedores.RData")
ls()   #revisar 

###Matrices de covariables 
#manejo
#uso
#vvegz
#tempz
#y
#nsite


##########################################################################
##### Model - Random effects community occupancy model####################
######Allows for estimation over less frequent species####################
############    sin DA + tiempo random + covariates y p and psi###########
##########################################################################

# Armar paquete de datos y resumirlos

str(roed.data <- list(y=y, manejo=manejo, uso=uso, vvegz=vvegz, tempz=tempz, nsite=nsite,
nrep=dim(y)[2], nspec=dim(y)[3], season=dim(y)[4]))

### JAGS code
sink("roedor3.jags")  
cat("
model{

# modelos para missing covariates
 for(j in 1:134){                                   #por sitio i
  for (l in 1:3){
    uso[j,l]~dbern(puso)
    vvegz[j,l]~dnorm(mu.vvegz, tau.vvegz)
    tempz[j,l]~dnorm(mu.tempz, tau.tempz)
    }}
  puso~dunif(0,1)
  tau.vvegz<-pow(sd.vvegz,-2)
  sd.vvegz~dunif(0,100)
  mu.vvegz~dnorm(0,0.0001)
  tau.tempz<-pow(sd.tempz,-2)
  sd.tempz~dunif(0,100)
  mu.tempz~dnorm(0,0.0001)

# Hyperpriors 
 mu.co.lpsi ~ dnorm(0, 0.01)
 mu.co.betalpsi1 ~ dnorm(0, 0.01)
 mu.co.betalpsi2 ~ dnorm(0, 0.01)
 mu.co.betalpsi3 ~ dnorm(0, 0.01)
 mu.co.lp ~ dnorm(0, 0.01)
 mu.co.betalp1 ~ dnorm(0, 0.01)
 tau.co.lpsi <- pow(sd.co.lpsi,-2)   
 sd.co.lpsi ~ dunif(0,8)
 tau.co.lp <- pow(sd.co.lp,-2)   
 sd.co.lp ~ dunif(0,8)
 tau.co.betalpsi1 <- pow(sd.co.betalpsi1,-2)
 sd.co.betalpsi1 ~ dunif(0,8)
 tau.co.betalpsi2 <- pow(sd.co.betalpsi2,-2)
 sd.co.betalpsi2 ~ dunif(0,8)
 tau.co.betalpsi3 <- pow(sd.co.betalpsi3,-2)
 sd.co.betalpsi3 ~ dunif(0,8)
 tau.co.betalp1 <- pow(sd.co.betalp1,-2)
 sd.co.betalp1 ~ dunif(0,8)

# priors for species specific effects in occupancy and detection 
# sp represents per species
# co represents community 
  for (k in 1:nspec){
    mu.sp.lpsi[k] ~ dnorm(mu.co.lpsi, tau.co.lpsi)  
    mu.sp.betalpsi1[k]~ dnorm(mu.co.betalpsi1, tau.co.betalpsi1)
    mu.sp.betalpsi2[k]~ dnorm(mu.co.betalpsi2, tau.co.betalpsi2)
    mu.sp.betalpsi3[k]~ dnorm(mu.co.betalpsi3, tau.co.betalpsi3)
    mu.sp.lp[k] ~ dnorm(mu.co.lp, tau.co.lp)
    mu.sp.betalp1[k]~ dnorm(mu.co.betalp1, tau.co.betalp1)
    tau.sp.lpsi[k] <- pow(sd.sp.lpsi[k],-2)   
    sd.sp.lpsi[k] ~ dunif(0,8)
    tau.sp.lp[k] <- pow(sd.sp.lp[k],-2)
    sd.sp.lp[k] ~ dunif(0,8)
    tau.sp.betalpsi1[k] <- pow(sd.sp.betalpsi1[k],-2)
    sd.sp.betalpsi1[k] ~ dunif(0,8)
    tau.sp.betalpsi2[k] <- pow(sd.sp.betalpsi2[k],-2)
    sd.sp.betalpsi2[k] ~ dunif(0,8)
    tau.sp.betalpsi3[k] <- pow(sd.sp.betalpsi3[k],-2)
    sd.sp.betalpsi3[k] ~ dunif(0,8)
    tau.sp.betalp1[k] <- pow(sd.sp.betalp1[k],-2)
    sd.sp.betalp1[k] ~ dunif(0,8)

# priors for species-specific time effects
  for (l in 1:season){ 
      lpsi[k,l] ~ dnorm(mu.sp.lpsi[k], tau.sp.lpsi[k])
      lp[k,l] ~ dnorm(mu.sp.lp[k], tau.sp.lp[k])
      betalpsi1[k,l] ~ dnorm (mu.sp.betalpsi1[k], tau.sp.betalpsi1[k])
      betalpsi2[k,l] ~ dnorm (mu.sp.betalpsi2[k], tau.sp.betalpsi2[k])
      betalpsi3[k,l] ~ dnorm (mu.sp.betalpsi3[k], tau.sp.betalpsi3[k])
      betalp1[k,l] ~ dnorm(mu.sp.betalp1[k],tau.sp.betalp1[k])
}
  
# Ecological model, process model (true occurrence at site i) 
  ## time as random effect
     
      for (i in 1:nsite) {                                           #loop sobre sitios 
        for (l in 1:season){                                         #loop sobre tiempo 
              logit(psi[i,k,l]) <- lpsi[k,l] + betalpsi1[k,l]*vvegz[i,l] + 
                    betalpsi2[k,l]*uso[i,l]+ betalpsi3[k,l]*manejo[i,l] 
              z[i,k,l] ~ dbern(psi[i,k,l])
   
# Observation model for site i, replicate nrep=j, nspec=k, season=l

            for (j in 1:nrep) {                                     #loop sobre replicates nrep   
               logit(p[i,j,k,l]) <-  lp[k,l] + betalp1[k,l]* tempz[i,l]
               mu.p[i,j,k,l] <- p[i,j,k,l]*z[i,k,l]
               y[i,j,k,l] ~ dbern(mu.p[i,j,k,l])
               } #nrep
           } #season
         } #nsite
       } #nspec

# Derived quantities

 for (i in 1:nsite) {                                           
   for (l in 1:season){                                         
    Nsite [i,l] <- sum (z[i,,])
   }
  }

} #model
",fill=TRUE)
sink()

# valores iniciales
  zst <- apply(y,c(1,3,4),max)
  zst[is.na(zst)]<-1
  inits <- function() list(z=zst,mu.sp.lpsi=rnorm(roed.data$nspec), mu.sp.betalpsi1=rnorm(roed.data$nspec),
  mu.sp.betalpsi2=rnorm(roed.data$nspec),mu.sp.betalpsi3=rnorm(roed.data$nspec),
  mu.sp.lp=rnorm(roed.data$nspec),mu.sp.betalp1=rnorm(roed.data$nspec),sd.sp.lpsi=runif(roed.data$nspec,0.1,5),
  sd.sp.lp=runif(roed.data$nspec,0.1,5),sd.sp.betalpsi1=runif(roed.data$nspec,0.1,5),sd.sp.betalpsi2=runif(roed.data$nspec,0.1,5),
  sd.sp.betalpsi3=runif(roed.data$nspec,0.1,5),sd.sp.betalp1=runif(roed.data$nspec,0.1,5))

# parametros a monitorear
  # species-specific
  params1 <- c("lp", "betalp1",  "lpsi", "betalpsi1", "betalpsi2", "betalpsi3", 
               "mu.sp.lpsi", "mu.sp.betalpsi1", "mu.sp.betalpsi2", "mu.sp.betalpsi3", "mu.sp.lp", "mu.sp.betalp1",
               "sd.sp.lpsi","sd.sp.lp", "sd.sp.betalpsi1", "sd.sp.betalpsi2","sd.sp.betalpsi3","sd.sp.betalp1")
  #community
   params2 <- c( "mu.co.lpsi", "mu.co.betalpsi1", "mu.co.betalpsi2", "mu.co.betalpsi3", "mu.co.lp", 
                 "mu.co.betalp1", "sd.co.lpsi", "sd.co.lp", "sd.co.betalpsi1","sd.co.betalpsi2", "sd.co.betalpsi3",
                  "sd.co.betalp1")
   params3 <- "Nsite"

# ajustes de MCMC
  ni <- 30000
  nt <- 10 
  nb <- 1000
  nc <- 3

# llamar JAGS desde R 
  library(jagsUI)    #paquete JAGS

out41 = jags(roed.data, inits, params1, "roedor3.jags", n.chains=nc, 
            n.iter=ni, n.burnin=nb, n.thin=nt, parallel=TRUE)
out42 = jags(roed.data, inits, params2, "roedor3.jags", n.chains=nc, 
            n.iter=ni, n.burnin=nb, n.thin=nt, parallel=TRUE)
out4N = jags(roed.data, inits, params3, "roedor3.jags", n.chains=nc, 
            n.iter=ni, n.burnin=nb, n.thin=nt, parallel=TRUE)

par(mfrow=c(3,2))
traceplot(out31)
traceplot(out41)

print(out41,dig=2)



ls()












o.temp <- seq(0,30,,500)
o.vveg <- seq(0,1.05,,500)
temp.pred <- (o.temp - mean(o.temp))/sd(o.temp)
vveg.pred <- (o.vveg - mean(o.vveg))/sd(o.vveg)

#Litado de especies para graficar
sp<-c('Aa', 'Cl', 'Cm', 'Cv', 'Mm', 'Of', 'Or')
season<-c('Spring','Summer','Fall')

summary(temp.pred)

#posterior means por especie and season

# DETECCION
lp <- apply(out41$sims.list$lp,c(2,3),mean)
cri.lp <- apply(out41$sims.list$lp,c(2,3),function(x) quantile(x,prob=c(0.025,0.975)))
betalp1 <- apply(out41$sims.list$betalp1,c(2,3),mean)
cri.betalp1 <- apply(out41$sims.list$betalp1,c(2,3),function(x) quantile(x,prob=c(0.025,0.975)))

# Plot covariate effects on species detection
par(mfrow=c(1,1), cex=1, cex.lab=1.2, cex.axis=0.9)
plot(betalp1,1:21, xlim=c(-7,7),xlab='Temperature', ylab='Species',yaxt="n",,pch=16,tck=-0.02)
axis(2,c(4.5,11.5,18.5),labels=season, xaxt="n", padj=-1.4, tick=FALSE) 
axis(2,1:21,labels=c(sp,sp,sp),xaxt="n", las=2, hadj=0.5, tck=-0.02) 
abline(a=7.5, b=0, lwd=1,,lty=2, col='black')
abline(a=14.5, b=0, lwd=1,,lty=2, col='black')
abline(v=0, lwd=2, col='black')
segments(cri.betalp1[1,,],1:21,cri.betalp1[2,,],1:21,col='grey', lwd=1)
sig2<- cri.betalp1[1,,]*cri.betalp1[2,,]>0    #mayor a cero son los dos + o dos - == efecto signif   
segments(cri.betalp1[1,,][sig2==1],(1:21)[sig2==1],cri.betalp1[2,,][sig2==1],(1:21)[sig2==1],lwd=1, col='red')
abline(v=out42$summary[6,1],lwd=1,col='red')
abline(v=out42$summary[6,c(3,7)],lwd=1,col='red',lty=2)

# OCUPACION
lpsi <- apply(out41$sims.list$lpsi,c(2,3),mean)
cri.lpsi <- apply(out41$sims.list$lpsi,c(2,3),function(x) quantile(x,prob=c(0.025,0.975)))
betalpsi1 <- apply(out41$sims.list$betalpsi1,c(2,3),mean)
cri.betalpsi1 <- apply(out41$sims.list$betalpsi1,c(2,3),function(x) quantile(x,prob=c(0.025,0.975)))
betalpsi2 <- apply(out41$sims.list$betalpsi2,c(2,3),mean)
cri.betalpsi2 <- apply(out41$sims.list$betalpsi2,c(2,3),function(x) quantile(x,prob=c(0.025,0.975)))
betalpsi3 <- apply(out41$sims.list$betalpsi3,c(2,3),mean)
cri.betalpsi3 <- apply(out41$sims.list$betalpsi3,c(2,3),function(x) quantile(x,prob=c(0.025,0.975)))

#Betalpsi1<-c(betalpsi1[,1],betalpsi1[,2],betalpsi1[,3])

# Plot covariate effects on species occupancy
par(mfrow=c(1,1), cex=1, cex.lab=1.2, cex.axis=0.9)

# Vegeteation Volume
plot(betalpsi1,1:21, xlim=c(-5,15),xlab='Vegetation Vol.', ylab='Species',yaxt="n",,pch=16,tck=-0.02)
axis(2,c(4.5,11.5,18.5),labels=season, xaxt="n", padj=-1.4, tick=FALSE) 
axis(2,1:21,labels=c(sp,sp,sp),xaxt="n", las=2, hadj=0.5, tck=-0.02) 
abline(a=7.5, b=0, lwd=1,,lty=2, col='black')
abline(a=14.5, b=0, lwd=1,,lty=2, col='black')
abline(v=0, lwd=2, col='black')
segments(cri.betalpsi1[1,,],1:21,cri.betalpsi1[2,,],1:21,col='grey', lwd=1)
sig2<- cri.betalpsi1[1,,]*cri.betalpsi1[2,,]>0    #mayor a cero son los dos + o dos - == efecto signif   
segments(cri.betalpsi1[1,,][sig2==1],(1:21)[sig2==1],cri.betalpsi1[2,,][sig2==1],(1:21)[sig2==1],lwd=1, col='red')
abline(v=out42$summary[2,1],lwd=1,col='red')
abline(v=out42$summary[2,c(3,7)],lwd=1,col='red',lty=2)

# Use 
plot(betalpsi2,1:21, xlim=c(-10,10),xlab='Use', ylab='Species',yaxt="n",,pch=16,tck=-0.02)
axis(2,c(4.5,11.5,18.5),labels=season, xaxt="n", padj=-1.4, tick=FALSE) 
axis(2,1:21,labels=c(sp,sp,sp),xaxt="n", las=2, hadj=0.5, tck=-0.02) 
abline(a=7.5, b=0, lwd=1,,lty=2, col='black')
abline(a=14.5, b=0, lwd=1,,lty=2, col='black')
abline(v=0, lwd=2, col='black')
segments(cri.betalpsi2[1,,],1:21,cri.betalpsi2[2,,],1:21,col='grey', lwd=1)
sig2<- cri.betalpsi2[1,,]*cri.betalpsi2[2,,]>0    #mayor a cero son los dos + o dos - == efecto signif   
segments(cri.betalpsi2[1,,][sig2==1],(1:21)[sig2==1],cri.betalpsi2[2,,][sig2==1],(1:21)[sig2==1],lwd=1, col='red')
abline(v=out42$summary[3,1],lwd=1,col='red')
abline(v=out42$summary[3,c(3,7)],lwd=1,col='red',lty=2)

# Organic
plot(betalpsi3,1:21, xlim=c(-9,9),xlab='Organic Management', ylab='Species',yaxt="n",,pch=16,tck=-0.02)
axis(2,c(4.5,11.5,18.5),labels=season, xaxt="n", padj=-1.4, tick=FALSE) 
axis(2,1:21,labels=c(sp,sp,sp),xaxt="n", las=2, hadj=0.5, tck=-0.02) 
abline(a=7.5, b=0, lwd=1,,lty=2, col='black')
abline(a=14.5, b=0, lwd=1,,lty=2, col='black')
abline(v=0, lwd=2, col='black')
segments(cri.betalpsi3[1,,],1:21,cri.betalpsi3[2,,],1:21,col='grey', lwd=1)
sig2<- cri.betalpsi3[1,,]*cri.betalpsi3[2,,]>0    #mayor a cero son los dos + o dos - == efecto signif   
segments(cri.betalpsi3[1,,][sig2==1],(1:21)[sig2==1],cri.betalpsi3[2,,][sig2==1],(1:21)[sig2==1],lwd=1, col='red')
abline(v=out42$summary[4,1],lwd=1,col='red')
abline(v=out42$summary[4,c(3,7)],lwd=1,col='red',lty=2)
