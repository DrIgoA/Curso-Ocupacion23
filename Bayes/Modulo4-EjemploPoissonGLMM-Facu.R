rm(list=ls())
library(jagsUI)

data <- read.csv("M4-GLMMPoisson-Facu.csv", header=T)
attach(data)
str(data)

plot(volveg,abund)

# sitios
site.list <- data$point2[1:128]
site.name <- seq(32)

# establecimientos
farm <- as.factor(data$farm)
farm <- levels(farm)

# año
year <- as.factor(data$year)
year <- levels(year)

# especies
spec <- as.factor(data$orden)
spec <- levels(spec)

nspec <- length(spec)
nsite <- length(site.name)
nfarm <- length(farm)
nyear <- length(year)

library(abind) 

### COVARIABLES
##  AÑO 1 
cova1 <- data[1:128,]

### Construir matrices individuales para cada covariable (sitio, año)
# FEBRERO separo 1ero x establecimiento
dh1 <- subset(cova1, cova1$farm == 'DH')
av1 <- subset(cova1, cova1$farm == 'AV')
lg1 <- subset(cova1, cova1$farm == 'LG')
ch1 <- subset(cova1, cova1$farm == 'CH')

# abundancia total de artrópodos
abund.dh1 <- dh1$abund 
abund.ch1 <- ch1$abund 
abund.lg1 <- lg1$abund 
abund.av1 <- av1$abund 

abund1 <- abind(abund.av1, abund.ch1, abund.dh1, abund.lg1, along = 2)
str(abund1)
dimnames(abund1) <- list(site=site.name, farm=farm)
abund1[,1]

# vveg - volumen vegetal
vveg.dh1 <- dh1$volveg 
vveg.dh1 <- as.matrix(scale(vveg.dh1))

vveg.ch1 <- ch1$volveg 
vveg.ch1 <- as.matrix(scale(vveg.ch1))

vveg.lg1 <- lg1$volveg 
vveg.lg1 <- as.matrix(scale(vveg.lg1))

vveg.av1 <- av1$volveg 
vveg.av1 <- as.matrix(scale(vveg.av1))

vveg1 <- abind(vveg.av1, vveg.ch1, vveg.dh1, vveg.lg1, along = 2)
str(vveg1)
dimnames(vveg1) <- list(site=site.name, farm=farm)
vveg1[,1]

#manejo = Orgánico(1)-Convencional(0)
manejo.dh1 <- as.matrix(dh1$org) 
manejo.av1 <- as.matrix(av1$org)
manejo.lg1 <- as.matrix(lg1$org)
manejo.ch1 <- as.matrix(ch1$org)

manejo1 <- abind(manejo.av1, manejo.ch1, manejo.dh1, manejo.lg1, along = 2)

dimnames(manejo1) <- list(site=site.name, estab=farm)
str(manejo1)

#TEMPORADA = VERANO(1)-PRIMAVERA(0)
season.dh1 <- as.matrix(dh1$ver) 
season.av1 <- as.matrix(av1$ver)
season.lg1 <- as.matrix(lg1$ver)
season.ch1 <- as.matrix(ch1$ver)

season1 <- abind(season.av1, season.ch1, season.dh1, season.lg1, along = 2)

dimnames(season1) <- list(site=site.name, estab=farm)
str(season1)

##  AÑO 2 
cova2 <- data[129:256,]

### Construir matrices individuales para cada covariable (sitio, año)
# FEBRERO separo 2ero x establecimiento
dh2 <- subset(cova2, cova2$farm == 'DH')
av2 <- subset(cova2, cova2$farm == 'AV')
lg2 <- subset(cova2, cova2$farm == 'LG')
ch2 <- subset(cova2, cova2$farm == 'CH')

# abundancia total de artrópodos
abund.dh2 <- dh2$abund 
abund.ch2 <- ch2$abund 
abund.lg2 <- lg2$abund 
abund.av2 <- av2$abund 

abund2 <- abind(abund.av2, abund.ch2, abund.dh2, abund.lg2, along = 2)
str(abund2)
dimnames(abund2) <- list(site=site.name, farm=farm)
abund2[,2]

# vveg - volumen vegetal
vveg.dh2 <- dh2$volveg 
vveg.dh2 <- as.matrix(scale(vveg.dh2))

vveg.ch2 <- ch2$volveg 
vveg.ch2 <- as.matrix(scale(vveg.ch2))

vveg.lg2 <- lg2$volveg 
vveg.lg2 <- as.matrix(scale(vveg.lg2))

vveg.av2 <- av2$volveg 
vveg.av2 <- as.matrix(scale(vveg.av2))

vveg2 <- abind(vveg.av2, vveg.ch2, vveg.dh2, vveg.lg2, along = 2)
str(vveg2)
dimnames(vveg2) <- list(site=site.name, farm=farm)
vveg2[,2]

#manejo = Orgánico(2)-Convencional(0)
manejo.dh2 <- as.matrix(dh2$org) 
manejo.av2 <- as.matrix(av2$org)
manejo.lg2 <- as.matrix(lg2$org)
manejo.ch2 <- as.matrix(ch2$org)

manejo2 <- abind(manejo.av2, manejo.ch2, manejo.dh2, manejo.lg2, along = 2)

dimnames(manejo2) <- list(site=site.name, estab=farm)
str(manejo2)

#TEMPORADA = VERANO(2)-PRIMAVERA(0)
season.dh2 <- as.matrix(dh2$ver) 
season.av2 <- as.matrix(av2$ver)
season.lg2 <- as.matrix(lg2$ver)
season.ch2 <- as.matrix(ch2$ver)

season2 <- abind(season.av2, season.ch2, season.dh2, season.lg2, along = 2)

dimnames(season2) <- list(site=site.name, estab=farm)
str(season2)

##########
abund <- abind(abund1,abund2, along = 3)
dimnames(abund) <- list(site=site.name, 
                        estab=farm,
                        year=year)

vveg <- abind(vveg1,vveg2, along = 3)
dimnames(vveg) <- list(site=site.name, 
                        estab=farm,
                        year=year)

manejo <- abind(manejo1,manejo2, along = 3)
dimnames(manejo) <- list(site=site.name, 
                       estab=farm,
                       year=year)

season <- abind(season1,season2, along = 3)
dimnames(season) <- list(site=site.name, 
                         estab=farm,
                         year=year)

###################################################################################################
# Modelo 2 = Volumen vegetal
win.data <- list(abund=abund, 
                 vveg=vveg,
                 manejo=manejo,
                 season=season,
                 nsite=dim(abund)[1],
                 nfarm=dim(abund)[2],
                 nyear=dim(abund)[3]
                 )
str(win.data)

# Specify model in BUGS language
cat(file = "M4-GLMMPoisson-Facu.txt","
model {
# modelos para missing covariates
    for(l in 1:nyear){
     for (e in 1:nfarm){ 
       for(j in 1:nsite) {
          
          vveg  [j,e,l] ~ dnorm(mu.vveg, tau.vveg)
          manejo[j,e,l] ~ dbern(pmanejo) 
          season[j,e,l] ~ dbern(pseason)
                          
                        }
                       }
                       }
                     
          tau.vveg <-pow(sd.vveg,-2)
          sd.vveg ~ dunif(0,1)
          mu.vveg ~ dnorm(0,1)
          
          pmanejo ~ dunif(0,1)
          
          pseason ~ dunif(0,1)
  
# Priors
   mu.alpha ~ dnorm(0, 0.001)  # Mean hyperparam
   tau.alpha <- pow(sd.alpha, -2)
   sd.alpha ~ dunif(0, 10)     # sd hyperparam
 
for(k in 1:5){
  alpha[k] ~ dunif(-10, 10)   # Regression params
}

# Likelihood
for(j in 1:nsite){
 
  alpha0[j] ~ dnorm(mu.alpha, tau.alpha) # Random effects and hyperparams
  re0[j] <- alpha0[j] - mu.alpha         # zero-centered random effects
  
for (e in 1:nfarm){ 
  for(l in 1:nyear) {   
 
  abund[j,e,l] ~ dpois(lambda[j,e,l])
  log(lambda[j,e,l]) <- alpha0[j] + alpha[1] * vveg  [j,e,l] 
                                  + alpha[2] * manejo[j,e,l]
                                  + alpha[3] * season[j,e,l]
                                  + alpha[4] * manejo[j,e,l]*vveg[j,e,l]
                                  + alpha[5] * manejo[j,e,l]*season[j,e,l]
                                  
 }
 }
  }
}")

# Other model run preparations
inits <- function() list(alpha0 = rnorm(1:32), alpha = rnorm(5)) # Inits
params <- c("mu.alpha", "sd.alpha", "alpha0", "alpha")           # Params
ni <- 300000 ; nt <- 25 ; nb <- 150000 ; nc <- 3                 # MCMC settings

# Call WinBUGS or JAGS from R (ART 6-7 min) and summarize posteriors

out <- jags(win.data, inits, params, "00_abund.txt", n.chains = nc, n.thin = nt,
             n.iter = ni, n.burnin = nb)

print(out, dig=3)
par(mfrow = c(3,2)) ; traceplot(out, c("mu.alpha", "sd.alpha", "alpha[1:5]"))

save(out, file='out.rda')

load('out.rda')

traceplot(out,c("mu.alpha","sd.alpha","alpha[1:5]"))

str(out)
# PREDICHOS
vveg.predo <- seq(0,1.97, , 20000) # Covariate values (VOLVEG) for prediction ORGÁNICO (MAX = 1.97)
vveg.predc <- seq(0,1.55, , 20000) # Covariate values (VOLVEG) for prediction CONVENTIONAL (MAX = 1.55)

# ORGÁNICO - PRIMAVERA
pred1 <- array(NA, dim = c(20000, 32))
for(i in 1:32){
  pred1[,i]<- exp(out$mean$alpha0[i] 
                + out$mean$alpha[1] * vveg.predo 
                + out$mean$alpha[2] * 1
                + out$mean$alpha[3] * 0
                #+ out$mean$alpha[4] * vveg.predo
                #+ out$mean$alpha[5] * vveg.predo * 0
                ) # Predictions for each site
}
# CONVENCIONAL - PRIMAVERA
pred2 <- array(NA, dim = c(20000, 32))
for(i in 1:32){
  pred2[,i]<- exp(out$mean$alpha0[i] 
                  + out$mean$alpha[1] * vveg.predc 
                  + out$mean$alpha[2] * 0 
                  + out$mean$alpha[3] * 0
                  #+ out$mean$alpha[4] * vveg.predc * 0
                  #+ out$mean$alpha[5] * vveg.predc * 0
                  ) # Predictions for each site
}
# ORGÁNICO - VERANO
pred3 <- array(NA, dim = c(20000, 32))
for(i in 1:32){
  pred3[,i]<- exp(out$mean$alpha0[i] 
                  + out$mean$alpha[1] * vveg.predo 
                  + out$mean$alpha[2] * 1 
                  + out$mean$alpha[3] * 1
                  #+ out$mean$alpha[4] * vveg.predo * 1
                  #+ out$mean$alpha[5] * vveg.predo * 1
                  ) # Predictions for each site
}
# CONVENCIONAL - VERANO
pred4 <- array(NA, dim = c(20000, 32))
for(i in 1:32){
  pred4[,i]<- exp(out$mean$alpha0[i] 
                  + out$mean$alpha[1] * vveg.predc 
                  + out$mean$alpha[2] * 0 
                  + out$mean$alpha[3] * 1
                  #+ out$mean$alpha[4] * vveg.predc * 0
                  #+ out$mean$alpha[5] * vveg.predc * 1
                  ) # Predictions for each site
}

cri<-function(x) quantile(x,prob=c(0.05,0.95))
#cri1<-apply(out$q2.5$alpha,1, mean)
#cri2<-apply(out$q97.5$alpha,1, mean)
cri1 <- apply(pred1[,], 1, cri)
cri2 <- apply(pred2[,], 1, cri)
cri3 <- apply(pred3[,], 1, cri)
cri4 <- apply(pred4[,], 1, cri)
str(cri1[1,])
str(cri1[2,])
str(vveg.predo)


### GRÁFICO 
tiff(file = "M4-GLMMPoisson-Facu.tiff",                #Guardados como Tiff
     width = 190, height = 150,   #width=140 height=150
     units = "mm", res = 300)
par(mfrow = c(1,2),mar=c(4,4,1,1))

# ORGÁNICO - PRIMAVERA
matplot(vveg.predo, cri1[1,], type = "l", lty = 0, lwd = 0, col = "blue", xlab = "",
        ylab = "Abundancia de artrópodos", frame.plot = F, ylim = c(0, 700)) # Fig. 5.17 (b)
mtext("Primavera",          at=1,side=1, line = -12, cex = 1, padj= -20.5, font=1)
mtext("(a)",             at=0.07,       line=-2,    cex = 1)
polygon(x= c(vveg.predo, rev(vveg.predo)), y= c(cri1[1,], rev(cri1[2,])), 
        col =  adjustcolor("green", alpha.f = 0.15), border = NA)

lines(vveg.predo, exp(out$mean$mu.alpha 
                      + out$mean$alpha[1] * vveg.predo
                      + out$mean$alpha[2] * 1
                      + out$mean$alpha[3] * 0
                      #+ out$mean$alpha[4] * vveg.predo
                      #+ out$mean$alpha[5] * vveg.predo * 0
                      ),
      col = "darkgreen", lwd = 5)
mtext("Volumen vegetal",       font=1, side=1, cex=1, at=1, line = 1, padj= 1.3, )

# CONVENCIONAL - PRIMAVERA
polygon(x= c(vveg.predc, rev(vveg.predc)), y= c(cri2[1,], rev(cri2[2,])), 
        col =  adjustcolor("magenta", alpha.f = 0.15), border = NA)
lines(vveg.predc, exp(out$mean$mu.alpha 
                      + out$mean$alpha[1] * vveg.predc 
                      + out$mean$alpha[2] * 0 
                      + out$mean$alpha[3] * 0
                      #+ out$mean$alpha[4] * vveg.predc * 0
                      #+ out$mean$alpha[5] * vveg.predc * 0
                      ),
      col = "darkmagenta", lwd = 5)

# ORGÁNICO - VERANO
matplot(vveg.predo, cri3[1,], type = "l", lty = 0, lwd = 0, col = "darkgreen", xlab = "",
        ylab = "", frame.plot = F, ylim = c(0, 250)) 
#mtext("Organic farming", at=1  ,side=1, line = -4.5, cex = 1, padj= -20.5, font=2)
mtext("Verano", at=0.8,side=1, line = -12, cex = 1, padj= -20.5, font=1)
mtext("(b)",    at=0.07,       line=-2,    cex = 1)
polygon(x= c(vveg.predo, rev(vveg.predo)), y= c(cri3[1,], rev(cri3[2,])), 
        col =  adjustcolor("green", alpha.f = 0.15), border = NA)
#lines(vveg.predo, cri3[2,], lty = 2, lwd = 2, col = "darkgreen")
lines(vveg.predo, exp(out$mean$mu.alpha 
                      + out$mean$alpha[1] * vveg.predo 
                      + out$mean$alpha[2] * 1 
                      + out$mean$alpha[3] * 1
                      #+ out$mean$alpha[4] * vveg.predo * 1
                      #+ out$mean$alpha[5] * vveg.predo * 1
                      ),
      col = "darkgreen", lwd = 5)
mtext("Volumen vegetal",       font=1, side=1, cex=1, at=1, line = 1, padj= 1.3, )

# CONVENCIONAL - VERANO
polygon(x= c(vveg.predc, rev(vveg.predc)), y= c(cri4[1,], rev(cri4[2,])), 
        col =  adjustcolor("magenta", alpha.f = 0.15), border = NA)
#lines(vveg.predc, cri4[2,], lty = 2, lwd = 2, col = "red")
lines(vveg.predc, exp(out$mean$mu.alpha 
                      + out$mean$alpha[1] * vveg.predc 
                      + out$mean$alpha[2] * 0 
                      + out$mean$alpha[3] * 1
                      #+ out$mean$alpha[4] * vveg.predc * 0
                      #+ out$mean$alpha[5] * vveg.predc * 1
                      ),
      col = "darkmagenta", lwd = 5)
legend(0.8,250,legend=c( "Convencional","Orgánico" ), 
       lty=c(1,1),cex=1, col=c("darkmagenta","darkgreen"), bty="n", lwd=1.5,
       text.font=1)
#####
dev.off()

