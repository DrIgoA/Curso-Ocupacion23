#comn los dos tiempos en el array de datos
#variables escaladas 
#escalas 150 200 300 400 500 600  - con % de tierra cultivada
#SIN ambientes, solo el %tC

############################################################
rm(list=ls(all=TRUE))

setwd("/home/vanesa/Documentos/Vanesa UNRC/Publicaciones/Paper complejidad/analisis/Analisis revision landscape 2")
setwd("C:\\Users\\andrea\\Downloads")

library(jagsUI)    #paquete JAGS
library(boot)
library(matrixStats)

data<-read.csv("datos7.csv",header = T)
str(data)
head(data)

attach(data)
str(roed.data <- list(y=y,J=J, IS150=IS150, sp =sp, season = season))

  ###Modelos solo PA
### JAGS code
sink("modelo-prueba.jags")  
cat("
    model{
    
    # Hyperpriors 
    lpsi ~ dnorm(0, 0.1)
    b1 ~   dnorm(0, 0.1)
    p ~   dunif(0, 1)
    
  # Likelihood
  for (i in 1:50){
    # Ecological model, process model (true occurrence at site i) 
      logit(psi[i]) <- lpsi[sp[i]] + b1[sp[i]] *IS150[i]
      z[i] ~ dbin(psi[i],1)
    
    ## Modelo para la deteccion
        #logit(p) <- lp[sp[i]]     
        tmp[i] <- p * z[i]
        y[i] ~ dbin(tmp[i], J[i])
  }
    }
    
    
    
      ",fill=TRUE)
# 
# model { 
#   for(i in 1:M){ 
#     z[i] ~ dbin(psi[i],1)      
#     logit(psi[i]) <- b0 + b1 * alti[i] + b2 * alti2[i] + b3 * bosq[i] 
#     tmp[i] <- z[i] * p             
#     y[i] ~ dbin(tmp[i], J[i])  
#   } 
#   
  
  sink()

# valores iniciales
#inits <- function()     list(z=zst  )

inits <- function() list(lpsi=runif(1), b1=runif(1))
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
out1= jags(roed.data,inits, params1, "modelo-prueba.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

