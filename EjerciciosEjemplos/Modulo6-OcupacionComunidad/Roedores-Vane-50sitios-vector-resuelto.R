##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########         Ejercicio Ocupacion          #########
###########          Roedores - Vectores         #########
##########################################################


# Los datos corresponden a dos muestreos (season) de roedores realizados en la provincia 
# de Cordoba. El diseño posee 50 sitios, en los cuales se coloco 1 linea con 
# 20 trampas de captura viva tipo Shermann.  
# Se capturaron 9 especies de roedores: Akodon azarae, A. dolores, Calomys 
# laucha, C. musculinus, C. venustus, Oligorizomys flavesces, Oxymicterus rufus, 
# Tylamys pallidor, Monodelphys dimidiata.
# 

rm(list=ls(all=TRUE))

library(jagsUI)    #paquete JAGS

data<-read.csv("datos_roedores.csv",header = T)
str(data)
head(data, 5)

data1 <- data[1:50,]

attach(data)
str(roed.data <- list(y=y, J=J,is = is, sp =sp, season = season))

### JAGS code
sink("modelo.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 1)
    mu.b1 ~   dnorm(0, 1)
    mu.lp ~   dnorm(0, 1)

      #  p ~   dunif(0,1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)

    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)

  for (k in 1:9) {
    lpsi[k]~ dnorm(mu.lpsi, tau.lpsi)
    b1[k]~ dnorm(mu.b1, tau.b1)
    lp[k]~ dnorm(mu.lp, tau.lp)
  }

  # Likelihood
  for (i in 1:900){
      # Ecological model, process model (true occurrence at site i) 
      logit(psi[i]) <- lpsi[sp[i]] + b1[sp[i]] * is[i]
    
        z[i] ~ dbin(psi[i],1)
    
    ## Modelo para la deteccion
        logit(p[i]) <- lp[sp[i]]     
        tmp[i] <- p[i] * z[i]
        y[i] ~ dbin(tmp[i], J[i])
  }
    }
    
    
      ",fill=TRUE)

  sink()

# valores iniciales

inits <- function() list(lpsi=runif(9), b1=runif(9),
                         mu.lpsi=runif(1),mu.b1=runif(1))

params1 <- c("lpsi", "b1", "lp","mu.lpsi","mu.lp","mu.b1")
params2 <- c("Nsite", "Nocc.fs")

# ajustes de MCMC
ni <- 60000
nt <- 10
nb <- 10000
nc <- 3

# llamar JAGS desde R 
out1= jags(roed.data,inits, params1, "modelo.jags", n.chains=nc, 
               n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace




library(denstrip)
plot(out1$sims.list$b1, xlim=c(-3, 3), ylim=c(1, 9), xlab="", ylab="", type="n", axes =
       F, main = "Density strip plots")
axis(1)
axis(2, at = 1:9, labels = c('Aa','Ad','Cm','Cl','Cv','Of', 'Or','Md', 'Tp'), las = 1)
abline(v = c(--1,-0.5,0.5,1), col = "grey") ; abline(v = 0)
for(k in 1:9){
  denstrip(unlist(out1$sims.list$b1[,k]), at = k,
           colmax = "#4292c6", colmin = "#f7fbff")
}

















sink("modelo-prueba.jags")  
cat("
    model{
    
    # Hyperpriors 
    mu.lpsi ~ dnorm(0, 0.1)
    mu.b1 ~   dnorm(0, 0.1)
    p ~   dunif(0,1)
    
    tau.lpsi <- pow(sd.lpsi,-2)   
    tau.lp <- pow(sd.lp,-2)   
    tau.b1 <- pow(sd.b1,-2)

    sd.lpsi ~ dunif(0,1)
    sd.lp ~ dunif(0,1)
    sd.b1 ~ dunif(0,1)

  for (k in 1:9) {
    lpsi[k]~ dnorm(mu.lpsi, tau.lpsi)
    b1[k]~ dnorm(mu.b1, tau.b1)
  #  lp[k]~ dnorm(mu.lp, tau.lp)
  }

  # Likelihood
  for (i in 1:900){
    # Ecological model, process model (true occurrence at site i) 
      logit(psi[i]) <- lpsi[sp[i]] + b1[sp[i]] *IS150[i]
            z[i] ~ dbin(psi[i],1)
    
    ## Modelo para la deteccion
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

inits <- function() list(lpsi=runif(9), b1=runif(9),
                         mu.lpsi=runif(1),mu.b1=runif(1))
# ,sd.lpsi=runif(roed.data$nspec,0.1,5),
# sd.lp=runif(roed.data$nspec,0.1,5), sd.b1=runif(roed.data$nspec,0.1,5))


params1 <- c("lpsi", "b1", "lp","mu.lpsi","mu.lp","mu.b1", "psi")
params2 <- c("Nsite", "Nocc.fs")

# ajustes de MCMC
ni <- 5000
nt <- 10
nb <- 2500
nc <- 3

# llamar JAGS desde R 
out1 = jags(roed.data,inits, params1, "modelo-prueba.jags", n.chains=nc, 
           n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace
