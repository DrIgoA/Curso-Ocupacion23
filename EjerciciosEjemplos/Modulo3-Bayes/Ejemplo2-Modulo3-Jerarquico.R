##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########           Ejemplo Jerarquico         #########
##########################################################
########              Ejemplo de:                   ######
########          Introduction to WinBUGS           ######
########             for ecologists                 ######
########    A bayesian approach to regression,      ######
########  ANOVA, mixed models an related analyses   ######
########            Marc Kery - 2010                ######
##########################################################

# Ejemplo media 

# Conteos del numero de ranas en 15 sitios

library(jagsUI)    #paquete JAGS

# Agrupar los datos en una lista para que pueda ser leida por JAGS
data <- list(y = c(6,3,1,2,1,7,1,5,2,8,2,3,5,9,11))

# Modelo
sink("pois-ranas.jags")
cat("
model
    {
     # Previa
      lambda ~ dnorm(0, 1.0E-6)
    
     #likelihood
      for(i in 1:15) {
          y[i]~ dpois(lambda)     
    }
}

",fill = TRUE)
sink()

#MCMC settings
ni <- 10000
nt <- 2
nb <- 1000
nc <- 3

#Parametros
parameters <- c("lambda")

# Initial values
inits <- function() list(lambda = 5)

#Run JAGS model
out <- jags(data, inits, parameters, "pois-ranas.jags", 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


# ---------------------------------------------------------------------------
# Jerarquico
#---------------------------------------------------------------------------

# Modelo
sink("pois-ranas2.jags")
cat("
model
    {
     # Previas
       mu~dnorm(0,1.0E-6)
       sd~dunif(0,10)
       tau <- 1/(sd*sd)
    
     #Likelihood  
      for(i in 1:15) {
        lambda[i] ~ dlnorm(mu, tau)
        y[i]~ dpois(lambda[i])     
    }
    
      median <- exp(mu)  # Se usa esta linea porque es la inversa del log normal
}

",fill = TRUE)
sink()

#MCMC settings
ni <- 100000
nt <- 2
nb <- 10000
nc <- 3

#Parametros
parameters <- c("lambda", "mu", "median")

# Initial values
inits <- function() list(lambda = c(5,2,1,3,4,2,5,2,1,3,4,2,5,2,1), mu=2,sd=5)

#Run JAGS model
 out <- jags(data, 
             inits,              parameters, 
             "pois-ranas2.jags", 
             n.chains = nc, 
             n.thin = nt, 
             n.iter = ni,
             n.burnin = nb)

#--------------------------------------------------------------------------- 
# Actividad
# 1. Que representran los 15 valores de lambda estimados en el segundo modelo?
# 2. Que representa mu en el segundo modelo?
# 3. Como modificaria la probabilidad previa para que sea informativa?
# --------------------------------------------------------------------------
 