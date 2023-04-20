##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########      Ejercicio  jerarquico           #########
##########################################################
########              Ejemplo de:                   ######
########          Introduction to WinBUGS           ######
########             for ecologists                 ######
########    A bayesian approach to regression,      ######
########  ANOVA, mixed models an related analyses   ######
########            Marc Kery - 2010                ######
##########################################################

#---------------------------------------------------------------------------
# Estimación de media con distribucion Poisson
# ---------------------------------------------------------------------------

# Investigadores en NY midieron la densidad de arboles en 10 cuadrantes (400m2 c/U)
# en Van Cortlandt Park

# Una distribución previa no informativa para la densidad media del red oak en cda cuadrante
# debe tomar valores positivos y tener una amplia distribucion.

#Una distribución lonormal que tenga media 0 y desvio standar 1000 para valores 
# log-transform puede funcionar.


# La presición que es la que se usa, es de 0.04 (1/varianza)

library(jagsUI)    #paquete JAGS

# Agrupar los datos en una lista para que pueda ser leida por JAGS
data <- list(y = c(6,0,1,2,1,7,1,5,2,0))

# Modelo
sink("media.jags")
cat("
model
{
  for (i in 1:10)
  {
    y[i] ~ dpois (m)
    }
  
  m ~ dlnorm(0.0, 1.0E-6)
}

",fill = TRUE)
sink()

#MCMC settings
ni <- 100000
nt <- 2
nb <- 10000
nc <- 3

#Parametros
parameters <- c("m")

inits <- function() list(m = 5)

#Run JAGS model
# out <- jags(data, 
#             inits=inits, 
#             parameters, 
#             "media-poisson.jags", 
#             n.chains = nc, 
#             n.thin = nt, 
#             n.iter = ni,
#             n.burnin = nb)
# 
# 


## Estimación de media con Poisson~~~~~~~~~~~~~~####
### Con variacion entre la densidad promedio de los caudarte###


# Investigadores en NY midieron la densidad de arboles en 10 cuadrantes (400m2 c/U)
# en Van Cortlandt Park

# Una distribución previa no informativa para la densidad media del red oak en cda cuadrante
# debe tomar valores positivos y tener una amplia distribucion.

#Una distribución lonormal que tenga media 0 y desvio standar 1000 para valores 
# log-transform puede funcionar.


# La presición que es la que se usa, es de 0.04 (1/varianza)
library(jagsUI)    #paquete JAGS

# Agrupar los datos en una lista para que pueda ser leida por JAGS
data <- list(y = c(6,0,1,2,1,7,1,5,2,0))

# Modelo
sink("media-poisson-variacion.jags")
cat("
model
{
  for (i in 1:10)        # for each of the 10 quadrants
  {
    mean[i] ~ dlnorm(m, tau)   # mean density drawn from lognormal
    y[i] ~ dpois (mean[i])     # no. of plants drawn from Poisson
    }
  
  m ~ dlnorm(0, 1.0E-6)       # mean of the log density of plants
  sd ~ dunif(0,10)            # sd of the log density of plants
  
  tau <- 1/(sd*sd)
  
}

",fill = TRUE)
sink()

#MCMC settings
ni <- 100000
nt <- 2
nb <- 10000
nc <- 3

#Parametros
parameters <- c("mean")

inits <- function() list(m = 5)

#Run JAGS model
# out <- jags(data,
#             inits=inits,
#             parameters,
#             "media-poisson-variacion.jags",
#             n.chains = nc,
#             n.thin = nt,
#             n.iter = ni,
#             n.burnin = nb)
# out
# 
# 


