## Ranas jerarquico 1~~~~~~~~~~~~~~####

# Diametro promedio de árboles en un parque de bosque de eucalipto

# Incluye previas no informativas, con incertidumbre en la varianza,
# a diferencia del ejemplo anterior donde se asume que es un valor conocido

# La presición que es la que se usa, es de 0.04 (1/varianza)
library(jagsUI)    #paquete JAGS

# Agrupar los datos en una lista para que pueda ser leida por JAGS



data <- list(y = c(6,3,1,2,1,7,1,5,2,8,2,3,5,9,11))

# Modelo
sink("pois-ranas.jags")
cat("
model
    {
      lambda ~ dnorm(0, 1.0E-6)
      
      for(i in 1:15)
    {
    y[i]~ dpois(lambda)     # tree diameter drawn from normal distribution
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
#out <- jags(data, inits, parameters, "pois-ranas.jags", 
         #   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


## Ranas jerarquico 2~~~~~~~~~~~~~~####

# Diametro promedio de árboles en un parque de bosque de eucalipto

# Incluye previas no informativas, con incertidumbre en la varianza,
# a diferencia del ejemplo anterior donde se asume que es un valor conocido

# La presición que es la que se usa, es de 0.04 (1/varianza)
library(jagsUI)    #paquete JAGS

# Agrupar los datos en una lista para que pueda ser leida por JAGS
data <- list(y = c(6,3,1,2,1,7,1,5,2,8,2,3,5,9,11))

# Modelo
sink("pois-ranas2.jags")
cat("
model
    {
      mu~dnorm(0,1.0E-6)
      sd~dunif(0,10)
      tau <- 1/(sd*sd)
      
      for(i in 1:15)
    {
    lambda[i] ~ dlnorm(mu, tau)
    y[i]~ dpois(lambda[i])     # tree diameter drawn from normal distribution
    }
    
    median <- exp(mu)
    
    
    ### Se usa esta linea porque es la inversa del log normal
}

",fill = TRUE)
sink()

# que representran los 15 valores de lambda
# cual es la median(
#   
# agregar un previa informativa


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
# out <- jags(data, 
#             inits, 
#             parameters, 
#             "pois-ranas2.jags", 
#             n.chains = nc, 
#             n.thin = nt, 
#             n.iter = ni,
#             n.burnin = nb)


## Estimación de media con Poisson~~~~~~~~~~~~~~####

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
sink("media-poisson.jags")
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


