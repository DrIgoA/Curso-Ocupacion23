
library(jagsUI)    #paquete JAGS

#Ejemplo tamaño medio de zorzales (n=10)
# media de una distribucion normal

data <- list(size = c(7.9,8.1,11,10.6,9.2,8,9.8,10.1,10.9,9))

sink("zorzales.jags")
cat("

model
{
# previas no informativas
mean ~ dnorm (0, 1.0E-6)       # tamaño medio de los zorzales
varianza ~ dlnorm (0 ,1.0E-6)  # varianza tamaño zorzales

prec <- 1/varianza             # pasar de varianza a precision 

for (i in 1:10)                # para cada uno de los Zorzales
{
size[i] ~ dnorm (mean, prec)   # tamaño zorzal trazado de una distribución normal
}
}

",fill = TRUE)
sink()

inits <- function() list(varianza=100, mean = 10)

# MCMC settings
ni <- 10000  # número de iteraciones
nt <- 2      # tasa de thining
nb <- 1000   # iteraciones para el burn in
nc <- 3      # número de cadenas que corremos

parameters <- c("mean", "varianza")

out <- jags(data, inits, parameters, "zorzales.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

out

str(out)

plot(out)

traceplot(out, parameters=c("mean", "varianza")) 

densityplot(out, parameters=c("mean", "varianza")) 

