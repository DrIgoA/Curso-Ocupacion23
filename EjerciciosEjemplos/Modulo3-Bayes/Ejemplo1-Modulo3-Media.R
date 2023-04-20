##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########           Ejemplo Bayes              #########
##########################################################

library(jagsUI)    #paquete JAGS

# Ejemplo tamaño medio de zorzales (n=10)
# media de una distribucion normal

# es importante que siempre especifiquemos los datos como una lista
data <- list(size = c(7.9,8.1,11,10.6,9.2,8,9.8,10.1,10.9,9))

# Modelo
sink("zorzales.jags")
cat("

model
{
  # previas no informativas
    mean ~ dnorm (0, 1.0E-6)       # tamaño medio de los zorzales
    varianza ~ dlnorm (0 ,1.0E-6)  # varianza tamaño zorzales

    prec <- 1/varianza             # pasar de varianza a precision 

  for (i in 1:10) {                # para cada uno de los Zorzales
     size[i] ~ dnorm (mean, prec)   # tama?o zorzal trazado de una distribuci?n normal  
 
   }
}

",fill = TRUE)
sink()

# Indicamos los valores iniciales 
inits <- function() list(varianza=100, mean = 10)

# MCMC settings
ni <- 10000     # numero de iteraciones
nt <- 2         # tasa de thining
nb <- 1000      # iteraciones para el burn in
nc <- 3         # numero de cadenas que corremos

# Parametros que se van monitorear 
parameters <- c("mean", "varianza")

# Implementamos el modelo con jags
out <- jags(data, inits, parameters, "zorzales.jags",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# Podemos observar la salida del modelo 
out

# Aquí nos muestra la estructura del resultado de modelo
str(out)

plot(out)

# Con la función traceplot (indicando los parametros que queremos visualizar) 
# podemos ver las cadenas, y ver graficamente si convergen
traceplot(out, parameters=c("mean", "varianza")) 

# En este gráfico vemos la distribucion de probabilidad para los parametros
densityplot(out, parameters=c("mean", "varianza")) 

# Grafico Distribucion
# Previa
x <- seq(0,12, 0.5)
prior_Non <- dlnorm(x, mean = 0, sd = 1.0E-6)

# Grafica la distribucion posterior
plot(density(out$sims.list$mean), main = "Previa No informativa",
     bty = 'l', col= "darkgreen", xlab = "")
lines(x, prior_Non,  col = 'red')                  # Agrega al grafico la previa del modelo
legend('topright', lty=1, c('Previa', "Posterior"), col = c("red", "darkgreen"),)

## Actividad ----
# 1. Como harias para agregar previas informativas a este ejercicio?
# 2. Cambian los resultados con lo propuesto anteriormente? 
# 3. Podes graficar los dos modelos con sus distribuciones previas y posteriores?
  