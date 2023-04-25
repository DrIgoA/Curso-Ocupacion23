##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########           Ejemplo Bayes              #########
###########                                      #########
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
    mean ~ dnorm (0, 1.0E-6)              # tamaño medio de los zorzales
    varianza ~ dlnorm (0 ,1.0E-6)         # varianza tamaño zorzales

  #likelihood
  for (i in 1:10) {                       # para cada uno de los Zorzales
     size[i] ~ dnorm (mean, 1/varianza)   # tama?o zorzal trazado de una distribuci?n normal  
 
   }
}

",fill = TRUE)
sink()

# Indicamos los valores iniciales 
inits <- function() list(varianza=100, mean = 10)

# MCMC settings
ni <- 100000     # numero de iteraciones
nt <- 2         # tasa de thining
nb <- 10000      # iteraciones para el burn in
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

## Graficos
## Grafico Media
# 1) Previa
     x <- seq(-5000,5000, 0.5)  #graficar solo la distribucion previa
     prior_Non <- dnorm(x, mean = 0, sd = sqrt(1/1.0E-6))  # sd=1000
     plot(x, prior_Non, main = "Previa No informativa",
     type = 'l', bty ="l", col= "darkgreen", xlab = "", xlim = c(-5000,5000))

     #Si miramos la previa en torno a los valores del tamano de los zorzales
     plot(x, prior_Non, main = "Previa No informativa",
          type = 'l', bty ="l", col= "darkgreen", xlab = "", xlim = c(0,20))

# 2) Posterior con previa no informativa
x <- seq(0,20, 0.1)
plot(density(out$sims.list$mean), main = "Posterior con previa No informativa",
     type = 'l', col= "darkgreen", xlab = "")

# 2) Graficar conjuntamente Posterior y previa no informativa
x <- seq(0,20, 0.1)  #graficar la distribucion previa en rango de los valores
prior_Non <- dnorm(x, mean = 0, sd = sqrt(1/1.0E-6))  # sd=1000

#grafica la distribucion posterior
plot(density(out$sims.list$mean), main = "Posterior con previa No informativa",
     type = 'l', col= "darkgreen", xlab = "")
lines(x, prior_Non,  col = 'red')                  # Agrega al grafico la previa del modelo
legend('topright', lty=1, c('Previa No informativa', "Posterior"), col = c("red", "darkgreen"),)



## Actividad ----
# 1) Como harias para agregar previas informativas a este ejercicio?
sink("zorzales-informativa.jags")
cat("

model
{
  # Previas INFORMATIVAS
    mean ~ dnorm (9.45, 1/(0.8)*(0.8))       # tamaño medio de los zorzales
    varianza ~ dlnorm (0.81, 1/(0.52*0.52))  # varianza tamaño zorzales

  # Likelihood
  for (i in 1:10) {                          # para cada uno de los Zorzales
     size[i] ~ dnorm (mean, 1/varianza)      # tama?o zorzal trazado de una distribuci?n normal  
 
   }
}

",fill = TRUE)
sink()

out_info <- jags(data, inits, parameters, "zorzales-informativa.jags",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

out_info

## Grafico Media
# 1) Previa Informativa
x <- seq(-5000,5000, 0.1)  #graficar solo la distribucion previa
prior_Inf <- dnorm(x, mean = 9.45, sd = 0.8)  # Va desvio 
plot(x, prior_Inf, main = "Previa informativa",
     type = 'l', bty ="l", col= "darkgreen", xlab = "", xlim = c(0,20))

# 2) Posterior con previa informativa
x <- seq(0,20, 0.1)
plot(density(out_info$sims.list$mean), main = "Posterior con previa informativa",
     type = 'l', col= "darkgreen", xlab = "")

# 2) Graficar conjuntamente Posterior y previa informativa
x <- seq(0,20, 0.1)  #graficar la distribucion previa en rango de los valores
prior_Inf <- dnorm(x, mean = 9.45, sd = 0.8)  # Va desvio 

#grafica la distribucion posterior
plot(density(out_info$sims.list$mean), main = "Posterior con previa informativa",
     type = 'l', col= "darkgreen", xlab = "")
lines(x, prior_Inf,  col = 'red')                  # Agrega al grafico la previa del modelo
legend('topright', lty=1, c('Previa informativa', "Posterior"), col = c("red", "darkgreen"),)


# Sintesis de graficos
tiff(filename = "previa-posterior_sintesis",
     res = 300, width = 180,height = 250, units = "mm")
     

layout(matrix(1:2, nrow = 2, ncol= 1))

plot(density(out$sims.list$mean), main = "Posterior con previa No informativa",
     type = 'l', col= "darkgreen", xlab = "", xlim = c(7,12))
lines(x, prior_Non,  col = 'red')                  # Agrega al grafico la previa del modelo
legend('topright', lty=1, c('Previa No informativa', "Posterior"), col = c("red", "darkgreen"),)

plot(density(out_info$sims.list$mean), main = "Posterior con previa informativa",
     type = 'l', col= "darkgreen", xlab = "", xlim = c(7,12))
lines(x, prior_Inf,  col = 'blue')                  # Agrega al grafico la previa del modelo
legend('topright', lty=1, c('Previa informativa', "Posterior"), col = c("blue", "darkgreen"),)

dev.off()

