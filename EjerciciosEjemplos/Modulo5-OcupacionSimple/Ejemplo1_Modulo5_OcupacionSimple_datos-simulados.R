############################################################
##### CURSO METODOS CUANTITATIVOS DETECCION IMPERFECTA #####
###########        UNRC - JUNIO JULIO 2016       ###########
###########             A. P. GOIJMAN            ###########
############################################################
###########       Ejemplo Ocupacion en JAGS      ###########
###########     Una Especie, Estacion Simple     ###########
############################################################
########              Ejemplos de:                  ########
########   Kéry, M., & Schaub, M. (2012). Bayesian  ########
########    population analysis using WinBUGS: a    ######## 
########          hierchical perspective.           ######## 
########              Academic Press.               ########
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

######################################################
########     Una Especie, Estacion Simple     ########
########             Datos simulados          ########
######################################################

########################################
### Simulo los datos 
########################################
# Muestreo
R <- 200    # número de sitios
T <- 3      # numero de repeticiones

# Determinar parametros
psi <- 0.8    # Probabilidad de Ocupancion
p <- 0.5      # Probabilidad de Deteccion

# Creo una estructura para contener los conteos
y <- matrix(NA, nrow = R, ncol = T)

# Proceso ecologico: Muestreo la ocurrencia verdadera (z, si/no) de una 
# Bernoulli (probabilidad de ocurrencia = psi)
z <- rbinom(n = R, size = 1, prob = psi)  # Estado Latente de ocurrencia

# Proceso de Observacion: Muestra de observaciones de detection/nodeteccion de una 
# Bernoulli(con p) si z=1
for (j in 1:T){
   y[,j] <- rbinom(n = R, size = 1, prob = z * p)
   }

# Datos generados
y

# Vemos la verdad y nuestras observaciones imperfectas
sum(z)                 # Ocupacion Real sobre 200 sitios muestreados
sum(apply(y, 1, max))  # Ocupacion Observada

########################################
### Analisis de los datos 
########################################
# especificar modelo en JAGS
sink("ocup-simp-ej.jags")
cat("
model {

# Priors
psi ~ dunif(0, 1)
p ~ dunif(0, 1)

# Likelihood
# Modelo Ecologico para ocurrencia real
for (i in 1:R) {
   z[i] ~ dbern(psi)
   p.eff[i] <- z[i] * p

   # Modelo de Observacion para observaciones replicadas deteccion/no deteccion
   for (j in 1:T) {
      y[i,j] ~ dbern(p.eff[i])    # y condicional a z=1 (observac condicional a la presencia)
      } #j
   } #i

# Cantidades Derivadas
occ.fs <- sum(z[])       # Numero de sitios ocupados sobre 200
}
",fill = TRUE)
sink()

# Junto los datos
win.data <- list(y = y, R = nrow(y), T = ncol(y))

# Valores Iniciales
zst <- apply(y, 1, max)		# Ocurrencia Observada como valores iniciales de z
inits <- function() list(z = zst)

# Parametros monitoreados
params <- c("psi", "p", "occ.fs")

# MCMC settings
ni <- 15000
nt <- 2
nb <- 5000
nc <- 3

# Llamo JAGS desde R (BRT < 1 min)
out <- jags(win.data, inits, params, "ocup-simp-ej.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

########################################
### Miro los Resultados 
########################################
# Summarize posteriors
print(out, dig = 2)

# Veo las cadenas MCMC
traceplot(out, parameters=NULL)

# Veo las cadenas y la densidad ("probabilidad") de los parametros
plot(out)

# Valores de los parametros y Intervalos de credibilidad
whiskerplot(out, c('psi', 'p'), quantiles=c(0.025,0.975), zeroline=TRUE)

######################################################
########     Una Especie, Estacion Simple     ########
########            con covariables           ########
########             Datos simulados          ########
######################################################

########################################
### Simulo los datos 
########################################
# Definir funcion para generar datos de distribucion de la especie
data.fn <- function(R = 200, T = 3, xmin = -3, xmax = 1, alpha.psi = -1, 
beta.psi = 3, alpha.p = 1, beta.p = 6) {

   y <- array(dim = c(R, T))	# Arreglo para los conteos

   # Proceso Ecologico
   # Valores de la Covariable
   X <- sort(runif(n = R, min = xmin, max = xmax))

   # Relacion entre ocurrencia esperada y covarible
   psi <- plogis(alpha.psi + beta.psi * X)	# Aplicar logit inverso

   # Agrego Ruido Bernoulli: sorteo occurrencia del indicador z desde Bernoulli(psi)
   z <- rbinom(n = R, size = 1, prob = psi)
   occ.fs <- sum(z)	# Ocupacion de una muestra Finita (see Royle and Kéry 2007)

   # Proceso de observacion
   # Relacion entre prob deteccion - covariable
   p <- plogis(alpha.p + beta.p * X)

   # Hacer 'censo'
   p.eff <- z * p
   for (i in 1:T){
      y[,i] <- rbinom(n = R, size = 1, prob = p.eff)
      }

   # Regresion Sencilla (p=1)
   naive.pred <- plogis(predict(glm(apply(y, 1, max) ~ X + I(X^2), family = binomial)))

   # Graficar caracteristicas del sistema simulado
   par(mfrow = c(2, 2))
   plot(X, psi, main = "Ocurrencia Esperada", xlab = "Covariable", ylab = "Probabilidad de Ocupacion", las = 1, type = "l", col = "red", lwd = 3, frame.plot = FALSE)
   plot(X, z, main = "Ocurrencia Real", xlab = "Covariable", ylab = "Occurrncia", las = 1, frame.plot = FALSE, col = "red",)
   plot(X, p, ylim = c(0,1), main = "Probabilidad de Deteccion", xlab = "Covariable", ylab = "p", type = "l", lwd = 3, col = "red", las = 1, frame.plot = FALSE)
   plot(X, naive.pred, main = "Observaciones Deteccion/nodeteccion \n y analsis convencional p=1", xlab = "Covariable", ylab = "Ocupacion Aparente", ylim = c(min(y), max(y)), type = "l", lwd = 3, lty = 2, col = "blue", las = 1, frame.plot = FALSE)
   points(rep(X, T), y)

   # Mostrar (Return) cosas
   return(list(R = R, T = T, X = X, alpha.psi = alpha.psi, beta.psi = beta.psi, alpha.p = alpha.p , beta.p = beta.p, psi = psi, z = z, occ.fs = occ.fs, p = p, y = y))
   }

sodata <- data.fn()
str(sodata)                 # Look at data

# Analisis con p=1
summary(glm(apply(y, 1, max) ~ X + I(X^2), family = binomial, data = sodata))

## Resumen de datos generados
head(sodata$y)   # historia de detecciones
head(sodata$X)   # covariable
nrow(sodata$y)   # Sitios
ncol(sodata$y)  # Repeticiones

########################################
### Analisis de los datos 
########################################

# Especificar modelo en lenguaje JAGS
sink("model.jags")
cat("
model {

# Priors
alpha.occ ~ dunif(-10, 10)
beta.occ ~ dunif(-10, 10)
alpha.p ~ dunif(-10, 10)
beta.p ~ dunif(-10, 10)

# Verosimilitud
for (i in 1:R) {
   # Modelo del Estado Verdadero para el estado verdadero parcialmente observado
   z[i] ~ dbern(psi[i])             # Ocupacion Real z en el sitio i
   logit(psi[i]) <- alpha.occ + beta.occ * X[i]

   for (j in 1:T) {
      # Modelo de Observacion para las observaciones actuales
      y[i,j] ~ dbern(p.eff[i,j])    # Deteccion-no deteccion en i  j
      p.eff[i,j] <- z[i] * p[i,j]
      logit(p[i,j]) <- alpha.p + beta.p * X[i]
      } #j
   } #i

# Cantidades Derivadas
occ.fs <- sum(z[])       # Numero de sitios ocupados entre los estudiados
}
",fill = TRUE)
sink()

# Juntar datos
win.data <- list(y = sodata$y, X = sodata$X, R = nrow(sodata$y), T = ncol(sodata$y))

# Valores iniciales
zst <- apply(sodata$y, 1, max)   #Valores iniciales buenos son esenciales
inits <- function(){list(z = zst, alpha.occ = runif(1, -3, 3), beta.occ = runif(1, -3, 3), alpha.p = runif(1, -3, 3), beta.p = runif(1, -3, 3))}

# Parameteros monitoreados
params <- c("alpha.occ", "beta.occ", "alpha.p", "beta.p", "occ.fs")

# Seteo de MCMC 
ni <- 10000
nt <- 8
nb <- 5000
nc <- 3

# Llamar JAGS desde R (BRT 1 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

TRUTH <- c(sodata$alpha.psi, sodata$beta.psi, sodata$alpha.p, sodata$beta.p, sum(sodata$z))
print(cbind(TRUTH, out$summary[c(1:4, 6), c(1,2,3,7)]), dig = 3)
sum(apply(sodata$y, 1, sum) > 0)# Numero aparente de sitios ocupados

naive.pred <- plogis(predict(glm(apply(sodata$y, 1, max) ~ X + I(X^2), family = binomial, data = sodata)))
lin.pred2 <- out$mean$alpha.occ + out$mean$beta.occ * sodata$X

plot(sodata$X, sodata$psi, ylim = c(0, 1), main = "Real (rojo), p<1 (negro), \n p=1 (azul)", ylab = "Probabilidad de Ocupacion", xlab = "Covariable", type = "l", lwd = 3, col = "red", las = 1, frame.plot = FALSE)
lines(sodata$X, naive.pred, ylim = c(0 ,1), type = "l", lty = 2, lwd = 3, col = "blue")
lines(sodata$X, plogis(lin.pred2), ylim = c(0, 1), type = "l", lty = 1, lwd = 2, col = "black")
 

######################################################
########     Una Especie, Estacion Simple     ########
########            con covariables           ########
########              Datos Reales            ########
######################################################

# Una especie amenazada, que pone huevos en la madera de arboles muertos. Una buena estrategia
# de busqueda de esta especie rara es buscar pilas de madera muerta.
# Los datos consisten en conteos replicados en 27 sitios en 2009, con hasta 6 replicas.
# Las covariables son si estaba cerca del borde, las fechas y hora del dia.

data <- read.table("bluebug.txt", header = TRUE)
data

# Colectar datos en las estructuras correspondientes
y <- as.matrix(data[,4:9])         # armar matriz "as.matrix" esencial para JAGS
y[y>1] <- 1                        # Reducir conteos a 0/1

edge <- data$forest_edge           #covariable de borde
dates <- as.matrix(data[,10:15])   #covariable de fechas
hours <- as.matrix(data[,16:21])   #covariable de horas

# Estandarizar covariables
mean.date <- mean(dates, na.rm = TRUE)
sd.date <- sd(dates[!is.na(dates)])
DATES <- (dates-mean.date)/sd.date     # Estandarizar fecha
DATES[is.na(DATES)] <- 0               # Agregar ceros a datos NA (media)

mean.hour <- mean(hours, na.rm = TRUE)
sd.hour <- sd(hours[!is.na(hours)])
HOURS <- (hours-mean.hour)/sd.hour      # Estandarizar hora
HOURS[is.na(HOURS)] <- 0                # Agregar ceros a datos NA (media)

########################################
### Analisis de los datos 
########################################
# Especificar modelo en JAGS
sink("model.jags")
cat("
model {

# Priors
alpha.psi ~ dnorm(0, 0.01)
beta.psi ~ dnorm(0, 0.01)
alpha.p ~ dnorm(0, 0.01)
beta1.p ~ dnorm(0, 0.01)
beta2.p ~ dnorm(0, 0.01)
beta3.p ~ dnorm(0, 0.01)
beta4.p ~ dnorm(0, 0.01)

# Likelihood
# Modelo Ecologico para el estado real parcialmente observado
for (i in 1:R) {
   z[i] ~ dbern(psi[i])                    # Ocurrencia Real z en sitio i
   psi[i] <- 1 / (1 + exp(-lpsi.lim[i]))
   lpsi.lim[i] <- min(999, max(-999, lpsi[i]))   #restringe los valores mas extremos en escala logit
   lpsi[i] <- alpha.psi + beta.psi * edge[i]

   # Modelo de Observacion
   for (j in 1:T) {
      y[i,j] ~ dbern(mu.p[i,j])	# Deteccion-no detection en i j
      mu.p[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(999, max(-999, lp[i,j]))
      lp[i,j] <- alpha.p + beta1.p * DATES[i,j] + beta2.p * pow(DATES[i,j], 2) + beta3.p * HOURS[i,j] + beta4.p * pow(HOURS[i,j], 2)
      } #j
   } #i

# Derived quantities
occ.fs <- sum(z[])                             # Numero de sitios ocupados
mean.p <- exp(alpha.p) / (1 + exp(alpha.p))    # Deteccion media
}
",fill = TRUE)
sink()

# Juntar datos
win.data <- list(y = y, R = nrow(y), T = ncol(y), edge = edge, DATES = DATES, HOURS = HOURS)

# Valores iniciales
zst <- apply(y, 1, max, na.rm = TRUE)	# Good starting values crucial
inits <- function(){list(z = zst, alpha.psi=runif(1, -3, 3), alpha.p = runif(1, -3, 3))}

# Parametros monitoreados
params <- c("alpha.psi", "beta.psi", "mean.p", "occ.fs", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta4.p")

# seteo MCMC
ni <- 30000
nt <- 10
nb <- 20000
nc <- 3

# Llama JAGS desde R (BRT < 1 min)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# Resumen de posterior
print(out, dig = 2)
traceplot (out)
# Valores de los parametros y Intervalos de credibilidad
whiskerplot(out, c('alpha.psi', 'mean.p'), quantiles=c(0.025,0.975), zeroline=TRUE)
whiskerplot(out, c('beta.psi'), quantiles=c(0.025,0.975), zeroline=TRUE)
whiskerplot(out, c('alpha.p', 'mean.p'), quantiles=c(0.025,0.975), zeroline=TRUE)

# Distribucion posterior del numero de sitios ocupados en la muestra actual
hist(out$sims.list$occ.fs, nclass = 30, col = "gray", main = "", xlab = "Numero de pilas de madera ocupadas (occ.fs)", xlim = c(9, 27))
abline(v = 10, lwd = 2) # Numero Observado
  
par(mfrow = c(2, 1))
hist(plogis(out$sims.list$alpha.psi), nclass = 40, col = "gray", main = "Interior del Bosque", xlab = "Probabilidad de Ocupacion", xlim = c(0, 1))
hist(plogis(out$sims.list$alpha.psi+ out$sims.list$beta.psi), nclass = 40, col = "gray", 
main = "Borde del Bosque", xlab = "Probabilidad de Ocupacion", xlim = c(0, 1))
 

# Predict effect of time of day with uncertainty
mcmc.sample <- out$mcmc.info$n.samples   #n de simulaciones

# Secuencia de la covariable
original.date.pred <- seq(0, 60, length.out = 30) 
original.hour.pred <- seq(180, 540, length.out = 30)

# Estandarizo la covariable p las estimaciones
date.pred <- (original.date.pred - mean.date)/sd.date
hour.pred <- (original.hour.pred - mean.hour)/sd.hour

# predicciones de deteccion en relacion a la covariable
p.pred.date <- plogis(out$mean$alpha.p + out$mean$beta1.p * date.pred + out$mean$beta2.p * date.pred^2 )
p.pred.hour <- plogis(out$mean$alpha.p + out$mean$beta3.p * hour.pred + out$mean$beta4.p * hour.pred^2 )

#arreglos para poner los datos
array.p.pred.hour <- array.p.pred.date <- array(NA, dim = c(length(hour.pred), mcmc.sample))
for (i in 1:mcmc.sample){
   array.p.pred.date[,i] <- plogis(out$sims.list$alpha.p[i] + out$sims.list$beta1.p[i] * date.pred + out$sims.list$beta2.p[i] * date.pred^2)
   array.p.pred.hour[,i] <- plogis(out$sims.list$alpha.p[i] + out$sims.list$beta3.p[i] * hour.pred + out$sims.list$beta4.p[i] * hour.pred^2)
   }

# Grafico algunas muestras a azar del MCMC
sub.set <- sort(sample(1:mcmc.sample, size = 20))

par(mfrow = c(2, 1))
plot(original.date.pred, p.pred.date, main = "", ylab = "Probabilidad de Deteccion", 
xlab = "Fecha (1 = 1 Julio)", ylim = c(0, 1), type = "l", lwd = 3, frame.plot = FALSE)
for (i in sub.set){
   lines(original.date.pred, array.p.pred.date[,i], type = "l", lwd = 1, col = "gray")
   }
lines(original.date.pred, p.pred.date, type = "l", lwd = 3, col = "blue")

plot(original.hour.pred, p.pred.hour, main = "", ylab = "Probabilidad de Deteccion",
 xlab = "Hora del dia (mins luego del mediodia)", ylim = c(0, 1), type = "l", lwd = 3, frame.plot = FALSE)
for (i in sub.set){
   lines(original.hour.pred, array.p.pred.hour[,i], type = "l", lwd = 1, col = "gray")
   }
lines(original.hour.pred, p.pred.hour, type = "l", lwd = 3, col = "blue")

