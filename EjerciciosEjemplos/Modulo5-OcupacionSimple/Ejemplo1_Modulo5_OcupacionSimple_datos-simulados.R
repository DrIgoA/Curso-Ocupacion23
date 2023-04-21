############################################################
##### CURSO Modelado y estimación de ocupación para  #######
#####  poblaciones y comunidades de especies bajo    #######
#####           enfoque Bayesiano.                   #######
#######      CCT Mendoza - ABRIL 2023                #######
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
setwd("")

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
