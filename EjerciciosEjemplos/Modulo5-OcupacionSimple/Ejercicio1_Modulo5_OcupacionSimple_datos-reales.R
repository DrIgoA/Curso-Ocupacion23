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
beta3.p ~ dnorm(0, 0.01)
beta4.p ~ dnorm(0, 0.01)

# Likelihood
# Modelo Ecologico para el estado real parcialmente observado
for (i in 1:R) {
   z[i] ~ dbern(psi[i])                    # Ocurrencia Real z en sitio i
   psi[i] <- 1 / (1 + exp(-lpsi[i]))
   lpsi[i] <- alpha.psi + beta.psi * edge[i]

   # Modelo de Observacion
   for (j in 1:T) {
      y[i,j] ~ dbern(mu.p[i,j])	# Deteccion-no detection en i j
      mu.p[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp[i,j]))
      lp[i,j] <- alpha.p + beta3.p * HOURS[i,j] + beta4.p * pow(HOURS[i,j], 2)
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
params <- c("alpha.psi", "beta.psi", "mean.p", "occ.fs", "alpha.p", "beta3.p", "beta4.p")

# seteo MCMC
ni <- 30000
nt <- 10
nb <- 20000
nc <- 3

# Llama JAGS desde R (BRT < 1 min)
out.ej <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# Resumen de posterior
print(out.ej, dig = 2)
traceplot (out)
# Valores de los parametros y Intervalos de credibilidad

par(mfrow = c(2, 1))
whiskerplot(out, c('alpha.psi', 'mean.p'), quantiles=c(0.025,0.975), zeroline=TRUE)
whiskerplot(out.ej, c('alpha.psi', 'mean.p'), quantiles=c(0.025,0.975), zeroline=TRUE)
whiskerplot(out, c('beta.psi'), quantiles=c(0.025,0.975), zeroline=TRUE)
whiskerplot(out.ej, c('beta.psi'), quantiles=c(0.025,0.975), zeroline=TRUE)
whiskerplot(out, c('alpha.p', 'mean.p'), quantiles=c(0.025,0.975), zeroline=TRUE)
whiskerplot(out.ej, c('alpha.p', 'mean.p'), quantiles=c(0.025,0.975), zeroline=TRUE)

whiskerplot(out, c('beta1.p','beta2.p','beta3.p','beta4.p'), quantiles=c(0.025,0.975), zeroline=TRUE)
whiskerplot(out.ej, c('beta3.p','beta4.p'), quantiles=c(0.025,0.975), zeroline=TRUE)



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

