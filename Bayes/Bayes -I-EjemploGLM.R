############################################################
###########          Ejemplo GLM y  GLMM         ###########
############################################################
########              Basado en e:                  ########
########   K?ry, M., & Schaub, M. (2012). Bayesian  ########
########    population analysis using WinBUGS: a    ######## 
########          hierarchical perspective.         ######## 
########              Academic Press.               ########
############################################################

######################################
### Borrar documentacion anterior de R
######################################
rm(list=ls(all=TRUE))

###################################### 
### Elegir directorio de trabajo
######################################
setwd("C:\\Users\\andrea\\Documents\\GitHub\\Curso-Ocupacion23\\Bayes")

#####################################
#Explorar distribuciones
#####################################
plot(density(rbeta(n=10^6, shape1=2, shape2 = 4)))
hist(rbeta(10^6, 2, 4), nclass=100, col="gray")

hist(rpois(10^6, 4), nclass=100, col="gray")
### ?Como se veria un histograma con poisson, y con binomial?
## ayuda... para ver los parametros de las distribuciones
?dpois
?dbinom
?rpois  #genero numeros al azar

## como generarian numeros al azar bajo una distribucion binomial? 

#####################################
# Modelo Lineal
#####################################

y <- c(25, 37, 68, 79, 86, 139, 49, 91, 111)
A <- factor(c(1,1,1,2,2,2,3,3,3))
x<- c(1, 14, 22, 2, 9,20, 2, 13, 22)

plot(x, y, col = c(rep("red",3),rep("blue",3),rep("green",3)),xlim = c(-1,25),ylim = c(0,140))


###################################################
### GLM Poisson JAGS 
### Conteos en el tiempo
###################################################

# Simulacion de datos
n = 40
alpha = 3.5576
beta1 = -0.0912 
beta2 = 0.0091 
beta3 = -0.00014
# n: Numero de anios
# alpha, beta1, beta2, beta3: coeficientes de un polinomio cubico de conteos anuales
# Genera valores de la covariable en el tiempo

year <- 1:n

# Parte del GLM
log.expected.count <- alpha + beta1 * year + beta2 * year^2 + beta3 * year^3
expected.count <- exp(log.expected.count)  #inversa del log para ver conteos en valores reales

# Ruido: genera parte aleatoria del GLM: Ruido Poisson alrededor de los puntos 
C <- rpois(n = n, lambda = expected.count)

# Grafico datos simulados
plot(year, C, type = "b", lwd = 2, col = "black", main = "", las = 1, ylab = "Tamaño poblacional",
     xlab = "Año", cex.lab = 1.2, cex.axis = 1.2)
lines(year, expected.count, type = "l", lwd = 3, col = "red")

plot(year, log.expected.count, type = "b", lwd = 2, col = "black", main = "", las = 1, ylab = "Log Tamaño poblacional",
     xlab = "Año", cex.lab = 1.2, cex.axis = 1.2)

##############################
#####  GLM Poisson      ######
##############################
fm <- glm(C ~ year + I(year^2) + I(year^3), family = poisson, data = data)
summary(fm)


##############################
#####  Analisis en JAGS ######
##############################

### Cargar Paquete 
library(jagsUI)    #paquete JAGS

# Especificar modelo en JAGS
sink("GLM_Poisson.jags")
cat("
model {

# Priors
alpha ~ dunif(-20, 20)
beta1 ~ dunif(-10, 10)
beta2 ~ dunif(-10, 10)
beta3 ~ dunif(-10, 10)

# Verosimilitud: Componentes clave de un GLM en cada linea
for (i in 1:n){                       #desde 1 hasta n observaciones (n=40)
   C[i] ~ dpois(lambda[i])          # 1. Distribucion de parte aleatoria
   log(lambda[i]) <- log.lambda[i]  # 2. Funcion de enlace (Link function)
   log.lambda[i] <- alpha + beta1 * year[i] + beta2 * pow(year[i],2) + 
   beta3 * pow(year[i],3)            # 3. Predictor Lineal
   } #i
}
",fill = TRUE)
sink()

# Unir datos
# Aca se utiliza la variable sin estandarizar, que muchas veces hace que el modelo no funcione, que no corra
win.data <- list(C = C, n = length(C), year = year)

# Unir datos
# con variable estandarizada!
mean.year <- mean(data$year)             # Mean of year covariate
sd.year <- sd(data$year)                 # SD of year covariate
win.data <- list(C = data$C, n = length(data$C), year = (data$year - mean.year) / sd.year)

# Valores iniciales
inits <- function() list(alpha = runif(1, -2, 2), beta1 = runif(1, -3, 3))

# Parametros monitoreados
params <- c("alpha", "beta1", "beta2", "beta3", "lambda")

# MCMC settings
ni <- 10000
nt <- 2
nb <- 1000
nc <- 3

# Call JAGS from R (BRT < 1 min)
out <- jags(data = win.data, inits = inits, parameters.to.save = params, 
            model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, 
            n.iter = ni, n.burnin = nb)
print(out, dig = 3)
plot(out)
traceplot(out)


##### Como se verían las cadenas MCMC sin convergencia?

# graficamos los resultados anteriores 
plot(1:40, data$C, type = "b", lwd = 2, col = "black", main = "", las = 1, ylab = "Population size", xlab = "Year")
R.predictions <- predict(glm(C ~ year + I(year^2) + I(year^3), family = poisson, data = data), type = "response")
lines(1:40, R.predictions, type = "l", lwd = 3, col = "green")

out$mean$lambda
JAGS.predictions <- out$mean$lambda
lines(1:40, JAGS.predictions, type = "l", lwd = 3, col = "blue", lty = 2)
cbind(R.predictions, JAGS.predictions)

###################################################
###TRABAJAR ESTE EJEMPLO EN CLASE??################
###agregar otros ejemplos?################
###################################################
### GLM Poisson GLM en JAGS y frecuentista
### Datos reales de Halcon peregrino 
### Conteos en el tiempo parejas de adultos, parejas reproductivas y jovenes 
### Modelo 3.3.2 pp. BPA
###################################################
# Read data
peregrine <- read.table("falcons.txt", header = TRUE)

attach(peregrine)
plot(Year, Pairs, type = "b", lwd = 2, main = "", las = 1, ylab = "Pair count", xlab = "Year", ylim = c(0, 200), pch = 16)

# Bundle data
mean.year <- mean(1:length(Year))        # Mean of year covariate
sd.year <- sd(1:length(Year))            # SD of year covariate
win.data <- list(C = Pairs, n = length(Pairs), year = (1: length(Year) - mean.year) / sd.year)

# Initial values
inits <- function() list(alpha = runif(1, -2, 2), beta1 = runif(1, -3, 3))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "lambda")

# MCMC settings
ni <- 5000
nt <- 2
nb <- 1000
nc <- 3

# Call JAGS from R (BRT < 1 min)
out1 <- jags(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.jags", n.chains = nc, 
n.thin = nt, n.iter = ni, n.burnin = nb)

# Summarize posteriors
print(out1, dig = 3) 

JAGS.predictions <- out1$BUGSoutput$mean$lambda
lines(Year, JAGS.predictions, type = "l", lwd = 3, col = "blue", lty = 2)


