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




##############################
#####  GLM Binomial     ######
##############################

# 3.5. Binomial GLM for modeling bounded counts or proportions
# 3.5.1. Generation and analysis of simulated data
nyears = 40 
alpha = 0 
beta1 = -0.1 
beta2 = -0.9
  # nyears: Number of years
  # alpha, beta1, beta2: coefficients
  
  # Generate untransformed and transformed values of time covariate
  year <- 1:nyears
  YR <- (year-round(nyears/2)) / (nyears / 2)
  
  # Generate values of binomial totals (N)
  N <- round(runif(nyears, min = 20, max = 100))
  
  # Signal: build up systematic part of the GLM
  exp.p <- plogis(alpha + beta1 * YR + beta2 * (YR^2))
  
  # Noise: generate random part of the GLM: Binomial noise around expected counts (which is N)
  C <- rbinom(n = nyears, size = N, prob = exp.p)
  
  # Plot simulated data
  plot(year, C/N, type = "b", lwd = 2, col = "black", main = "", las = 1, 
       ylab = "Proporci�n de pares exitosos", xlab = "A�o", ylim = c(0, 1))
  points(year, exp.p, type = "l", lwd = 3, col = "red")
  
  
##### Ejemplo

nyears = 40
alpha = 1
beta1 = -0.03
beta2 = -0.9


###################### Binomial model

# Specify model in BUGS language
sink("GLM_Binomial.txt")
cat("
model {

# Priors
alpha ~ dnorm(0, 0.001)
beta1 ~ dnorm(0, 0.001)
beta2 ~ dnorm(0, 0.001)

# Likelihood
for (i in 1:nyears){
   C[i] ~ dbin(p[i], N[i])          # 1. Distribution for random part
   logit(p[i]) <- alpha + beta1 * year[i] + beta2 * pow(year[i],2) # link function and linear predictor
   }
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = C, N = N, nyears = length(C), year = YR)

# Initial values
inits <- function() list(alpha = runif(1, -1, 1), beta1 = runif(1, -1, 1), beta2 = runif(1, -1, 1))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "p")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call WinBUGS from R (BRT < 1 min)
out <- jags(data = win.data, inits = inits, parameters.to.save = params, 
            model.file = "GLM_Binomial.txt", n.chains = nc, n.thin = nt, n.iter = ni,
            n.burnin = nb)

# Plot predictions
predictions <- out$mean$p
lines(1:length(C), predictions, type = "l", lwd = 3, col = "blue", lty = 2)




###### completar ejemplos GLM

###### ejemplos GLMM random slopes and intercepts


https://github.com/mikemeredith/AHM_code




### GLMM normal

## generar datos
n.groups<-56
n.sample<-10
n<-n.groups*n.sample

pop<-gl(n=n.groups, k=n.sample)

#longitud del cuerpo (cm)
largo.original<-runif(n, 45,70)
mn<-mean(largo.original)
sd<-sd(largo.original)
largo<-(largo.original-mn)/sd
hist(length)

Xmat<-model.matrix(~pop*largo-1-largo)

intercept.mean<-230
intercept.sd<-20
slope.mean<-60
slope.sd<-30

intercept.effects<-rnorm(n=n.groups, mean=intercept.mean, sd=intercept.sd)
slope.effects<-rnorm(n=n.groups, mean=slope.mean, sd=slope.sd)
all.effects<-c(intercept.effects,slope.effects)

lin.pred<-Xmat[,]%*%all.effects
eps<-rnorm(n=n, mean=0, sd=30)
mass<-lin.pred+eps

hist(mass)
library(lattice)
xyplot(mass~largo|pop)


### RANDOM INTERCEPTS, COMMON SLOPE
sink("lme.model1.txt")
cat("
model{

mu.int~dnorm(0,0.001)
tau.int<-1/(sigma.int*sigma.int)
sigma.int~dunif(0,100)

beta~dnorm(0,0.001)
tau<-1/(sigma*sigma)
sigma~dunif(0,100)

#previas
for (i in 1:ngroups){
alpha[i]~dnorm(mu.int, tau.int)
}

#likelihood
for(i in 1:n){
mass[i]~dnorm(mu[i], tau)
mu[i]<- alpha[pop[i]]+beta*largo[i]
}
}
    ", fill=TRUE)
sink()

# datos para jags
jags.data<-list(mass=as.numeric(mass), pop=as.numeric(pop), 
                largo=largo,ngroups=max(as.numeric(pop)), n=n)

# inits 
inits<-function(){list(alpha=rnorm(n.groups,0,2), beta=rnorm(1,1,1), 
                       mu.int=rnorm(1,0,1), sigma.int=rlnorm(1), sigma=rlnorm(1))}

parameters<-c("alpha","beta","mu.int","sigma.int","sigma")

ni<-2000
nb<-500
nt<-2
nc<-3

out<-jags(jags.data, inits, parameters, "lme.model1.txt", n.thin = nt, n.chains = nc,
          n.burnin = nb, n.iter = ni)

largo.pred<-seq(-2,2,,1000)
pred<-array(NA,dim=c(1000,56))

for(i in 1:56){
  pred[,i]<-out$mean$alpha[i]+out$mean$beta*largo.pred
}

matplot(largo.pred,pred, col = "grey", xlab = "largo", ylab = "masa",type="l",lty = 1)
lines(largo.pred, out$mean$mu.int+out$mean$beta*largo.pred, col="black", lwd=3)


### RANDOM INTERCPTS, RANDOM SLOPES
sink("lme.model2.txt")
cat("
model{

mu.int~dnorm(0,0.001)
tau.int<-1/(sigma.int*sigma.int)
sigma.int~dunif(0,100)

mu.slope~dnorm(0,0.001)
tau.slope<-1/(sigma.slope*sigma.slope)
sigma.slope~dunif(0,100)

tau<-1/(sigma*sigma)
sigma~dunif(0,100)

#previas
for (i in 1:ngroups){
alpha[i]~dnorm(mu.int, tau.int)
beta[i]~ dnorm(mu.slope, tau.slope)
}

#likelihood
for(i in 1:n){
mass[i]~dnorm(mu[i], tau)
mu[i]<- alpha[pop[i]]+beta[pop[i]]*largo[i]
}
}
    ", fill=TRUE)
sink()

# datos para jags
jags.data<-list(mass=as.numeric(mass), pop=as.numeric(pop), 
                largo=largo,ngroups=max(as.numeric(pop)), n=n)

# inits 
inits<-function(){list(alpha=rnorm(n.groups,0,2), beta=rnorm(n.groups,1,1), 
                       mu.int=rnorm(1,0,1), sigma.int=rlnorm(1), sigma=rlnorm(1))}

parameters<-c("alpha","beta","mu.int","mu.slope","sigma.int","sigma.slope")

ni<-2000
nb<-500
nt<-2
nc<-3

out<-jags(jags.data, inits, parameters, "lme.model2.txt", n.thin = nt, n.chains = nc,
          n.burnin = nb, n.iter = ni)

largo.pred<-seq(-2,2,,1000)
pred<-array(NA,dim=c(1000,56))

for(i in 1:56){
  pred[,i]<-out$mean$alpha[i]+out$mean$beta[i]*largo.pred
}

matplot(largo.pred,pred, col = "grey", xlab = "largo", ylab = "masa",type="l",lty = 1)
lines(largo.pred, out$mean$mu.int+out$mean$mu.slope*largo.pred, col="black", lwd=3)



