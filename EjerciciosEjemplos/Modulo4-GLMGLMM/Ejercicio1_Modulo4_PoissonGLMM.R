##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
##########       Ejercicio GLMM Poisson              #####
##########################################################

# Los datos presentados corresponden a conteo de artropodos obtenidos en bordes de lotes agricolas con el objetivo de 
# conocer como afectan los diferentes manejos agrícolas (orgánico / Convencional) a la abundancia de artropodos. Se realizaron muestreos
# durante dos años consecutivos. En cada año se relevaron dos estaciones (Primavera / Verano)

rm(list=ls()) #limpio el ambiente de R

library(jagsUI)

data <- read.csv("Datos_Ejercicio1_Modulo4_PoissonGLMM.csv", header=T)
attach(data)
str(data)

# Grafico exploratorio (abundancia de artropodos y volumen vegetal)
plot(volveg,abund)

# seteo sitios
site.list <- data$point2[1:128]
site.name <- seq(32)

# seteo establecimientos. 4 establecimientos: DH, AV, CH, LG
farm <- as.factor(data$farm)
farm <- levels(farm) 

# seteo años. 2 años
#year <- as.factor(data$year)
#year <- levels(year)

nsite <- length(site.name)
nfarm <- length(farm)
#nyear <- length(year)

library(abind) # matrices

###################################################################################################
# Ejercicio
# Paso 1: correr el modelo 1. El mismo incorpora al volumen vegetal como una covariable, llamada "vveg"
# en este caso, las matrices van a ser de dos dimensiones: sitio y establecimiento.
# Análisis del año 1
##  AÑO 1 - lectura de los datos del año 1
data1 <- data[1:128,] # en el csv los datos están ordenados por años, por eso el año 1 va desde la fila 1 a la 204. Fila 1 porque no se cuenta el encabezado

### Construir matrices individuales
# separo cada año x establecimiento
# en este caso, los establecimientos del año 1
dh1 <- subset(data1, data1$farm == 'DH')
av1 <- subset(data1, data1$farm == 'AV')
lg1 <- subset(data1, data1$farm == 'LG')
ch1 <- subset(data1, data1$farm == 'CH')

# abundancia total de artrópodos
abund.dh1 <- dh1$abund 
abund.ch1 <- ch1$abund 
abund.lg1 <- lg1$abund 
abund.av1 <- av1$abund 

abund <- abind(abund.av1, abund.ch1, abund.dh1, abund.lg1, along = 2)
str(abund)
dimnames(abund) <- list(site=site.name, farm=farm)
abund[,1]

# vveg - volumen vegetal (covariable)
vveg.dh1 <- dh1$volveg 
vveg.dh1 <- as.matrix(scale(vveg.dh1))

vveg.ch1 <- ch1$volveg 
vveg.ch1 <- as.matrix(scale(vveg.ch1))

vveg.lg1 <- lg1$volveg 
vveg.lg1 <- as.matrix(scale(vveg.lg1))

vveg.av1 <- av1$volveg 
vveg.av1 <- as.matrix(scale(vveg.av1))

vveg <- abind(vveg.av1, vveg.ch1, vveg.dh1, vveg.lg1, along = 2)
str(vveg)
dimnames(vveg) <- list(site=site.name, farm=farm)
vveg[,1]

#
win.data <- list(abund=abund, 
                 vveg=vveg,
                 nsite=dim(abund)[1],
                 nfarm=dim(abund)[2]
                 )
str(win.data)

# Modelo 1
cat(file = "M4-GLMMPoisson-Ejercicio-1-Facu.txt","
model {
# modelos para missing covariates
    
     for (e in 1:nfarm){ 
       for(j in 1:nsite) {
          
          vveg  [j,e] ~ dnorm(mu.vveg, tau.vveg)
                         }
                       }
                     
                     
          tau.vveg <-pow(sd.vveg,-2)
          sd.vveg ~ dunif(0,1)
          mu.vveg ~ dnorm(0,1)
          
          
# Priors
   mu.alpha ~ dnorm(0, 0.001)  # Mean hyperparam
   tau.alpha <- pow(sd.alpha, -2)
   sd.alpha ~ dunif(0, 10)     # sd hyperparam
 
for(k in 1:1){
  alpha[k] ~ dunif(-10, 10)   # Regression params
}

# Likelihood
for(j in 1:nsite){
 
  alpha0[j] ~ dnorm(mu.alpha, tau.alpha) # Random effects and hyperparams
  re0[j] <- alpha0[j] - mu.alpha         # zero-centered random effects
  
for (e in 1:nfarm){ 
     
 
  abund[j,e] ~ dpois(lambda[j,e])
  log(lambda[j,e]) <- alpha0[j] + alpha[1] * vveg  [j,e] 
                                  
 }
 }
  
}")

# Other model run preparations
inits <- function() list(alpha0 = rnorm(1:32), alpha = rnorm(1)) # Inits
params <- c("mu.alpha", "sd.alpha", "alpha0", "alpha")           # Params
ni <- 300000 ; nt <- 25 ; nb <- 150000 ; nc <- 3                 # MCMC settings

# Call WinBUGS or JAGS from R (6-7 min) and summarize posteriors
out <- jags(win.data, inits, params, "M4-GLMMPoisson-Ejercicio-1-Facu.txt", n.chains = nc, n.thin = nt,
             n.iter = ni, n.burnin = nb)

str(out)
# checkeo de convergencia
print(out, dig=3) 
par(mfrow = c(3,2)) ; traceplot(out, c("mu.alpha", "sd.alpha", "alpha[1:1]"))
traceplot(out,c("mu.alpha","sd.alpha","alpha[1:1]"))
# guardar la salida
save(out, file='out.rda')

# cargar la salida
load('out.rda')

###################################################################################################
# Paso 2: graficar la relación entre el volumen vegetal ("vveg") y la abundancia de artropodos a partir de la salida (out) del GLMM. 

###################################################################################################
# Paso 3: Incorporar otra covariable. En este caso, resulta de interés agregar el manejo agrícola (Orgánico / Convencional) 
# y la estación/temporada de muestreo (Primavera / Verano). Plantee interacciones entre covariables de interés. 

###################################################################################################
# Paso 4: Incorporar el año como efecto aleatorio. Debe agregar una tercera dimension que corresponda al año. 

