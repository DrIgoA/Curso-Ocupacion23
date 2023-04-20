##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########       Ejemplo Ocupación en JAGS      #########
###########              Comunidades             #########
##########################################################
########              Ejemplos de:                  ######
########  Applied hierarchical modeling in ecology  ######
########  Modeling distribution, abundance and      ######
########  species richness using R and BUGS         ######
########  Volume 1: Prelude and Static models       ######
########      Marc Kéry & J. Andy Royle             ######
############################################################

######################################
### Borrar documentacion anterior de R
######################################
rm(list=ls(all=TRUE))

###################################### 
### Elegir directorio de trabajo
######################################
setwd("C:\\Users\\andrea\\Documents\\GitHub\\Curso-Ocupacion23\\EjerciciosEjemplos\\Modulo6-OcupacionComunidad")

######################################
### Cargar Paquetes 
######################################
library(jagsUI)    #paquete JAGS

######################################
### Procesamiento de Datos 
######################################
## Datos del Swiss breeding bird survey MHB
# Leemos los datos y los revisamos
data <- read.csv("MHB_2014.csv", header = T, sep = ",", dec= ".")

str(data)   
head(data)
dim(data)

# en este caso podemos decir que los datos estan ingresados en la BD de manera longitudinal los sitios y 
# las especies, pero las repeticiones están por columnas. Cada especie esta una abajo de la otra

###### Aqui comienzan a procesar los datos y reorganizarlos

# Acá crean varias listas de especies (basadas en los nombres en inglés y ordenes)
species.list <- levels(as.factor(data$engname))   # lista por orden alfabético
spec.name.list <- tapply(data$specid, data$engname, mean) # ID de las especies
spec.id.list <- unique(data$specid)    # lista ID
ordered.spec.name.list <- spec.name.list[order(spec.name.list)] # lista de ID en orden

DET <- cbind(data$count141, data$count142, data$count143)  # Conteos
DET[DET > 1] <- 1       # Todo lo que es mayor a 1, lo convierto en 1 (datos para detección/no detección)

# Ponen los datos de deteccion en en arreglo 3D: sitio x rep x especie
nsite <- 267                    # numero de sitios en Swiss MHB
nrep <- 3                       # numero de repeticiones por temporada
nspec <- length(species.list)   # 67 especies

# armo el array para meter el loop
Y <- array(NA, dim = c(nsite, nrep, nspec))

# Esto es un loop para cada conjunto de datos
# lo que hace es rellenar para cada una de las especies con las filas que correspondan
# podemos probar con i=1 o 158 ((i-1)*nsite+1)) se agrega el +1 porque si no tomaría datos 
# de la especie anterior

for(i in 1:nspec){
   Y[,,i] <- DET[((i-1)*nsite+1):(i*nsite),]
}
dimnames(Y) <- list(NULL, NULL, names(ordered.spec.name.list))


# Chequear los datos de una especie y "llamarlos" del array 3D
which(names(ordered.spec.name.list) == "Common Rosefinch")
tmp <- Y[,,which(names(ordered.spec.name.list) == "Common Rosefinch")]
dim(tmp)

# Distribucion de frecuencias de numero de muestreos por sitio
# aqui veo que 1 sitio no tuvo muestreos
 table(nsurveys <- apply(Y[,,1], 1, function(x) sum(!is.na(x))))

# Con esto identifico cual es el sitio que no tuvo muestreo
(NAsites <- which(nsurveys == 0) )


# Numero de especies observadas por sitio 
tmp <- apply(Y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
sort(C <- apply(tmp, 1, sum))     # Compute and print sorted species counts

par(mfrow=c(1,1))

plot(table(C), xlim = c(0, 40), xlab = "Numero de especies observadas", 
     ylab = "Numero de cuadrantes", frame = F)
abline(v = mean(C, na.rm = TRUE), col = "blue", lwd = 3)

######################################################
########     Modelo de Dorazio-Royle (DR)     ########
########    para ocupacion de comunidades     ########
########        con aumento de datos          ########
########           sin covariables            ########
######################################################

# Colapsar los datos 3D de deteccion/no deteccion a frecuencias de deteccion 2D
Ysum <- apply(Y, c(1,3), sum, na.rm = T) # 
Ysum[NAsites,] <- NA                     # los que tienen sitios NA, ponerles NA

# Aumentar set de datos (DA - data augmentation)
nz <- 30                # Numero de especies potenciales en la superpoblacion
M <- nspec + nz          # Tamaño del data set aumentado ('superpoblacion')
Yaug <- cbind(Ysum, array(0, dim=c(nsite, nz))) # Agregar historias con ceros

# unir y llamar al set de datos
str( win.data <- list(Yaug = Yaug, nsite = nrow(Ysum), nrep = data$nsurvey[1:nsite], M = M, nspec = nspec, nz = nz) )

# Especificar el modelo en lenguaje JAGS
sink("modelDA.txt")
cat("
model {

# Previas de la comunidad
for(k in 1:M){                  # Loop sobre especies
  lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
  lp[k] ~ dnorm(mu.lp, tau.lp)
}

# Hiperprevias de community
omega ~ dunif(0,1)              
mu.lpsi ~ dnorm(0,0.001)        
mu.lp ~ dnorm(0,0.001)          
tau.lpsi <- pow(sd.lpsi, -2)
sd.lpsi ~ dunif(0,5)            
tau.lp <- pow(sd.lp, -2)
sd.lp ~ dunif(0,5)              

# proceso de superpoblacion
for(k in 1:M){
  w[k] ~ dbern(omega)           # indicador de pertenencia a la metacomunidad
}                               

# Modelo ecologico
for(k in 1:M){
  mu.psi[k] <- w[k] * psi[k]    # las especies que no son parte de la comunidad obtienen un cero w
  logit(psi[k]) <- lpsi[k]
  for (i in 1:nsite) {
    z[i,k] ~ dbern(mu.psi[k])
  }
}

# Modelo de obsevacion
for(k in 1:M){
  logit(p[k]) <- lp[k]
  for (i in 1:nsite) {
    mu.p[i,k] <- z[i,k] * p[k]  # las especies que no ocurren tienen 0 en z
    Yaug[i,k] ~ dbin(mu.p[i,k], nrep[i])
  }
}

# Derived quantities
for(k in 1:M){
   Nocc.fs[k] <- sum(z[,k])     # Number of occupied sites among the 267
}
for (i in 1:nsite) {
   Nsite[i] <- sum(z[i,])       # Number of occurring species at each site
}
n0 <- sum(w[(nspec+1):(nspec+nz)]) # Number of unseen species in metacommunity
Ntotal <- sum(w[])              # Total metacommunity size (= nspec + n0)
}
",fill = TRUE)
sink()

# Valores iniciales
wst <- rep(1, nspec+nz)                   # Setear todos a ocurrencia 1
zst <- array(1, dim = c(nsite, nspec+nz)) # lo mimso para z
inits <- function() list(z = zst, w = wst, lpsi = rnorm(n = nspec+nz), lp = rnorm(n = nspec+nz))

# Parametros monitoreados
params <- c("mu.lpsi", "mu.lp", "psi", "p", "Nsite", "Ntotal", "omega", "n0")

# seteos de MCMC
ni <- 2000   ;   nt <- 10  ;   nb <- 1000   ;   nc <- 3;   na <- 2000

# Llamar JAGS de R, chequear convergencia y resumir posteriores
outDA <- jags(win.data, inits, params, "modelDA.txt", n.chains = nc, n.thin = nt, 
              n.iter = ni, n.burnin = nb,n.adapt = na, parallel = TRUE)


par(mfrow = c(2,2)) ; traceplot(outDA, c('mu.lpsi', 'mu.lp'))


print(outDA, dig = 3)

# Guardar los datos de la corrida (no usar asi no se sobre escribe la corrida completa)
# save(outDA, file='outDA.rda')

# Llamar a los datos de la corrida completa
load('outDA.rda')


# Graficar la distribucion posteriors de la riqueza sitio especifica (Nsite)
# Graficar para una selección de sitios
par(mfrow = c(3,3), mar = c(5,4,3,2))
for(i in c(9, 32, 162, 12, 27, 30, 118, 159, 250)){
   plot(table(outDA$sims.list$Nsite[,i]), main = paste("Quadrat", i), 
   xlab = "Riqueza local de especies", ylab = "", frame = F, 
   xlim = c((min(C[i], outDA$sims.list$Nsite[,i], na.rm = T)-2),
   max(outDA$sims.list$Nsite[,i]) ))
   abline(v = C[i], col = "grey", lwd = 4)
}

# Graficar la distribucion posteriors de la riqueza total (Ntotal)
par(mfrow = c(1,1), mar = c(5,4,3,2))
plot(table(outDA$sims.list$Ntotal), main = "", ylab = "", xlab = "Metacomunidad de aves", 
     frame = F) #, xlim = c(144, 245))
abline(v = nspec, col = "grey", lwd = 4)


######################################################
########     Modelo de Dorazio-Royle (DR)     ########
########    para ocupacion de comunidades     ########
########        con aumento de datos          ########
########           con covariables            ########
######################################################

######################################
### Covariables 
######################################

# Llamar a las covariables y estandarizarlas
# Elevacion y cobertura de bosque
orig.ele <- data$elev[1:nsite]
(mean.ele <- mean(orig.ele, na.rm = TRUE))
(sd.ele <- sd(orig.ele, na.rm = TRUE))
ele <- (orig.ele - mean.ele) / sd.ele
orig.forest <- data$forest[1:nsite]
(mean.forest <- mean(orig.forest, na.rm = TRUE))
(sd.forest <- sd(orig.forest, na.rm = TRUE))
forest <- (orig.forest - mean.forest) / sd.forest

# Get survey date and survey duration and standardise both
# Survey date (this is Julian date, with day 1 being April 1)
orig.DAT <- cbind(data$date141, data$date142, data$date143)[1:nsite,]
(mean.date <- mean(orig.DAT, na.rm = TRUE))
(sd.date <- sd(c(orig.DAT), na.rm = TRUE))
DAT <- (orig.DAT - mean.date) / sd.date      # scale
DAT[is.na(DAT)] <- 0                         # impute missings
# Survey duration (in minutes)
orig.DUR <- cbind(data$dur141, data$dur142, data$dur143)[1:nsite,]
(mean.dur <- mean(orig.DUR, na.rm = TRUE))
(sd.dur <- sd(c(orig.DUR), na.rm = TRUE))
DUR <- (orig.DUR - mean.dur) / sd.dur        # scale
DUR[is.na(DUR)] <- 0                         # mean impute missings

# Augment data set: choose one of two different priors on Ntotal

nz <- 30                 
Yaug <- array(0, dim=c(nsite, nrep, nspec+nz)) # array with only zeroes
Yaug[,,1:nspec] <- Y      # copy into it the observed data

# Create same NA pattern in augmented species as in the observed species
missings <- is.na(Yaug[,,1]) # e.g., third survey in high-elevation quads
for(k in (nspec+1):(nspec+nz)){
   Yaug[,,k][missings] <- NA
}

# Bundle and summarize data
str(win.data <- list(Y = Yaug, nsite = dim(Y)[1], nrep = dim(Y)[2], nspec = dim(Y)[3], 
                     nz = nz, M = nspec + nz, ele = ele, forest = forest, DAT = DAT ))


# Specify model in BUGS language
sink("modelDAcov.txt")
cat("
model {

# Priors
omega ~ dunif(0,1)
# Priors for species-specific effects in occupancy and detection
for(k in 1:M){
  lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)    # Hyperparams describe community
  betalpsi1[k] ~ dnorm(mu.betalpsi1, tau.betalpsi1)
  betalpsi2[k] ~ dnorm(mu.betalpsi2, tau.betalpsi2)
  betalpsi3[k] ~ dnorm(mu.betalpsi3, tau.betalpsi3)
  lp[k] ~ dnorm(mu.lp, tau.lp)
  betalp1[k] ~ dnorm(mu.betalp1, tau.betalp1)
  betalp2[k] ~ dnorm(mu.betalp2, tau.betalp2)
}

# Hyperpriors
# For the model of occupancy
mu.lpsi ~ dnorm(0,0.01)
tau.lpsi <- pow(sd.lpsi, -2)
sd.lpsi ~ dunif(0,8)   # as always, bounds of uniform chosen by trial and error
mu.betalpsi1 ~ dnorm(0,0.1)
tau.betalpsi1 <- pow(sd.betalpsi1, -2)
sd.betalpsi1 ~ dunif(0, 4)
mu.betalpsi2 ~ dnorm(0,0.1)
tau.betalpsi2 <- pow(sd.betalpsi2, -2)
sd.betalpsi2 ~ dunif(0,2)
mu.betalpsi3 ~ dnorm(0,0.1)
tau.betalpsi3 <- pow(sd.betalpsi2, -2)
sd.betalpsi3 ~ dunif(0,2)

# For the model of detection
mu.lp ~ dnorm(0,0.1)
tau.lp <- pow(sd.lp, -2)
sd.lp ~ dunif(0, 2)
mu.betalp1 ~ dnorm(0,0.1)
tau.betalp1 <- pow(sd.betalp1, -2)
sd.betalp1 ~ dunif(0,1)
mu.betalp2 ~ dnorm(0,0.1)
tau.betalp2 <- pow(sd.betalp2, -2)
sd.betalp2 ~ dunif(0,1)


# Superpopulation process: Ntotal species sampled out of M available
for(k in 1:M){
   w[k] ~ dbern(omega)
}

# Ecological model for true occurrence (process model)
for(k in 1:M){
  for (i in 1:nsite) {
    logit(psi[i,k]) <- lpsi[k] + betalpsi1[k] * ele[i] + 
      betalpsi2[k] * pow(ele[i],2) + betalpsi3[k] * forest[i]
    mu.psi[i,k] <- w[k] * psi[i,k]
    z[i,k] ~ dbern(mu.psi[i,k])
  }
}

# Observation model for replicated detection/nondetection observations
for(k in 1:M){
  for (i in 1:nsite){
    for(j in 1:nrep){
      logit(p[i,j,k]) <- lp[k] + betalp1[k] * DAT[i,j] + 
        betalp2[k] * pow(DAT[i,j],2) 
      mu.p[i,j,k] <- z[i,k] * p[i,j,k]
      Y[i,j,k] ~ dbern(mu.p[i,j,k])
    }
  }
}

# Derived quantities
#for(k in 1:M){
#   Nocc.fs[k] <- sum(z[,k])       # Number of occupied sites among the 267
#}
for (i in 1:nsite){
   Nsite[i] <- sum(z[i,])          # Number of occurring species at each site
}
n0 <- sum(w[(nspec+1):(nspec+nz)]) # Number of unseen species
Ntotal <- sum(w[])                 # Total metacommunity size

# Vectors to save (S for ‘save’; discard posterior samples for 
# all minus 1 of the potential species to save disk space)
# we do this for nz = 250 (i.e., M = 395)
lpsiS[1:(nspec+1)] <- lpsi[1:(nspec+1)]
betalpsi1S[1:(nspec+1)] <- betalpsi1[1:(nspec+1)]
betalpsi2S[1:(nspec+1)] <- betalpsi2[1:(nspec+1)]
betalpsi3S[1:(nspec+1)] <- betalpsi3[1:(nspec+1)]
lpS[1:(nspec+1)] <- lp[1:(nspec+1)]
betalp1S[1:(nspec+1)] <- betalp1[1:(nspec+1)]
betalp2S[1:(nspec+1)] <- betalp2[1:(nspec+1)]
}
",fill = TRUE)
sink()


# Initial values
wst <- rep(1, nspec+nz)                   # Simply set everybody at occurring
zst <- array(1, dim = c(nsite, nspec+nz)) # ditto
inits <- function() list(z = zst, w = wst, lpsi = rnorm(n = nspec+nz), 
                         betalpsi1 = rnorm(n = nspec+nz), betalpsi2 = rnorm(n = nspec+nz), 
                         betalpsi3 = rnorm(n = nspec+nz), lp = rnorm(n = nspec+nz), 
                         betalp1 = rnorm(n = nspec+nz), betalp2 = rnorm(n = nspec+nz))

# Set 1
params1 <- c("omega", "mu.lpsi", "mu.betalpsi1", "mu.betalpsi2","mu.betalpsi3", "mu.lp", 
             "mu.betalp1", "mu.betalp2", "Ntotal", "Nsite")

# MCMC settings
ni <- 40000   ;   nt <- 10   ;   nb <- 10000   ;   nc <- 3

# Run JAGS, check convergence and summarize posteriors
outDAcov <- jags(win.data, inits, params1, "modelDAcov.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

par(mfrow = c(2, 2))
traceplot(outDAcov, c(c("omega", "mu.lpsi", "mu.betalpsi1", 
                      "mu.betalpsi2", "mu.betalpsi3", "mu.lp", 
                      "mu.betalp1", "mu.betalp2","Ntotal")) )

# Set 2
params2 <- c("mu.lpsi", "sd.lpsi", "mu.betalpsi1", "mu.betalpsi2", "mu.betalpsi3", 
             "lpsi", "betalpsi1", "betalpsi2", "betalpsi3", "lp", "sd.lp", "betalp1", "betalp2", 
              "z", "w")
ni <- 400   ;   nt <- 20   ;   nb <- 100   ;   nc <- 3
outDAcov2 <- jags(win.data, inits, params2, "modelDAcov.txt", n.chains = nc, n.thin = nt, 
                  n.iter = ni, n.burnin = nb, parallel = TRUE)

# Guardar los datos de la corrida (no usar asi no se sobre escribe la corrida completa)
# save(outDAcov, file='outDAcov.rda')
save(outDAcov2, file='outDAcov.rda')

# Llamar a los datos de la corrida completa
load('outDAcov.rda')
load('outDAcov2.rda')

#### AGREGAR COVARIABLE DUR EN DETECCION? COMO EJERCICIO

# Visualize covariate mean relationships for the average species
o.ele <- seq(200, 2500,,500)               # Get covariate values for prediction
o.for <- seq(0, 100,,500)
o.dat <- seq(15, 120,,500)
o.dur <- seq(100, 420,,500)
ele.pred <- (o.ele - mean.ele) / sd.ele
for.pred <- (o.for - mean.forest) / sd.forest
dat.pred <- (o.dat - mean.date) / sd.date
dur.pred <- (o.dur - mean.dur) / sd.dur

# Predict occupancy for elevation and forest and detection for date and duration
# Put all fourpredictions into a single
str( tmp <- outDA$sims.list )              # grab MCMC samples
nsamp <- length(tmp[[1]])    # number of mcmc samples
predC <- array(NA, dim = c(500, nsamp, 3)) # "C" for 'community mean'


for(i in 1:nsamp){
   predC[,i,1] <- plogis(tmp$mu.lpsi[i] + tmp$mu.betalpsi1[i] * ele.pred + 
     tmp$mu.betalpsi2[i] * ele.pred^2 )
   #predC[,i,2] <- plogis(tmp$mu.lpsi[i] + tmp$mu.betalpsi3[i] * for.pred)
   #predC[,i,3] <- plogis(tmp$mu.lp[i] + tmp$mu.betalp1[i] * dat.pred + 
    # tmp$mu.betalp2[i] * dat.pred^2 )
   }

# Get posterior means and 95% CRIs and plot (Fig. 11–17)
pmC <- apply(predC, c(1,3), mean)
criC <- apply(predC, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975)))

par(mfrow = c(2, 2))
plot(o.ele, pmC[,1], col = "blue", lwd = 3, type = 'l', lty = 1, frame = F, ylim = c(0, 0.05), xlab = "Elevation (m a.s.l)", ylab = "Community mean occupancy")
matlines(o.ele, t(criC[,,1]), col = "grey", lty = 1)
plot(o.for, pmC[,2], col = "blue", lwd = 3, type = 'l', lty = 1, frame = F, ylim = c(0, 0.05), xlab = "Forest cover", ylab = "Community mean occupancy")
matlines(o.for, t(criC[,,2]), col = "grey", lty = 1)
plot(o.dat, pmC[,3], col = "blue", lwd = 3, type = 'l', lty = 1, frame = F, ylim = c(0.2, 0.8), xlab = "Survey date", ylab = "Community mean detection")
matlines(o.dat, t(criC[,,3]), col = "grey", lty = 1)
plot(o.dur, pmC[,4], col = "blue", lwd = 3, type = 'l', lty = 1, frame = F, ylim = c(0.2, 0.8), xlab = "Survey duration", ylab = "Community mean detection")
matlines(o.dur, t(criC[,,4]), col = "grey", lty = 1)


# Plot posterior distribution of site-specific species richness (Nsite)
par(mfrow = c(3,3), mar = c(5,4,3,2))
for(i in 1:267){
   plot(table(out10$sims.list$Nsite[,i]), main = paste("Quadrat", i), 
   xlab = "Local species richness", ylab = "", frame = F, 
   xlim = c((min(C[i], out10$sims.list$Nsite[,i], na.rm = T)-2),
   max(out10$sims.list$Nsite[,i]) ))
   abline(v = C[i], col = "grey", lwd = 4)
   browser()
}

# Plot it only for a selection of sites (Fig. 11-18)
par(mfrow = c(3,3), mar = c(5,4,3,2))
for(i in c(9, 32, 162, 12, 27, 30, 118, 159, 250)){
   plot(table(out10$sims.list$Nsite[,i]), main = paste("Quadrat", i), 
   xlab = "Local species richness", ylab = "", frame = F, 
   xlim = c((min(C[i], out10$sims.list$Nsite[,i], na.rm = T)-2),
   max(out10$sims.list$Nsite[,i]) ))
   abline(v = C[i], col = "grey", lwd = 4)
}

# Plot Nsite estimates under models 9 & 10 vs. elevation (Fig. 11-19)
offset <- 30    # Set off elevation for better visibility
plot(elev, out9$mean$Nsite, xlab = "Elevation (metres)", ylab = "Community size estimate (Nsite)", frame = F, ylim = c(0,60), pch = 16) # black: model 9
lines(smooth.spline(out9$mean$Nsite ~ elev), lwd = 3)
points(elev+offset, out10$mean$Nsite, pch = 16, col = "blue") # red: model 10
lines(smooth.spline(out10$mean$Nsite ~ elev), lwd = 3, col = "blue")


str(all10)                    # look at the MCMC output
pm <- apply(all10, 2, mean)    # Get posterior means and 95% CRIs
cri <- apply(all10, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs


# Effects of date (linear and quadratic) and of duration on detection
#par(mfrow = c(1,3), cex.lab = 1.3, cex.axis = 1.3) # Can put all three in one
par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
# Date linear (Fig. 11 – 20 left)
plot(pm[1:145], 1:145, xlim = c(-1.5, 1.5), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of date (linear) on detection", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 1:145], 1:145, cri[2, 1:145], 1:145, col = "grey", lwd = 1)
sig1 <- (cri[1, 1:145] * cri[2, 1:145]) > 0
segments(cri[1, 1:145][sig1 == 1], (1:145)[sig1 == 1], cri[2, 1:145][sig1 == 1], (1:145)[sig1 == 1], col = "blue", lwd = 2)
abline(v = out101$summary[11,1], lwd = 3, col = "red")
abline(v = out101$summary[11,c(3,7)], lwd = 2, col = "red", lty = 2)



# Date quadratic (not shown)
plot(pm[216:360], 1:145, xlim = c(-1.5, 1.5), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of date (quadratic) on detection", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 216:360], 1:145, cri[2, 216:360], 1:145, col = "grey", lwd = 1)
sig2 <- (cri[1, 216:360] * cri[2, 216:360]) > 0
segments(cri[1, 216:360][sig2 == 1], (1:145)[sig2 == 1], cri[2, 216:360][sig2 == 1], (1:145)[sig2 == 1], col = "blue", lwd = 2)
abline(v = out101$summary[13,1], lwd = 3, col = "red")
abline(v = out101$summary[13, c(3,7)], lwd = 3, col = "red", lty = 2)


# Survey duration (Fig. 11-20 right)
plot(pm[431:575], 1:145, xlim = c(-0.5, 1), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of survey duration on detection", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 431:575], 1:145, cri[2, 431:575], 1:145, col = "grey", lwd = 1)
sig3 <- (cri[1, 431:575] * cri[2, 431:575]) > 0
segments(cri[1, 431:575][sig3 == 1], (1:145)[sig3 == 1], cri[2, 431:575][sig3 == 1], (1:145)[sig3 == 1], col = "blue", lwd = 2)
abline(v = out101$summary[15,1], lwd = 3, col = "red")
abline(v = out101$summary[15, c(3,7)], lwd = 3, col = "red", lty = 2)


# Effects of elevation (linear and quadratic) and of forest on occupancy
# par(mfrow = c(1,3), cex.lab = 1.3, cex.axis = 1.3) # can do all in one
# Effect of elevation (linear) on occupancy probability (Fig. 11-21)
plot(pm[646:790], 1:145, xlim = c(-8, 8), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of elevation (linear) on occupancy", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 646:790], 1:145, cri[2, 646:790], 1:145, col = "grey", lwd = 1)
sig4 <- (cri[1, 646:790] * cri[2, 646:790]) > 0
segments(cri[1, 646:790][sig4 == 1], (1:145)[sig4 == 1], cri[2, 646:790][sig4 == 1], (1:145)[sig4 == 1], col = "blue", lwd = 2)
abline(v = out101$summary[3,1], lwd = 3, col = "red")
abline(v = out101$summary[3,c(3,7)], lwd = 3, col = "red", lty = 2)


# Effect of elevation (quadratic) on occupancy probability (Fig. 11-22)
plot(pm[861:1005], 1:145, xlim = c(-4, 2), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of elevation (quadratic) on occupancy", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 861:1005], 1:145, cri[2, 861:1005], 1:145, col = "grey", lwd=1)
sig5 <- (cri[1, 861:1005] * cri[2, 861:1005]) > 0
segments(cri[1, 861:1005][sig5 == 1], (1:145)[sig5 == 1], cri[2, 861:1005][sig5 == 1], (1:145)[sig5 == 1], col = "blue", lwd = 2)
abline(v = out101$summary[5,1], lwd = 3, col = "red")
abline(v = out101$summary[5,c(3,7)], lwd = 3, col = "red", lty = 2)


# Effect of forest (linear) on occupancy probability (Fig. 11-23)
plot(pm[1076:1220], 1:145, xlim = c(-3, 4), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of forest cover on occupancy", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 1076:1220], 1:145, cri[2, 1076:1220],1:145, col = "grey", lwd=1)
sig6 <- (cri[1, 1076:1220] * cri[2, 1076:1220]) > 0
segments(cri[1, 1076:1220][sig6 == 1], (1:145)[sig6 == 1], cri[2, 1076:1220][sig6 == 1], (1:145)[sig6 == 1], col = "blue", lwd = 2)
abline(v = out101$summary[7,1], lwd = 3, col = "red")
abline(v = out101$summary[7,c(3,7)], lwd = 3, col = "red", lty = 2)
negsig6 <- (cri[1, 1076:1220] < 0 & cri[2, 1076:1220] < 0) == 1 # sig negative
possig6 <- (cri[1, 1076:1220] > 0 & cri[2, 1076:1220] > 0) == 1 # sig positive


# Predict detection for date and duration and occupancy for elevation and forest
# for each of the 145 observed species
predS <- array(NA, dim = c(500, nspec, 4))   # covariate value x species x response, "S" for 'species'
p.coef <- cbind(lp=pm[1292:1436], betalp1 = pm[1:145], betalp2 = pm[216:360], betalp3 = pm[431:575])
psi.coef <- cbind(lpsi=pm[1507:1651], betalpsi1 = pm[646:790], betalpsi2 = pm[861:1005], betalpsi3 = pm[1076:1220])

for(i in 1:nspec){          # Loop over 145 observed species
   predS[,i,1] <- plogis(p.coef[i,1] + p.coef[i,2] * dat.pred + 
     p.coef[i,3] * dat.pred^2 )     # p ~ date
   predS[,i,2] <- plogis(p.coef[i,1] + p.coef[i,4] * dur.pred) # p ~ duration
   predS[,i,3] <- plogis(psi.coef[i,1] + psi.coef[i,2] * ele.pred + 
     psi.coef[i,3] * ele.pred^2 )     # psi ~ elevation
   predS[,i,4] <- plogis(psi.coef[i,1] + psi.coef[i,4] * for.pred) # psi ~ forest
}

# Plots for detection probability and survey date and duration (Fig. 11-24)
par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
plot(o.dat, predS[,1,1], lwd = 3, type = 'l', lty = 1, frame = F, 
   ylim = c(0, 1), xlab = "Survey date (1 = 1 April)", 
   ylab = "Detection probability")
for(i in 2:145){
   lines(o.dat, predS[,i,1], col = i, lwd = 3)
}

plot(o.dur, predS[,1,2], lwd = 3, type = 'l', lty = 1, frame = F, 
   ylim = c(0, 1), xlab = "Survey duration (min)", 
   ylab = "Detection probability")
for(i in 2:145){
   lines(o.dur, predS[,i,2], col = i, lwd = 3)
}


# Plots for occupancy probability and elevation and forest cover (Fig. 11-25)
par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
plot(o.ele, predS[,1,3], lwd = 3, type = 'l', lty = 1, frame = F, 
   ylim = c(0, 1), xlab = "Elevation (m a.s.l.)", 
   ylab = "Occupancy probability")
for(i in 2:145){
   lines(o.ele, predS[,i,3], col = i, lwd = 3)
}

plot(o.for, predS[,1,4], lwd = 3, type = 'l', lty = 1, frame = F, 
   ylim = c(0, 1), xlab = "Forest cover (%)", ylab = "Occupancy probability")
for(i in 2:145){
   lines(o.for, predS[,i,4], col = i, lwd = 3)
}

