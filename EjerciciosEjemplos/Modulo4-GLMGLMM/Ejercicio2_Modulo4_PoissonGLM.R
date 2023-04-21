rm(list=ls())

#### Ejemplo de Poisson GLM

library(jagsUI)
data <- read.csv("datos_poisson_GLM.csv", header=T)
attach(data)
str(data)

# grafico exploratorio: volumen vegetal ("volveg") en relación a la abundancia total de insectos ("total")
plot(volveg,total)

# Manejo agrícola: efecto fijo
manejo <- as.numeric(manejo) # 1= manejo orgánico / 2= manejo convencional
table(manejo)
str(manejo)

######################################################################################
# Modelo 1 = Intercepto VARIABLE / Pendiente FIJA  ---- No hay interacción con manejo
win.data <- list(total=total, M=length(total), volveg=volveg, manejo=manejo, e=0.0001)
str(win.data)
cat(file="af.txt", "
model {
#Priors
for (k in 1:2){
  alpha[k] ~ dnorm(0,1.0E-06)
  beta [k] ~ dnorm(0,1.0E-06)
}
#Likelihood
for (i in 1:M){
total[i] ~ dpois(lambda[i])
log(lambda[i]) <- alpha[manejo[i]] + beta[1]*volveg[i]
resi[i] <- (total[i]-lambda[i])/(sqrt(lambda[i])+e)
}
}
    
    ")

inits <- function() list(alpha = rnorm(2,,3), beta = rnorm(2,,3))

# Parameters
params <- c("alpha", "beta", "lambda", "resi")
# MCMC settings
ni <- 10000 ; nt <- 1 ; nb <- 1000 ; nc <- 3
# Call JAGS from R and summarize posteriors
out_af <- jags(win.data, inits, params, "af.txt", 
               n.chains = nc, n.thin = nt, n.iter= ni, n.burnin = nb)
print(out_af)

#save(out_af, file='out_af.rda')
str(out_af)

#
hist(out_af$summary[278:550,1], xlab="Residuos Pearson", col="grey", breaks=50, main = "")
abline(v=0,col="red",lwd=2)
plot(out_af$summary[5:277,1],out_af$summary[278:550,1], xlab="Valores predichos", ylab="Residuos Pearson", main = "")
abline(h=0,col="red",lwd=2)

######################################################################################
# Modelo 2 = Intercepto ALEATORIO / Pendiente ALEATORIA  --- hay interacción con manejo
win.data <- list(total=total, M=length(total), volveg=volveg, manejo=manejo, e=0.0001)
str(win.data)
cat(file="aa.txt", "
model {
#Priors
for (k in 1:2){
  alpha[k] ~ dnorm(0,1.0E-06)
  beta [k] ~ dnorm(0,1.0E-06)
}
#Likelihood
for (i in 1:M){
total[i] ~ dpois(lambda[i])
log(lambda[i]) <- alpha[manejo[i]] + beta[manejo[i]]*volveg[i]
resi[i] <- (total[i]-lambda[i])/(sqrt(lambda[i])+e)
}
}
    
    ")

inits <- function() list(alpha = rnorm(2,,3), beta = rnorm(2,,3))

# Parameters
params <- c("alpha", "beta", "lambda", "resi")
# MCMC settings
ni <- 10000 ; nt <- 1 ; nb <- 1000 ; nc <- 3
# Call JAGS from R and summarize posteriors
out_aa <- jags(win.data, inits, params, "aa.txt", 
               n.chains = nc, n.thin = nt, n.iter= ni, n.burnin = nb)
print(out_aa)

#
hist(out_aa$summary[278:550,1], xlab="Residuos Pearson", col="grey", breaks=50, main = "")
abline(v=0,col="red",lwd=2)
plot(out_aa$summary[5:277,1],out_aa$summary[278:550,1], xlab="Valores predichos", ylab="Residuos Pearson", main = "")
abline(h=0,col="red",lwd=2)

#save(out_aa, file='out_aa.rda')


######################################################################################
# Modelo 4 = Intercepto FIJO / Pendiente FIJA
win.data <- list(total=total, M=length(total), volveg=volveg, manejo=manejo, e=0.0001)
str(win.data)
cat(file="ff.txt", "
model {
#Priors
for (k in 1:2){
  alpha[k] ~ dnorm(0,1.0E-06)
  beta [k] ~ dnorm(0,1.0E-06)
}
#Likelihood
for (i in 1:M){
total[i] ~ dpois(lambda[i])
log(lambda[i]) <- alpha[1] + beta[1]*volveg[i]
resi[i] <- (total[i]-lambda[i])/(sqrt(lambda[i])+e)
}
}
    
    ")

inits <- function() list(alpha = rnorm(2,,3), beta = rnorm(2,,3))

# Parameters
params <- c("alpha", "beta", "lambda", "resi")
# MCMC settings
ni <- 10000 ; nt <- 1 ; nb <- 1000 ; nc <- 3
# Call JAGS from R and summarize posteriors
out_ff <- jags(win.data, inits, params, "ff.txt", 
               n.chains = nc, n.thin = nt, n.iter= ni, n.burnin = nb)
print(out_ff)

#save(out_ff, file='out_ff.rda')

#
hist(out_ff$summary[278:550,1], xlab="Residuos Pearson", col="grey", breaks=50, main = "")
abline(v=0,col="red",lwd=2)
plot(out_ff$summary[5:277,1],out_ff$summary[278:550,1], xlab="Valores predichos", ylab="Residuos Pearson", main = "")
abline(h=0,col="red",lwd=2)

######################################################################################
# Gráficos - modificar de acuerdo a modelo
#tiff(file = "ff.tiff",                #Guardados como Tiff
#     width = 90, height = 100,
#     units = "mm", res = 600)

par(mfrow=c(2,2),mar=c(2.5,2.5,0.2,0.2))

sorted.volveg1 <- sort(volveg[manejo==1])
sorted.y1 <- out_ff$summary[5:277,][manejo==1,][order(volveg[manejo==1]),]
sorted.volveg2 <- sort(volveg[manejo==2])
sorted.y2 <- out_ff$summary[5:277,][manejo==2,][order(volveg[manejo==2]),]

plot(volveg[manejo==1],total[manejo==1], ylab="", 
     xlab="", ylim=c(0,500), xlim=c(0,1), yaxt="n", xaxt="n", col="darkgreen")
points(volveg[manejo==2],total[manejo==2], col="darkred")
axis(1,c(0,0.5,1), padj=-2, tck=0.02, cex = 0.8, cex.axis=0.7)
axis(2,c(0,100,200,300,400,500),las=1, hadj=0.5, tck=0.02, cex.axis=0.7)
mtext("Abundancia total de artropodos", font=1,   at=250, side=2, line = 1.5, cex=1, padj= -0.2)
mtext("VolVeg", at=1, side=1, line = 1, cex = 1, padj= 0.1, font=1)
lines(sorted.volveg1,sorted.y1[,1], col="forestgreen", lwd=2) #Media posterior
lines(sorted.volveg2,sorted.y2[,1], col="red", lwd=2) #Media posterior

#dev.off()


sorted.volveg1 <- sort(volveg[manejo==1])
sorted.y1 <- out_af$summary[5:277,][manejo==1,][order(volveg[manejo==1]),]
sorted.volveg2 <- sort(volveg[manejo==2])
sorted.y2 <- out_af$summary[5:277,][manejo==2,][order(volveg[manejo==2]),]

plot(volveg[manejo==1],total[manejo==1], ylab="", 
     xlab="", ylim=c(0,500), xlim=c(0,1), yaxt="n", xaxt="n", col="darkgreen")
points(volveg[manejo==2],total[manejo==2], col="darkred")
axis(1,c(0,0.5,1), padj=-2, tck=0.02, cex = 0.8, cex.axis=0.7)
axis(2,c(0,100,200,300,400,500),las=1, hadj=0.5, tck=0.02, cex.axis=0.7)
mtext("Abundancia total de artrópodos", font=1,   at=250, side=2, line = 1.5, cex=1, padj= -0.2)
mtext("VolVeg", at=1, side=1, line = 1, cex = 1, padj= 0.1, font=1)
lines(sorted.volveg1,sorted.y1[,1], col="forestgreen", lwd=2) #Media posterior
lines(sorted.volveg2,sorted.y2[,1], col="red", lwd=2) #Media posterior


sorted.volveg1 <- sort(volveg[manejo==1])
sorted.y1 <- out_aa$summary[5:277,][manejo==1,][order(volveg[manejo==1]),]
sorted.volveg2 <- sort(volveg[manejo==2])
sorted.y2 <- out_aa$summary[5:277,][manejo==2,][order(volveg[manejo==2]),]

plot(volveg[manejo==1],total[manejo==1], ylab="", 
     xlab="", ylim=c(0,500), xlim=c(0,1), yaxt="n", xaxt="n", col="darkgreen")
points(volveg[manejo==2],total[manejo==2], col="darkred")
axis(1,c(0,0.5,1), padj=-2, tck=0.02, cex = 0.8, cex.axis=0.7)
axis(2,c(0,100,200,300,400,500),las=1, hadj=0.5, tck=0.02, cex.axis=0.7)
mtext("Abundancia total de artrópodos", font=1,   at=250, side=2, line = 1.5, cex=1, padj= -0.2)
mtext("VolVeg", at=1, side=1, line = 1, cex = 1, padj= 0.1, font=1)
lines(sorted.volveg1,sorted.y1[,1], col="forestgreen", lwd=2) #Media posterior
lines(sorted.volveg2,sorted.y2[,1], col="red", lwd=2) #Media posterior