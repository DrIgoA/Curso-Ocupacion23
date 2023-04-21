##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
##########        Ejercicio GLM Poisson              #####
##########               Opcional                    #####
##########################################################

# Los datos presentados corresponden a conteo de artropodos obtenidos en bordes de lotes agricolas con el objetivo de 
# conocer como afectan los diferentes manejos agrícolas (orgánico / Convencional) a la abundancia de artropodos. Se realizaron muestreos
# durante dos años consecutivos. En cada año se relevaron dos estaciones (Primavera / Verano)

rm(list=ls()) #limpio el ambiente de R

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

###################################################################################################
# Ejercicio
# Paso 1: correr el modelo 1. El mismo incorpora al volumen vegetal como una covariable, llamada "volveg". 
# El manejo agrícola (organico / convencional) está puesto como efecto aleatorio. 
######################################################################################
# Modelo 1 = Intercepto Aleatorio / Pendiente Fija
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

save(out_af, file='out_af.rda')
str(out_af)

# Grafico
hist(out_af$summary[278:550,1], xlab="Residuos Pearson", col="grey", breaks=50, main = "")
abline(v=0,col="red",lwd=2)
plot(out_af$summary[5:277,1],out_af$summary[278:550,1], xlab="Valores predichos", ylab="Residuos Pearson", main = "")
abline(h=0,col="red",lwd=2)

###################################################################################################
# Paso 2: modificar el modelo 1 para incorporar el manejo agricola como efecto aleatorio en intercepto y pendiente.
###################################################################################################
# Paso 3: modificar el modelo 1 para obtener el intercepto y la pendiente como fijas.
######################################################################################
# Paso 4: Graficar cada salida de cada modelo. Modificar script de acuerdo a cada modelo
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
