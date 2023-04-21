##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########       Ejemplo Ocupación en JAGS      #########
###########              1 especie               #########
##########################################################
#############################################################
#         Ejercicio basado en                               #  
# 	S Dardanelli; NC Calamari; SB Canavelli; FR Barzan;     #
#   AP Goijman; L Lezana. Vegetation structure and livestock# 
#   grazing intensity affect ground-foraging birds in       #
#   xerophytic forests of Central-East Argentina. Forest    # 
#   Ecology and Management. 2022 Volume 521.                #
#   https://doi.org/10.1016/j.foreco.2022.120439            #
#  *** Aclaracion:  Los datos originales fueron levemente   #
#                   modificados                             #
#############################################################

rm(list=ls(all=TRUE))
library(jagsUI)

setwd('C:/Users/andrea/Documents/GitHub/Curso-Ocupacion23/EjerciciosEjemplos/Modulo5-OcupacionSimple')

#30 predios, 2 años, 6 puntos por predio y 2 repeticiones
data <- read.csv("SICFLA.csv",sep ="," , header = T)
str(data)
head(data)

#cargo los paquetes que voy a usar
library(reshape)
library(plyr)

# reacomodar los datos
data.melt=melt(data, id.var=c("Predio", "Anio", "NPunto", "Rep"), measure.var="SICFLA")
y=cast(data.melt, Predio ~ Anio ~ NPunto ~ Rep)

# Verifico las dimensiones
dim(y)

# Si trabajo con occupancy paso todos los numero a 0 o 1
y<-aaply(y,1,function(x) {x[x>1]<-1;x})

###################################################################################################
###################################################################################################
##  COVARIABLES no anuales

cova <- read.csv("cova_sicfla.csv", sep= ",", header = T)
head(cova)
str(cova)

# estandarizo covariables
cova$herbz<-scale(cova$C_Herb)
cova$darbz<-scale(cova$D_Arb)
cova$darbuz<-scale(cova$D_Arbust)

covs<-cbind(cova$herbz, cova$darbz, cova$darbuz)
pairs(covs)

# Cobertura Herbacea
herbz.melt=melt(cova, id.var=c("Predio", "NPunto"), measure.var="herbz")
herb=cast(herbz.melt, Predio ~ NPunto)
herb<- as.matrix(herb)

# Densidad Arbol
darbz.melt=melt(cova, id.var=c("Predio", "NPunto"), measure.var="darbz")
darb=cast(darbz.melt, Predio ~ NPunto)
darb<- as.matrix(darb)

# Densidad Arbust
darbuz.melt=melt(cova, id.var=c("Predio", "NPunto"), measure.var="darbuz")
darbu=cast(darbuz.melt, Predio ~ NPunto)
darbu<- as.matrix(darbu)


########################
##  COVARIABLES anuales - fueron medidas los dos anios

covat <- read.csv("cova_sicflayear.csv", sep= ",", header = T)
head(covat)
str(covat)

# estandarizo esta cova
covat$ipz<-scale(covat$IP)

# Intensidad de pastoreo
ip.melt=melt(covat, id.var=c("Predio", "NPunto", "Anio"), measure.var="ipz")
ip=cast(ip.melt, fun.aggregate = mean, Predio ~ Anio ~ NPunto)


###################################################################################################
###################################################################################################

dim(y)

str(win.data <- list(y = y, 
                     cov1 = ip,                                      
                     cov2 = herb,
                     cov4 = darb, 
                     nsite = dim(y)[3], 
                     nrep =  dim(y)[4], 
                     nyear=  dim(y)[2],
                     npredio= dim(y)[1]))

#modelo en lenguaje BUGS
sink("singlesp_bal_nested_time_det.txt")
cat("
    model {

    # modelos para missing covariates
    for (e in 1:npredio){                                       
    for (i in 1:nsite) {                                       
    cov2  [e,i] ~ dnorm(mu.cov2,   tau.cov2)
    cov4  [e,i] ~ dnorm(mu.cov4,   tau.cov4)
    }
    for (t in 1:nyear){                                         
    for (i in 1:nsite) { 
    cov1  [e,t,i] ~ dnorm(mu.cov1,   tau.cov1)
    }
    }
    }
    
    tau.cov1 <-pow(sd.cov1,-2)
    sd.cov1 ~ dunif(0,10)
    mu.cov1 ~ dnorm(0,0.01)
    
    tau.cov2 <-pow(sd.cov2,-2)
    sd.cov2 ~ dunif(0,10)
    mu.cov2 ~ dnorm(0,0.01)
    
    tau.cov4 <-pow(sd.cov4,-2)
    sd.cov4 ~ dunif(0,10)
    mu.cov4 ~ dnorm(0,0.01)
    
    #Hyperprevias
    mu.lp           ~ dnorm(0, 1/2.25^2)
    mu.lpsi      ~ dnorm(0, 1/2.25^2)
    betalpsi1 ~ dnorm (0, 1/2.25^2)    
    betalpsi2 ~ dnorm (0, 1/2.25^2)    
    betalpsi4 ~ dnorm (0, 1/2.25^2)    

    tau.lp        <- pow(sd.lp, -2)
    tau.lpsi      <- pow(sd.lpsi, -2)
    sd.lp  ~ dunif(0,8)
    sd.lpsi  ~ dunif(0,8)
    
    # priors efecto tiempo
    for (t in 1:nyear){ 
    lp [t] ~ dnorm (mu.lp, tau.lp)
    mu.t.lpsi [t] ~ dnorm (mu.lpsi, tau.lpsi)
    tau.t.lpsi[t] <- pow(sd.t.lpsi[t],-2)
    sd.t.lpsi [t] ~ dunif(0,1)
    }
    
    # priors para el efecto del predio en el tiempo
    for (e in 1:npredio){ 
    for (t in 1:nyear){ 
    lpsi  [e,t] ~ dnorm (mu.t.lpsi [t], tau.t.lpsi [t])
    } 
    } 
    
    # Ecological model, process model (true occurrence at site i) 
    ## establishment as random effect
    for (e in 1:npredio){                                       
    for (t in 1:nyear){                                         
    for (i in 1:nsite) {                                       
    
    logit(psi[e,t,i]) <- lpsi[e,t] + betalpsi1 * cov1 [e,t,i] + betalpsi2 * cov2 [e,i]
    + betalpsi4 *cov4 [e,i]

        z[e,t,i] ~ dbern(psi[e,t,i])
    zSim[e,t,i] ~ dbern(psi[e,t,i]) 
    } #nsite
    } #nyear
    } #npredio
    
    # Observation model 
    
    for (e in 1:npredio) {                                    
    for (t in 1:nyear)   {                                    
    for (i in 1:nsite)   {                                       
    for (k in 1:nrep)    {        
    logit(p[e,t,i,k]) <- lp[t] 
    mu.p[e,t,i,k] <- z[e,t,i] * p[e,t,i,k]  
    y[e,t,i,k] ~ dbern(mu.p[e,t,i,k])        
    
# GOF assessment
    Tobs0[e,t,i,k] <- (sqrt(y[e,t,i,k]) - sqrt(p[e,t,i,k]*z[e,t,i]))^2  # FT discrepancy for observed data
    ySim[e,t,i,k] ~ dbern(p[e,t,i,k]*zSim[e,t,i])
    Tsim0[e,t,i,k] <- (sqrt(ySim[e,t,i,k]) - sqrt(p[e,t,i,k]*zSim[e,t,i]))^2  # ...and for simulated data

    }
    }
    }
    } 
# GOF assessment
  Tobs <- sum(Tobs0)
  Tsim <- sum(Tsim0)
#  Derived variable
  N <- sum(z)

    } #model
    ",fill=TRUE)
sink()


# valores iniciales
zst <- apply(y,c(1,2,3),max)
zst[is.na(zst)]<-1
str(zst)

inits <- function(){list(z=zst)}

# parametros a monitorear
# species-specific

params1 <- c("lpsi","lp", "betalpsi1","betalpsi2","betalpsi4","N", "Tobs", "Tsim")

# version larga que ya  esta guardada, pueden correr menos para probar
ni <- 10000 #50000      
nt <- 10  
nb <- 1000 
na<- 5000 #20000     
nc <- 3

out = jags(win.data, inits, params1, "singlesp_bal_nested_time_det.txt", n.chains=nc, 
           n.iter=ni, n.burnin=nb, n.thin=nt, n.adapt = na)

# posterior predictive check
# se hace para evaluar si el modelo es correcto y hubo bondad de ajuste
# no esta del todo aceptado por algunos, pero no hay muchas otras opciones
# se comparan datos simulados versus los observados
pp.check(out, "Tobs", "Tsim", main="Freeman-Tukey discrepancy")

plot(out)

#save(out, file='SICFLA.rda')

load('SICFLA.rda')
print(out, dig=3)

# opcional para guardar la tabla resultados
#write.csv(out$summary, "SICFL_summary.csv")

plot(out)
pp.check(out, "Tobs", "Tsim") #, main="Freeman-Tukey discrepancy")

##############################################################################
#            GRAFICOS
##############################################################################
#cov1 = ip,
#cov1 = herb,                 
#cov3 = darbu,
#cov4 = darb,
#cov5 = dap,

# crear funcion para los CRI 95%
cri<-function(x) quantile(x,prob=c(0.025,0.975))

# means and CRI of each coefficient
betas<-c(mean(out$sims.list$betalpsi1),mean(out$sims.list$betalpsi2),
         mean(out$sims.list$betalpsi4))
cri<-c(cri(out$sims.list$betalpsi1),cri(out$sims.list$betalpsi2),
       cri(out$sims.list$betalpsi4))
lowcri<-c(cri[1],cri[3],cri[5])
upcri<-c(cri[2],cri[4],cri[6])

# compile parameters & names
parameters=c("betalpsi1","betalpsi2","betalpsi4")
names= c('ip', 'herb', 'darb')

# plot Beta coefficients
# tiff(file = "sicla_coeff.tiff", width = 200, height = 200, units = "mm", res = 300)
par(mfrow=c(1,1))
par(mai=c(0.8,0.8,0.3,0.15))
boxplot(betas ~ parameters, ylab = "Beta", main = "SICFL", ylim=c(-4.5, 3.5),xaxt="n")
segments(1:3, lowcri, 1:3, upcri, col='grey20', lwd=2)
abline(h=0, lwd=1,lty=2 ,col='grey30')
axis(1, at=1:3, labels=names, cex=0.7) 
sig1<- lowcri*upcri>0    #mayor a cero son los dos + o dos - == efecto signif
segments((1:3)[sig1==1],lowcri[sig1==1],(1:3)[sig1==1],upcri[sig1==1],lwd=2, col='red')
#dev.off()

## Grafico la prediccion de ocupacion
###############################
# valores observados y predichos
o.ip <- seq(min(covat$IP, na.rm = TRUE),max(covat$IP, na.rm = TRUE),,400)   # con los valores umbrales
pred.ip <- scale(o.ip)
o.herb <- seq(min(cova$C_Herb, na.rm = TRUE),max(cova$C_Herb, na.rm = TRUE),,400)   # con los valores umbrales
pred.herb <- scale(o.herb)
o.darb <- seq(min(cova$D_Arb, na.rm = TRUE),max(cova$D_Arb, na.rm = TRUE),,400)   # con los valores umbrales
pred.darb <- scale(o.darb)

#tiff(file = "Figsicfla.tiff", width = 190, height = 70, units = "mm", res = 300)
par(mfrow=c(1,3), cex=0.8, cex.lab=0.7, cex.axis=0.7)
par(mar=c(1.8,1.8,0.2,0),  xpd=TRUE)
par(mgp=c(0.8,0,0))
par(oma=c(0,0,0,2))

############################################################################################
### GRAZING INTENSITY####
nsamp<-out$mcmc.info$n.samples
cri<-function(x) quantile(x,prob=c(0.05,0.95))
tmp<-out$sims.list
lpsi<-apply(tmp$lpsi[,,],1,mean)
lp<-apply(tmp$lp[,],1,mean) 

pred<- array(NA, dim=c(400,nsamp,3))
for(i in 1:nsamp){
  pred[,i,3]<-plogis(lpsi[i]+tmp$betalpsi1[i]*pred.ip)}
cri3<-apply(pred[,,3],1,cri)
plot(pred.ip, apply(pred[,,3],1,mean), ylab="Ocupacion", xlab="Intensidad pastoreo (kg LW/kg DM)",
     type="l",lty=1,lwd=3,col="dodgerblue",ylim=c(0,1),frame=F,xaxt="n" )
axis(1, at=pred.ip,labels=round(o.ip,1), cex=1, tck=0)
polygon(x= c(pred.ip, rev(pred.ip)), y= c(cri3[1,], rev(cri3[2,])), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.10), border = NA)

# leyenda
legend(x=8, y=1.00, pch=15,cex=0.6,border=NA,bty="n",
       legend= substitute(paste(italic('Intensidad pastoreo')))
       ,col=c("dodgerblue"), xpd=NA)
legend(x=8, y=0.95, pch=15,cex=0.6,border=NA,bty="n",
       legend=substitute(paste(italic('Cobertura herbacea'))),
       ,col=c("mediumpurple4"), xpd=NA)
legend(x=8, y=0.9, pch=15,cex=0.6,border=NA,bty="n",
       legend= substitute(paste(italic('Densidad arboles')))
       ,col=c("plum"), xpd=NA)

############################################################################################
### COBERTURA HERBACEA####
pred<- array(NA, dim=c(400,nsamp,3))
for(i in 1:nsamp){
  pred[,i,3]<-plogis(lpsi[i]+tmp$betalpsi2[i]*pred.herb)}
cri<-function(x) quantile(x,prob=c(0.05,0.95))
cri3<-apply(pred[,,3],1,cri)
plot(pred.herb, apply(pred[,,3],1,mean), ylab="Ocupacion", xlab="Cobertura herbacea",
     type="l",lty=1,lwd=3,col="mediumpurple4",ylim=c(0,1),frame=F,xaxt="n" )
axis(1, at=pred.herb,labels=round(o.herb,1), cex=1, tck=0)
polygon(x= c(pred.herb, rev(pred.herb)), y= c(cri3[1,], rev(cri3[2,])), 
        col =  adjustcolor("mediumpurple4", alpha.f = 0.10), border = NA)

#########################################################
### DENSIDAD ARBOLES####
pred<- array(NA, dim=c(400,nsamp,3))
for(i in 1:nsamp){
  pred[,i,3]<-plogis(lpsi[i]+tmp$betalpsi4[i]*pred.darb)}
cri3<-apply(pred[,,3],1,cri)
plot(pred.darb, apply(pred[,,3],1,mean), ylab="Occupancy", xlab="Densidad arbol",
     type="l",lty=1,lwd=3,col="plum",ylim=c(0,1),frame=F,xaxt="n" )
axis(1, at=pred.darb,labels=round(o.darb,2), cex=1, tck=0)
polygon(x= c(pred.darb, rev(pred.darb)), y= c(cri3[1,], rev(cri3[2,])), 
        col =  adjustcolor("plum", alpha.f = 0.10), border = NA)



# ------------------------------------------------------------------------
# Tareas 
# -------------------------------------------------------------------------
# 1- Agregar arbustos como covariable de ocupacion
# 2- Agregar viento como covariable de detección (no olvidar que se midio para los dos años...)
# 3- Graficar con las covariables agregadas