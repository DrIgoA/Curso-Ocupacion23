###################################
rm(list=ls(all=TRUE))
library(jagsUI)

setwd('C:\\Users\\andrea\\Documents\\PROYECTOS y COLABORACIONES\\FONTAGRO\\Analisis')

data <- read.csv("FONTAGROAvesBD.csv",sep ="," , header = T)
str(data)
head(data)

#cargo los paquetes que voy a usar
library(reshape)
library(plyr)

# selecciono una especie
SICFL <- data$SICFL

# colectar informacion de cada punto
data.info<-data[1:720,1:8]  
head(data.info)

#armar la planilla con los registros 
data.SICFL<-data.frame(data.info, SICFL)
str(data.SICFL)
head(data.SICFL)

data.melt=melt(data.SICFL, id.var=c("Predio", "Anio", "NPunto", "Rep"), measure.var="SICFL")
y=cast(data.melt, Predio ~ Anio ~ NPunto ~ Rep)

# Verifico las dimensiones
dim(y)

# Si trabajo con occupancy paso todos los numero a 0 o 1
y<-aaply(y,1,function(x) {x[x>1]<-1;x})

###################################################################################################
###################################################################################################
##  COVARIABLES no anuales

cova <- read.csv("cova_noyearNA_V4.csv", sep= ",", header = T)
head(cova)
str(cova)

cova$herbz<-scale(cova$C_Herb)
cova$darbz<-scale(cova$D_Arb)
cova$darbuz<-scale(cova$D_Arbust)
cova$dapz<-scale(cova$DAP)

covs<-cbind(cova$herbz, cova$darbz, cova$darbuz, cova$dapz)
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

# DAP
dapz.melt=melt(cova, id.var=c("Predio", "NPunto"), measure.var="dapz")
dap=cast(dapz.melt, Predio ~ NPunto)
dap<- as.matrix(dap)


########################
##  COVARIABLES anuales

covat <- read.csv("cova_year_V4.csv", sep= ",", header = T)
head(covat)
str(covat)

covat$ipz<-scale(covat$IP)

# Intensidad de pastoreo
ip.melt=melt(covat, id.var=c("Predio", "NPunto", "Anio"), measure.var="ipz")
ip=cast(ip.melt, fun.aggregate = mean, Predio ~ Anio ~ NPunto)

# Observador x anio pero mismo x rep
obs.melt=melt(covat, id.var=c("Predio", "NPunto", "Anio"), measure.var="Observador")
obs=cast(obs.melt, fun.aggregate = mean, Predio ~ Anio ~ NPunto)
dim(obs)

# al observador 2 le ponemos un "0" para el modelo
obs[obs == 2] <- 0

# Momento del dia
hora.melt=melt(covat, id.var=c("Predio", "NPunto", "Anio", "Rep"), measure.var="Tiempo_dia")
hora=cast(hora.melt, fun.aggregate = mean, Predio ~ Anio ~ NPunto ~ Rep)
dim(hora)

# a la hora del dia 2 le ponemos un "0" para el modelo
hora[hora == 2] <- 0


###################################################################################################
###################################################################################################
###################################################################################################

#cambio el orden de las dimensiones
#y1<-aperm(y,c(4,3,1,2))
#y<-y1

dim(y)

str(win.data <- list(y = y, 
                     cov1 = ip,                                      
                     cov2 = herb,
                     cov3 = darbu,                 
                     cov4 = darb, 
                     cov5 = dap,
                     obs = obs,
                     hora = hora,
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
    cov3  [e,i] ~ dnorm(mu.cov3,   tau.cov3)
    cov4  [e,i] ~ dnorm(mu.cov4,   tau.cov4)
    cov5  [e,i] ~ dnorm(mu.cov5,   tau.cov5)
    }
    for (t in 1:nyear){                                         
    for (i in 1:nsite) { 
    cov1  [e,t,i] ~ dnorm(mu.cov1,   tau.cov1)
    obs   [e,t,i] ~ dbern(pobs)
    for (k in 1:nrep){
    hora   [e,t,i,k] ~ dbern(phora)
    }
    }
    }
    }
    tau.cov1 <-pow(sd.cov1,-2)
    sd.cov1 ~ dunif(0,10)
    mu.cov1 ~ dnorm(0,0.01)
    
    tau.cov2 <-pow(sd.cov2,-2)
    sd.cov2 ~ dunif(0,10)
    mu.cov2 ~ dnorm(0,0.01)
    
    tau.cov3 <-pow(sd.cov3,-2)
    sd.cov3 ~ dunif(0,10)
    mu.cov3 ~ dnorm(0,0.01)
    
    tau.cov4 <-pow(sd.cov4,-2)
    sd.cov4 ~ dunif(0,10)
    mu.cov4 ~ dnorm(0,0.01)
    
    tau.cov5 <-pow(sd.cov5,-2)
    sd.cov5 ~ dunif(0,10)
    mu.cov5 ~ dnorm(0,0.01)
    
    pobs ~ dunif(0,1)
    phora ~ dunif(0,1)
    
    #Hyperprevias
    mu.lp           ~ dnorm(0, 1/2.25^2)
    mu.lpsi      ~ dnorm(0, 1/2.25^2)
    betap1 ~ dnorm(0, 1/2.25^2)
    betap2 ~ dnorm(0, 1/2.25^2)
    betalpsi1 ~ dnorm (0, 1/2.25^2)    
    betalpsi2 ~ dnorm (0, 1/2.25^2)    
    betalpsi3 ~ dnorm (0, 1/2.25^2)
    betalpsi4 ~ dnorm (0, 1/2.25^2)    
    betalpsi5 ~ dnorm (0, 1/2.25^2)    
    
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
    
    logit(psi[e,t,i]) <- lpsi[e,t] + betalpsi1 * cov1 [e,t,i] + betalpsi2 * cov2 [e,i] + betalpsi3 *cov3 [e,i]
    + betalpsi3 *cov3 [e,i] + betalpsi4 *cov4 [e,i] + betalpsi5 *cov5 [e,i] 
    
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
    logit(p[e,t,i,k]) <- lp[t] + betap1 * obs[e,t,i] + betap2 * hora[e,t,i,k]
    mu.p[e,t,i,k] <- z[e,t,i] * p[e,t,i,k]  
    y[e,t,i,k] ~ dbern(mu.p[e,t,i,k])        #line 69
    
# GOF assessment
    Tobs0[e,t,i,k] <- (sqrt(y[e,t,i,k]) - sqrt(p[e,t,i,k]*z[e,t,i]))^2  # FT discrepancy for observed data
    ySim[e,t,i,k] ~ dbern(p[e,t,i,k]*zSim[e,t,i])
    Tsim0[e,t,i,k] <- (sqrt(ySim[e,t,i,k]) - sqrt(p[e,t,i,k]*zSim[e,t,i]))^2  # ...and for simulated data

    }
    }
    }
    } 
#119
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

params1 <- c("lpsi","lp", "betap1","betap2","betalpsi1","betalpsi2","betalpsi3","betalpsi4",
             "betalpsi5" , "N", "Tobs", "Tsim")

ni <- 70000      
nt <- 10  
nb <- 10000 
na<- 20000     
nc <- 3

out = jags(win.data, inits, params1, "singlesp_bal_nested_time_det.txt", n.chains=nc, 
           n.iter=ni, n.burnin=nb, n.thin=nt, n.adapt = na)

pp.check(out, "Tobs", "Tsim", main="Freeman-Tukey discrepancy")
#plot(out)

save(out, file='SICFL7_global.rda')

load('SICFL7_global.rda')
print(out, dig=3)

# guardar resultados
#write.csv(out$summary, "SICFL3_summary.csv")

plot(out)
#pp.check(out, "Tobs", "Tsim") #, main="Freeman-Tukey discrepancy")



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
         mean(out$sims.list$betalpsi3),mean(out$sims.list$betalpsi4),
         mean(out$sims.list$betalpsi5))
cri<-c(cri(out$sims.list$betalpsi1),cri(out$sims.list$betalpsi2),
       cri(out$sims.list$betalpsi3),cri(out$sims.list$betalpsi4),
       cri(out$sims.list$betalpsi5))
lowcri<-c(cri[1],cri[3],cri[5],cri[7],cri[9])
upcri<-c(cri[2],cri[4],cri[6],cri[8],cri[10])


# compile parameters & names
parameters=c("betalpsi1","betalpsi2","betalpsi3","betalpsi4","betalpsi5")
names= c('ip', 'herb', 'darbu', 'darb','dap')

# plot Beta coefficients
#tiff(file = "ammhu3_coeff.tiff", width = 200, height = 200, units = "mm", res = 300)
par(mfrow=c(1,1))
par(mai=c(0.8,0.8,0.3,0.15))
boxplot(betas ~ parameters, ylab = "Beta", main = "SICFL", ylim=c(-4.5, 3.5),xaxt="n")
segments(1:5, lowcri, 1:5, upcri, col='grey20')
abline(h=0, lwd=1,lty=2 ,col='grey30')
axis(1, at=1:5, labels=names, cex=0.7) 
sig1<- lowcri*upcri>0    #mayor a cero son los dos + o dos - == efecto signif
segments((1:5)[sig1==1],lowcri[sig1==1],(1:5)[sig1==1],upcri[sig1==1],lwd=1, col='red')
#dev.off()
