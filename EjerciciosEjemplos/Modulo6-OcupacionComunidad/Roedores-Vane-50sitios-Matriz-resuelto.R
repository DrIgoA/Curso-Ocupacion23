##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########         Ejercicio Ocupacion          #########
###########          Roedores - Matriz           #########
##########################################################


# Los datos corresponden a dos muestreos (season) de roedores realizados en la provincia 
# de Cordoba. El diseño posee 50 sitios, en los cuales se coloco 1 linea con 
# 20 trampas de captura viva tipo Shermann.  
# Se capturaron 9 especies de roedores: Akodon azarae, A. dolores, Calomys 
# laucha, C. musculinus, C. venustus, Oligorizomys flavesces, Oxymicterus rufus, 
# Tylamys pallidor, Monodelphys dimidiata.
# 

# Acá borramos la memoria de R, asi nos aseguramos que no haya objetos cargados
# con anterioridad
rm(list=ls(all=TRUE)) 

library(jagsUI)    #paquete JAGS

data<-read.csv("datos_roedores.csv",header = T)
str(data)
head(data)


# Para poder ordenar los datos en la matriz es importante que esten todos los sitios 
# por mas que no hayan capturado nada y todas las potenciales sp por sitios
# Las matrices tienen el mismo numero de filas y columnas

# Para utilizar el siguiente script, los datos deben estar ordenados primero 
# por sesion y despues por especie
data.season.1 <-data[1:450,1:7]      # datos del primer muestreo
data.season.2 <-data[451:900,1:7]    # datos del segundo muestreo

# Las listas con los nombres nos van a permitir organizar luego el array
species.list <- list("Aa","Ad","Cl","Cm","Cv","Md","Of","Or","Tp")   # alphabetic list
season.list<- list("1", "2")
site.list<- list("1", "2","3","4","5","6","7","8","9","10",
                 "11","12","13","14","15","16","17","18","19","20",
                 "21","22","23","24","25","26","27","28","29","30",
                 "31","32","33","34","35","36","37","38","39","40",
                 "41","42","43","44","45","46","47","48","49","50")
rep.list<- list("1", "2","3","4")

#seleccionamos solo las columnas que corresponden a los datos de conteos 
# Season 1 
counts.season.1 <- cbind(data.season.1$X1, data.season.1$X2, data.season.1$X3,data.season.1$X4)
DET.season.1 <- counts.season.1
DET.season.1[DET.season.1 > 1] <- 1       # transformamos los valores en datos deteccion/no deteccion

# Season 2 
counts.season.2 <- cbind(data.season.2$X1, data.season.2$X2, data.season.2$X3,data.season.2$X4)  
DET.season.2 <- counts.season.2
DET.season.2[DET.season.2 > 1] <- 1       # transformamos los valores en datos deteccion/no deteccion

nsite <- 50                     # numero de sitios
nrep <- 4                       # number of replicate surveys per season  - cada una de las 4 noches
nspec <- length(species.list)   # 9 species en data.1
nseason<-2                      # muestreos

# Creo el array vacio para luego completarlo con los datos del muestreo
Y.season.1 <- array(NA, dim = c(nsite, nrep, nspec))

# Nombre a cada una de las dimensiones del array
dimnames(Y.season.1) <- list(site=site.list,rep= rep.list, especie=species.list)
dim(Y.season.1)
Y.season.1

# Para completar el array hago un loop para cada conjunto de datos
# Lo que hace es completar cada array de especies con las filas que corresponden
# Esto es posible por como estan ordenados los datos, repitiendo todos los sitios
# por cada especie
#  i corresponden a las especies. Pruebo con i=1 para entender como funciona

for(i in 1:nspec){
  Y.season.1[,,i] <- DET.season.1[((i-1)*nsite+1):(i*nsite),]     #i es de la especie asi te rellena para cada sp con 50 sitios
}

Y.season.2 <- array(NA, dim = c(nsite, nrep, nspec))
dimnames(Y.season.2) <- list(site=site.list,rep= rep.list, especie=species.list)
dim(Y.season.2)
Y.season.2

for(i in 1:nspec){
  Y.season.2[,,i] <- DET.season.2[((i-1)*nsite+1):(i*nsite),]     #i es de la especie asi te rellena para cada sp con 50 sitios
}

Y.season.2

library(abind) #necesaria para combinar las matrices

y<-abind(Y.season.1,Y.season.2,along = 4)  # es 4 porque es en una dimension mas de las que ya tienen (filas, columnas y especies)
dimnames(y) <- list(site=site.list,rep= rep.list, especie=species.list,season=season.list)  #nombres de las dimensiones
y

# Escalo la covariable
is <- (scale(data$is))

### Construyo matrices de las covariables
is<- array(is, dim = c(nsite))
dimnames(is) <- list(site=site.list)  #nombres de las dimensiones

J1<- data$J[1:nsite]
J2<- data$J[451:500]

J<-abind(matrix(J1),matrix(J2))

# ----------------
# PARTE 1
# ---------------
# creo una lista con los datos que utilizara el modelo
  str(roed.data <- list(y=y, nsite=dim(y)[1],nspec=dim(y)[3],nrep=dim(y)[2],
                      nseason=dim(y)[4], is = is))

## JAGS code
sink("mod_is.jags")  
cat("
    model{
    
    # Hyperpriors 
      mu.lpsi ~ dnorm(0, 0.1)
      mu.b1 ~   dnorm(0, 0.1)
      mu.lp ~   dnorm(0, 0.1)
    
      tau.lpsi <- pow(sd.lpsi,-2)   
      tau.lp <- pow(sd.lp,-2)   
      tau.b1 <- pow(sd.b1,-2)

      sd.lpsi ~ dunif(0,1)
      sd.lp ~ dunif(0,1)
      sd.b1 ~ dunif(0,1)
    
    # priors for species-specific effects
      for (k in 1:nspec){
        lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
        b1[k] ~ dnorm (mu.b1, tau.b1)
    
      for (j in 1:nrep) {  
        lp[j,k] ~ dnorm (mu.lp, tau.lp)    
       } #nrep
  
  # Ecological model, process model (true occurrence at site i) 
    for (i in 1:nsite) {                                           #loop sobre sitios 
      for (l in 1:nseason){                                         #loop sobre tiempo 
          logit(psi[i,k,l]) <- lpsi[k] + b1[k]*is[i]
    
          z[i,k,l] ~ dbern(psi[i,k,l])
      } #season
    } #site
    
    # Observation model for site i, replicate nrep=j, nspec=k, season=l
    for (l in 1:nseason){                                         #loop sobre tiempo
      for (j in 1:nrep) {  
         logit(p[j,k,l]) <- lp[j,k]     #poner por fuera de i
       } #season
     }

    for (i in 1:nsite) {                                           #loop sobre sitios 
      for (j in 1:nrep) {                                           #loop sobre noche 
        for (l in 1:nseason){                                         #loop sobre tiempo 
          mu.p[i,j,k,l] <- p[j,k,l]*z[i,k,l]
          y[i,j,k,l] ~ dbern(mu.p[i,j,k,l])
        } #nseason
      } #nrep
    } #nsite
    } #nspec
  
      #Derived quantities
      for (i in 1:nsite) {
        for(l in 1:nseason){
          Nsite [i,l] <- sum (z[i,,l])
    }
    }
    
      for (k in 1:nspec){
         for(l in 1:nseason){
           Nocc.fs [k,l] <- sum (z[,k,l])
    }   
    }    
  
    
}
",fill = TRUE)
sink()
zst <- array(1,c(nsite, nspec,2))

# valores iniciales
inits <- function() list(z=zst,lpsi=rnorm(roed.data$nspec), b1=rnorm(roed.data$nspec),
                         mu.lpsi=rnorm(1),mu.b1=rnorm(1), mu.lp=rnorm(1))

inits <- function() list(z=zst,lpsi=runif(9), b1=runif(9),
                         mu.lpsi=runif(1),mu.b1=runif(1))

params1 <- c("lpsi", "b1", "lp","mu.lpsi","mu.lp","mu.b1")
params2 <- c("Nsite", "Nocc.fs")

# ajustes de MCMC
ni <- 10000
nt <- 10
nb <- 1000
nc <- 3

# llamar JAGS desde R 
out1 = jags(roed.data,inits, params1, "mod_is.jags", n.chains=nc, 
                  n.iter=ni, n.burnin=nb, n.thin=nt) # con PARALLEL = TRUE esto corren cada cadena en cada nucleo entonces hace
out2 = jags(roed.data,inits, params2, "mod_is.jags", n.chains=nc, 
                  n.iter=ni, n.burnin=nb, n.thin=nt) #con PARALLEL =TRUE esto corren cada cadena en cada nucleo entonces hace

# Guardar los objetos 
# save(out1, file='out1-roed-matrix.rda') 
# save(out2, file='out2-roed-matrix.rda') 

### Actividad
### 1. Seleccione los parametros que desea monitorear para evaluar como el indice
###    de Shannon afecta a cada una de las especies y a la comunidad de roedores. 
###    Corra el modelo y observe si convergen las cadenas
### 2. Incorpore al modelo anterior las covariables natural y cultivo, 
###    las cuales corresponden al tipo de ambiente en el cual se encontraba ubicada
###    la línea de trampas de captura viva.
### 3. Grafique el efecto del Indice de Shannon sobre la ocupacion de cada especie

library(denstrip)
plot(out1$sims.list$b1, xlim=c(-3, 3), ylim=c(1, 9), xlab="", ylab="", type="n", axes =
       F, main = "Density strip plots")
axis(1)
axis(2, at = 1:9, labels = c('Aa','Ad','Cm','Cl','Cv','Of', 'Or','Md', 'Tp'), las = 1)
abline(v = c(--1,-0.5,0.5,1), col = "grey") ; abline(v = 0)
for(k in 1:9){
  denstrip(unlist(out1$sims.list$b1[,k]), at = k,
           colmax = "#4292c6", colmin = "#f7fbff")
}



# ----------------
# PARTE 2
# ---------------

# Miodel Regression para Riqueza e IS
is <- (scale(data$is[1:50]))

#promedio de riqueza y intervalos para todas las escalas
N.pm <-round(apply(out2$mean$Nsite,1,mean))
N.psd <- apply(out2$sd$Nsite,1,mean)
N.cri <- apply(out2$sims.list$Nsite, 2,function(x) quantile(x,prob=c(0.125,0.875)))

is.pred<- seq(-1.436,2.6,length=500)

str(win.data<- list(is=is[1:50],N=N.pm, psd=N.psd,
                    n= length(N.pm), pred.is=is.pred, npred=length(is.pred)))

sink("meta.analysis.150_revision.txt")
cat("

    model{
    
    #Priors
    for (v in 1:3) {     #previas for intercep and plonomial coefficients
    beta[v] ~dnorm (0,0.0001)
    }
    tau.site<-pow(sd.site, -2)
    sd.site~dunif(0,10)
    
    # Likelihood
    for (i in 1:n) {
    N[i]~dnorm(muN[i], tau.psd[i])
    tau.psd[i]<-pow(psd[i],-2)
    muN[i]<- beta[1]+ beta[2]*is[i] +beta[3]*pow(is[i],2) + eps.site[i]
    eps.site [i] ~ dnorm (0, tau.site)
    }
    
    #get predictions for plot
    for (i in 1:npred){
    Npred[i] <- beta[1] + beta[2]*pred.is[i] + beta[3]*pow(pred.is[i],2) 
    
    }
    }
", fill= TRUE)
sink()

inits<- function() list(beta=rnorm(3))
params<- c("beta", "sd.site","Npred")
ni<- 12000 ; nt<- 10 ; nb <- 2000 ; nc<- 3

outR<- jags(win.data, params,inits = inits, "meta.analysis.150_revision.txt", n.chains = nc,
           n.thin = nt, n.iter = ni, n.burnin = nb)

# Guardar una figura en formato tiff
tiff(file = "fig 5b_nuevo_estimada.tiff",                #Guardados como Tiff
     width = 84, height = 70,
     units = "mm", res =1200)

par(mai=c(0.35,0.4,0.1,0.3),xpd=F)
#plot estimares as a function of is
plot(is, N.pm,  ylab= "Estimated species richness", 
     ylim = c(0,6.5), xlim=c(-1.8,2.7),pch= 8, bty="l",cex= 0.8,cex.lab=0.9,
     xaxt="n",yaxt="n", col="grey10")
axis(2,at=0:6,ylab= "Estimated species richness",labels = c(seq(0,6.5,1)), las=1,
     cex.axis=0.65,tck=-0.02, mgp=c(3,0.4,0))
axis(1,at=-2:2.8,labels = c(seq(-2,2.8)), las=1,
     cex.axis=0.6,tck=-0.02, mgp=c(3,-0.1,0))
mtext(1,line=0.55,text="Shannon Habitat Index",cex = 0.65)
mtext(2,line=1,text="Estimated species richness",cex = 0.65)

lines(is.pred, outR$mean$Npred, col="grey10", bty="l", lwd=1.2)

dev.off()


