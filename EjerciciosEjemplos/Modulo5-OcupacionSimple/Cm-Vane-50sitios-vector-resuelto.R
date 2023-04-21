##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########   Ejercicio Ocupacion - vectores     #########
###########        Una especie - Una estacion    #########
##########################################################


# Los datos corresponden a las capturas de Calomys musculinus en un muestreo
# de 50 sitios, en los cuales se coloco 1 linea con 20 trampas de captura viva 
# tipo Shermann.  

# Acá borramos la memoria de R, asi nos aseguramos que no haya objetos cargados
# con anterioridad
rm(list=ls(all=TRUE)) 

library(jagsUI)    #paquete JAGS

data<-read.csv("datos_cm.csv",header = T)
str(data)
head(data)
 attach(data)

# creo una lista con los datos que utilizara el modelo
  str(roed.data <- list(y=y, nsite=50, is = is, J=J))

## JAGS code
sink("mod_cm.jags")  
cat("
    model{
    
    # Priors 
      lpsi ~ dnorm(0, 0.1)
      b1 ~   dnorm(0, 0.1)
      p ~   dunif(0, 1)
    
    # Ecological model, process model (true occurrence at site i) 
    for (i in 1:nsite) {                                           #loop sobre sitios 
          logit(psi[i]) <- lpsi + b1*is[i]
    
          z[i] ~ dbin(psi[i],1)

          tmp[i] <- p * z[i]
          y[i] ~ dbin(tmp[i], J[i])
     } #nsite
    
    }
    
  
",fill = TRUE)
sink()

# valores iniciales
inits <- function() list(lpsi=runif(1), b1=runif(1))
                         

params1 <- c("lpsi", "b1", "p")

# ajustes de MCMC
ni <- 10000
nt <- 10
nb <- 1000
nc <- 3

# llamar JAGS desde R 
out1 = jags(roed.data,inits, params1, "mod_cm.jags", n.chains=nc, 
                  n.iter=ni, n.burnin=nb, n.thin=nt) # con PARALLEL = TRUE esto corren cada cadena en cada nucleo entonces hace


traceplot(out1, parameters=c("lpsi", "b1")) 
densityplot(out1, parameters=c("lpsi", "b1", "p")) 


