##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########           Ejemplo Muestreo           #########
###########                                      #########
##########################################################
########              Ejemplo de:                   ######
########          Introduction to WinBUGS           ######
########             for ecologists                 ######
########    A bayesian approach to regression,      ######
########  ANOVA, mixed models an related analyses   ######
########            Marc Kery - 2010                ######
##########################################################

# Ejercicio Halcones peregrinos
# ----------------------------------------------------------------------------
# En este ejercicio trabajaremos con la masa de pelegrinos machos.
# Los pelegrinos del Oeste Europeo pesan en promedio 600g y
# Monneret (2006) establece como rango de peso entre 500-680g.
# Asumiento una distribucion Nomral para la masa, implica un desvio
# estandar de alrededor de 30g.

# Genero los datos de peso para 1000 machos
set.seed(1234)
y1000 <- rnorm(n = 1000, mean = 600, sd = 30)
                            
# Plot data
hist(y1000, col = 'grey', xlim = c(450,750), main = ' Body mass (g) of
1000 male peregrines')

# 1)  Implemente un modelo lineal con estadistica frecuentista (utilizando 
#     la funcion lm) para estimar la la masa media de los halcones machos.

# SOLUCION
summary(lm(y1000 ~ 1))

# 2) Implemente el modelo con estadistica bayesiana, considerando
#    una previa no informativa con una distribcion uniforme = (0,5000)

# SOLUCION
library(jagsUI) 

# Data
data <- list(mass = y1000, nobs = length(y1000))

#Modelo con previas no informativas
sink("peregrinos-no-informativa.jags")
cat("

model
{
  # Previas No INFORMATIVAS
    mean ~ dunif (0, 5000)      
    varianza ~ dlnorm (0 ,1.0E-6)  

    prec <- 1/varianza             

  #Likelihood
  for (i in 1:nobs) {                # para cada uno de los individuos
     mass[i] ~ dnorm (mean, prec)   
 
   }
}

",fill = TRUE)
sink()

# Indicamos los valores iniciales 
inits <- function() list(varianza = 100, mean = rnorm(1,600))

# MCMC settings
ni <- 10000     # numero de iteraciones
nt <- 2         # tasa de thining
nb <- 1000      # iteraciones para el burn in
nc <- 3         # numero de cadenas que corremos

# Parametros que se van monitorear 
parameters <- c("mean", "varianza")

# Implementamos el modelo con jags
out <- jags(data, inits, parameters, "peregrinos-no-informativa.jags",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

out 

# ------------------------------------------------------------------------
# Preguntas
# -------------------------------------------------------------------------
# 1. Disgnostique visualmente si las cadenas convergen
# 
# 2. Que relacion observa entre la media estimada con el modelo lineal (lm) 
#     utilizando estadistica frecuentista y bayesiana?
# 
# 3. Grafica de manera conjunta las distribuciones previas y posteriores 
#   del modelo
# 
# 4. Que sucede si indicamos valores iniciales por fuera de la distribucion
#  previa indicada para el parametro?
#  
# 5. En el ejercicio utilizaron previas no informativas, lo que implica que su 
#     efecto en la distribucion posterior sea minimo. Podemos incorporar previas
#     informativas con valores que vayan entre 500 y 590g. Incorpora estos 
#     nuevos valores de distribuciones previas e implementa el nuevo modelo.
#     Como cambian los parametros y su incertidumbre? Esto es positivo o negativo?
#     
#     Compara graficamente las distribuciones previa y poseterior.
#     

#Modelo con previas informativas
sink("peregrinos-informativa.jags")
cat("

model
{
  # Previas INFORMATIVAS
    mean ~ dunif (500, 590)      
    varianza ~ dlnorm (0 ,1.0E-6)  

    prec <- 1/varianza             

  #Likelihood
  for (i in 1:nobs) {                # para cada uno de los individuos
     mass[i] ~ dnorm (mean, prec)   
 
   }
}

",fill = TRUE)
sink()

# Indicamos los valores iniciales 
inits <- function() list(varianza = 100, mean = rnorm(1,600))

# MCMC settings
ni <- 10000     # numero de iteraciones
nt <- 2         # tasa de thining
nb <- 1000      # iteraciones para el burn in
nc <- 3         # numero de cadenas que corremos

# Parametros que se van monitorear 
parameters <- c("mean", "varianza")

# Implementamos el modelo con jags
out <- jags(data, inits, parameters, "peregrinos-informativa.jags",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

out 

