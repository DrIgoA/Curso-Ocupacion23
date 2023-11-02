##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########         Rodrigo J. Alonso            #########
###########                                      #########
##########################################################

#Los datos son de muestreos de pequeños mamíferos, llevados a cabo de forma estacional durante 2019 y 2021. Se muestrearon, en total, 8 establecimientos dedicados a la producción de carne a corral (feedlots) y 8 a la producción de leche (tambos), totalizando 16 establecimientos. En cada año se muestrearon 4 establecimientos de cada tipo. Dentro de cada establecimiento, se dispusieron trampas de captura viva, tipo sherman y tipo jaula, a lo largo de transectas. Las transectas se ubicaron en distintos ambientes dentro de los establecimientos, por ejemplo zanjas, comederos, galpones de alimentos, etc. En total, se capturaron 7 especies de pequeños mamíferos (sp_1 a sp_7). 

#Limpio el environment
rm(list=ls()) 

#Cargo el paquete JAGS, abind, jagshelper y beepr
library(jagsUI)    
library(abind)
library(jagshelper)
library(beepr)

#cargo el working directory
setwd("E:/Rodrigo/Publicaciones/JAM 2023/Alonso Rodrigo/evaluacion_Alonso")


#Si no quieren correr el script, pueden cargar el environment:

#Cargo mis datos
data_ok<-read.csv("datos_propios.csv",header = T, sep =";")

#Armo las matrices. Empiezo dividiendo la información por año:
data1 <- subset(data_ok, anio == 2019)
names(data1)
data1 <- data1[,1:17]
data2 <- subset(data_ok, anio == 2021)
data2 <- data2[,1:17]

#Para ordenar el array por especies, hago listas con las diferentes dimensiones de las matrices:

#Lista de especies
levels(as.factor(data1$especie))
species.list <- list("sp_1","sp_2","sp_3","sp_4","sp_5","sp_6","sp_7")

#Lista de años
year.list<- list("1", "2") # 

#Lista de estacion del año
season.list <- list("1","2","3","4")

#Lista de sitios
site.list<- list("1", "2","3","4","5","6","7","8","9","10",
                 "11","12","13","14","15","16","17","18","19","20",
                 "21","22","23","24","25")

#Lista de repeticiones
rep.list<- list("1", "2","3")

#Tomo los datos de conteos y los transformo en detección/no detección, para cada año:
#año 2019
counts.year.1 <- cbind(data1$count1, data1$count2, data1$count3)
DET.year.1 <- counts.year.1
DET.year.1[DET.year.1 > 1] <- 1 

# año 2021 
counts.year.2 <- cbind(data2$count1, data2$count2, data2$count3)
DET.year.2 <- counts.year.2
DET.year.2[DET.year.2 > 1] <- 1 

#Para armar los array de matrices voy a usar un loop. Tengo que especificar algunos límites de los loops:
nsite <- 25 # numero de sitios
nrep <- 3 # número de repeticiones
nspec <- length(species.list)   # 7 especies
nyears <- 2 # años
nseasons <- 4

#Creo el array vacio
Y.year.1 <- array(NA, dim = c(nsite, nrep, nspec, nseasons))
Y.year.2 <- array(NA, dim = c(nsite, nrep, nspec, nseasons))

# Le doy un nombre a cada una de las dimensiones del array
#2019
dimnames(Y.year.1) <- list(site=site.list,rep= rep.list, especie=species.list, season=season.list)
dim(Y.year.1)
Y.year.1 #El array está vacio. Lo tego que llenar con los datos.

#For para completar el array. El sub indice indica qué estoy llenando y responde al orden que le di en el array. En este caso, las especies estan en el 4to lugar.

for(i in 1:nspec){
  Y.year.1[,,i,] <- DET.year.1[((i-1)*nsite+1):(i*nsite),]     #i es de la especie asi te rellena para cada sp con 25 sitios
}
Y.year.1 #El array del año 1 está completo

#2021
dimnames(Y.year.2) <- list(site=site.list,rep= rep.list, especie=species.list, season=season.list)
dim(Y.year.2)
Y.year.2 #El array está vacio. Lo tego que llenar con los datos.


for(i in 1:nspec){
  Y.year.2[,,i,] <- DET.year.2[((i-1)*nsite+1):(i*nsite),]     #i es de la especie asi te rellena para cada sp con 25 sitios
}
Y.year.2 #El array del año 2 está completo

#Una vez que tengo las dos matrices de los años, tengo que unir los rearreglos. Para eso uso el comando "abind"
y<-abind(Y.year.1,Y.year.2,along = 5)  # es 5 porque es en una dimension mas de las que ya tienen (filas(sitios), columnas(rep),  especies, estacion), le agrego año.

dimnames(y) <- list(site=site.list,rep= rep.list, especie=species.list, season=season.list, year=year.list)  #nombres de las dimensiones

#La matriz de "variable respuesta" (detecciones):
y


#### Variables explicativas.
#Empiezo por las coavariables más fáciles para aprender a armar la matriz correctamente:

#Variables de ocupación
#tipo de establecimiento, 1: tambo, 0: feedlot
tipo <- data_ok$tipo_t
tipo<- array(tipo, dim = c(nsite,nyears))
dimnames(tipo) <- list(site=site.list, year=year.list)  #nombres de las dimensiones

tipo #Generé una matriz con el tipo de establecimiento para los 25 sitios de cada año. Esta variable, entonces, depende del sitio y del año.

#Estacion del año
summer <- data_ok$summer
summer<- array(summer, dim = c(nsite, nyears, nseasons))
dimnames(summer) <- list(site=site.list, year=year.list, season=season.list)
dim(summer)

fall <- data_ok$fall
fall<- array(fall, dim = c(nsite, nyears, nseasons))
dimnames(fall) <- list(site=site.list, year=year.list, season=season.list)

winter <- data_ok$winter
winter<- array(winter, dim = c(nsite, nyears, nseasons))
dimnames(fall) <- list(site=site.list, year=year.list, season=season.list)

spring <- data_ok$spring
spring<- array(spring, dim = c(nsite, nyears, nseasons))
dimnames(spring) <- list(site=site.list, year=year.list, season=season.list)

#Ambientes/habitats
zanja <- data_ok$zanja
zanja<- array(zanja, dim = c(nsite, nyears))
dimnames(zanja) <- list(site=site.list, year=year.list)

galpon <- data_ok$galpon
galpon <- array(galpon, dim = c(nsite, nyears))
dimnames(galpon) <- list(site=site.list, year=year.list)

#es el ambiente "tambo" no el tipo de sistema. Es donde se ordeñan las vacas.
tambo <- data_ok$tambo
tambo <- array(tambo, dim = c(nsite, nyears))
dimnames(tambo) <- list(site=site.list, year=year.list)

vegetacion <- data_ok$vegetacion
vegetacion <- array(vegetacion, dim = c(nsite, nyears))
dimnames(vegetacion) <- list(site=site.list, year=year.list)

comedero <- data_ok$comedero
comedero <- array(comedero, dim = c(nsite, nyears))
dimnames(comedero) <- list(site=site.list, year=year.list)

silo <- data_ok$silo
silo <- array(silo, dim = c(nsite, nyears))
dimnames(silo) <- list(site=site.list, year=year.list)

foza <- data_ok$foza
foza <- array(foza, dim = c(nsite, nyears))
dimnames(foza) <- list(site=site.list, year=year.list)

#Ahora, pruebo con las variables de detección. Tengo que pensar cómo meter la temperatura (diaria) y las noches (diaria).

#Saco las temperaturas por año y estacion, las voy a unir luego

#verano 2019
data1_bis <- subset(data_ok, anio == 2019)
data1_bis <- subset(data1_bis , especie == "sp_1")
data1_bis <- subset(data1_bis, estacion == "verano")

temp1_1 <- data1_bis$temp_media_y1
#temp1_1 <- scale(temp1_1)

temp2_1 <- data1_bis$temp_media_y2
#temp2_1 <- scale(temp2_1)

temp3_1 <- data1_bis$temp_media_y3
#temp3_1 <- scale(temp3_1)

tempx_1<- abind(temp1_1, temp2_1, temp3_1)
temp_2019v<- array(tempx_1, dim = c(nsite, nrep))

#Otoño 2019
data1_bis <- subset(data_ok, anio == 2019)
data1_bis <- subset(data1_bis , especie == "sp_1")
data2_bis <- subset(data1_bis, estacion == "otonio")

temp1_o <- data2_bis$temp_media_y1
#temp1_1 <- scale(temp1_1)

temp2_o <- data2_bis$temp_media_y2
#temp2_1 <- scale(temp2_1)

temp3_o <- data2_bis$temp_media_y3
#temp3_1 <- scale(temp3_1)

temp_o_2019<- abind(temp1_o, temp2_o, temp3_o)
temp_2019o<- array(temp_o_2019, dim = c(nsite, nrep))


#Invierno 2019
data1_bis <- subset(data_ok, anio == 2019)
data1_bis <- subset(data1_bis , especie == "sp_1")
data2_bis <- subset(data1_bis, estacion == "invierno")

temp1_i <- data2_bis$temp_media_y1
#temp1_1 <- scale(temp1_1)

temp2_i <- data2_bis$temp_media_y2
#temp2_1 <- scale(temp2_1)

temp3_i <- data2_bis$temp_media_y3
#temp3_1 <- scale(temp3_1)

temp_i_2019<- abind(temp1_i, temp2_i, temp3_i)
temp_2019i<- array(temp_i_2019, dim = c(nsite, nrep))

#Primavera 2019
data1_bis <- subset(data_ok, anio == 2019)
data1_bis <- subset(data1_bis , especie == "sp_1")
data2_bis <- subset(data1_bis, estacion == "primavera")

temp1_p <- data2_bis$temp_media_y1
#temp1_1 <- scale(temp1_1)

temp2_p <- data2_bis$temp_media_y2
#temp2_1 <- scale(temp2_1)

temp3_p <- data2_bis$temp_media_y3
#temp3_1 <- scale(temp3_1)

temp_p_2019<- abind(temp1_p, temp2_p, temp3_p)
temp_2019p<- array(temp_i_2019, dim = c(nsite, nrep))

#uno todo 2019
temp_2019<-abind(temp_2019v,temp_2019o,temp_2019i,temp_2019p,along = 3)
dimnames(temp_2019) <- list(site=site.list, rep=rep.list, season=season.list)



#verano 2021
data1_bis <- subset(data_ok, anio == 2021)
data1_bis <- subset(data1_bis , especie == "sp_1")
data1_bis <- subset(data1_bis, estacion == "verano")

temp1_1 <- data1_bis$temp_media_y1
#temp1_1 <- scale(temp1_1)

temp2_1 <- data1_bis$temp_media_y2
#temp2_1 <- scale(temp2_1)

temp3_1 <- data1_bis$temp_media_y3
#temp3_1 <- scale(temp3_1)

tempx_1<- abind(temp1_1, temp2_1, temp3_1)
temp_2021v<- array(tempx_1, dim = c(nsite, nrep))

#Otoño 2021
data1_bis <- subset(data_ok, anio == 2021)
data1_bis <- subset(data1_bis , especie == "sp_1")
data2_bis <- subset(data1_bis, estacion == "otonio")

temp1_o <- data2_bis$temp_media_y1
#temp1_1 <- scale(temp1_1)

temp2_o <- data2_bis$temp_media_y2
#temp2_1 <- scale(temp2_1)

temp3_o <- data2_bis$temp_media_y3
#temp3_1 <- scale(temp3_1)

temp_o_2021<- abind(temp1_o, temp2_o, temp3_o)
temp_2021o<- array(temp_o_2021, dim = c(nsite, nrep))


#Invierno 2021
data1_bis <- subset(data_ok, anio == 2021)
data1_bis <- subset(data1_bis , especie == "sp_1")
data2_bis <- subset(data1_bis, estacion == "invierno")

temp1_i <- data2_bis$temp_media_y1
#temp1_1 <- scale(temp1_1)

temp2_i <- data2_bis$temp_media_y2
#temp2_1 <- scale(temp2_1)

temp3_i <- data2_bis$temp_media_y3
#temp3_1 <- scale(temp3_1)

temp_i_2021<- abind(temp1_i, temp2_i, temp3_i)
temp_2021i<- array(temp_i_2021, dim = c(nsite, nrep))

#Primavera 2021
data1_bis <- subset(data_ok, anio == 2021)
data1_bis <- subset(data1_bis , especie == "sp_1")
data2_bis <- subset(data1_bis, estacion == "primavera")

temp1_p <- data2_bis$temp_media_y1
#temp1_1 <- scale(temp1_1)

temp2_p <- data2_bis$temp_media_y2
#temp2_1 <- scale(temp2_1)

temp3_p <- data2_bis$temp_media_y3
#temp3_1 <- scale(temp3_1)

temp_p_2021<- abind(temp1_p, temp2_p, temp3_p)
temp_2021p<- array(temp_i_2021, dim = c(nsite, nrep))

#uno todo 2021
temp_2021<-abind(temp_2021v,temp_2021o,temp_2021i,temp_2021p,along = 3)
dimnames(temp_2021) <- list(site=site.list, rep=rep.list, season=season.list)

#junto todo en temp
temp<-abind(temp_2019, temp_2021, along = 4)
dimnames(temp) <- list(site=site.list, rep=rep.list, season=season.list, year = year.list)

#Decidí no escalar la temperatura, porque me quedaban medio raros los valores del escalado.

#Creo a la noche como variable (para ver si el hecho de remover los bichos afecta la detectabilidad)

data1_bis <- subset(data_ok, anio == 2021)
data1_bis <- subset(data1_bis , especie == "sp_1")
data2_bis <- subset(data1_bis, estacion == "primavera")

noches <- cbind(data2_bis$noche_1,data2_bis$noche_2,data2_bis$noche_3)
dimnames(noches) <- list(site=site.list, rep=rep.list)

#Chequeo dimensiones de y
dim(y)

#Ahora voy a armar el modelo propiamente dicho:
# Junto los datos en una lista:
str(my.data <- list(y=y,
                    tipo = tipo,
                    temp = temp, 
                    summer = summer, 
                    fall = fall, 
                    winter = winter, 
                    spring = spring, 
                    noches = noches,
                    zanja = zanja,
                    galpon = galpon,
                    tambo = tambo,
                    vegetacion = vegetacion,
                    comedero = comedero,
                    silo = silo,
                    foza = foza,
                    nrep = nrep,
                    nsite = nsite, 
                    nseasons = nseasons,
                    nspec = nspec,
                    nyears = nyears))
### código de JAGS
sink("mimodelo.jags")  
cat("
model{

# Hyperpriors 
 mu.co.lpsi ~ dnorm(0, 0.01) #queda
 mu.co.betalpsi1 ~ dnorm(0, 0.01)
 mu.co.betalpsi2 ~ dnorm(0, 0.01)
 mu.co.betalpsi3 ~ dnorm(0, 0.01)
 mu.co.betalpsi4 ~ dnorm(0, 0.01)
 mu.co.betalpsi5 ~ dnorm(0, 0.01)
 mu.co.betalpsi6 ~ dnorm(0, 0.01)
 mu.co.betalpsi7 ~ dnorm(0, 0.01)
 mu.co.betalpsi8 ~ dnorm(0, 0.01)
 mu.co.betalpsi9 ~ dnorm(0, 0.01)
 mu.co.betalpsi10 ~ dnorm(0, 0.01)
 mu.co.betalpsi11 ~ dnorm(0, 0.01)
 mu.co.betalpsi12 ~ dnorm(0, 0.01)
 mu.co.betalpsi13 ~ dnorm(0, 0.01)
 mu.co.lp ~ dnorm(0, 0.01)
 mu.co.betalp1 ~ dnorm(0, 0.01)
 mu.co.betalp2 ~ dnorm(0, 0.01)
 
 tau.co.lpsi <- pow(sd.co.lpsi,-2)   
 sd.co.lpsi ~ dunif(0,8)
 tau.co.betalpsi1 <- pow(sd.co.betalpsi1,-2)
 sd.co.betalpsi1 ~ dunif(0,8)
 tau.co.betalpsi2 <- pow(sd.co.betalpsi2,-2)
 sd.co.betalpsi2 ~ dunif(0,8)
 tau.co.betalpsi3 <- pow(sd.co.betalpsi3,-2)
 sd.co.betalpsi3 ~ dunif(0,8)
 tau.co.betalpsi4 <- pow(sd.co.betalpsi4,-2)
 sd.co.betalpsi4 ~ dunif(0,8)
 tau.co.betalpsi5 <- pow(sd.co.betalpsi5,-2)
 sd.co.betalpsi5 ~ dunif(0,8)
 tau.co.betalpsi6 <- pow(sd.co.betalpsi6,-2)
 sd.co.betalpsi6 ~ dunif(0,8)
 tau.co.betalpsi7 <- pow(sd.co.betalpsi7,-2)
 sd.co.betalpsi7 ~ dunif(0,8)
 tau.co.betalpsi8 <- pow(sd.co.betalpsi8,-2)
 sd.co.betalpsi8 ~ dunif(0,8)
 tau.co.betalpsi9 <- pow(sd.co.betalpsi9,-2)
 sd.co.betalpsi9 ~ dunif(0,8)
 tau.co.betalpsi10 <- pow(sd.co.betalpsi10,-2)
 sd.co.betalpsi10 ~ dunif(0,8)
 tau.co.betalpsi11 <- pow(sd.co.betalpsi11,-2)
 sd.co.betalpsi11 ~ dunif(0,8)
 tau.co.betalpsi12 <- pow(sd.co.betalpsi12,-2)
 sd.co.betalpsi12 ~ dunif(0,8)
 tau.co.betalpsi13 <- pow(sd.co.betalpsi13,-2)
 sd.co.betalpsi13 ~ dunif(0,8)
 tau.co.lp <- pow(sd.co.lp,-2)   
 sd.co.lp ~ dunif(0,8)

 tau.co.betalp1 <- pow(sd.co.betalp1,-2)
 sd.co.betalp1 ~ dunif(0,8)
 tau.co.betalp2 <- pow(sd.co.betalp2,-2)
 sd.co.betalp2 ~ dunif(0,8)

# priors de los efectos especie-específicos en la ocupacion y detecc
# sp representa el parámetro por especie
# co representa el parámetro de la comunidad 
  for (k in 1:nspec){
    mu.sp.lpsi[k] ~ dnorm(mu.co.lpsi, tau.co.lpsi)  #queda
    mu.sp.betalpsi1[k] ~ dnorm(mu.co.betalpsi1, tau.co.betalpsi1)
    mu.sp.betalpsi2[k] ~ dnorm(mu.co.betalpsi2, tau.co.betalpsi2)
    mu.sp.betalpsi3[k] ~ dnorm(mu.co.betalpsi3, tau.co.betalpsi3)
    mu.sp.betalpsi4[k] ~ dnorm(mu.co.betalpsi4, tau.co.betalpsi4)
    mu.sp.betalpsi5[k] ~ dnorm(mu.co.betalpsi5, tau.co.betalpsi5)
    mu.sp.betalpsi6[k] ~ dnorm(mu.co.betalpsi6, tau.co.betalpsi6)
    mu.sp.betalpsi7[k] ~ dnorm(mu.co.betalpsi7, tau.co.betalpsi7)
    mu.sp.betalpsi8[k] ~ dnorm(mu.co.betalpsi8, tau.co.betalpsi8)
    mu.sp.betalpsi9[k] ~ dnorm(mu.co.betalpsi9, tau.co.betalpsi9)
    mu.sp.betalpsi10[k] ~ dnorm(mu.co.betalpsi10, tau.co.betalpsi10)
    mu.sp.betalpsi11[k] ~ dnorm(mu.co.betalpsi11, tau.co.betalpsi11)
    mu.sp.betalpsi12[k] ~ dnorm(mu.co.betalpsi12, tau.co.betalpsi12)
    mu.sp.betalpsi13[k] ~ dnorm(mu.co.betalpsi13, tau.co.betalpsi13)

    mu.sp.lp[k] ~ dnorm(mu.co.lp, tau.co.lp)
    mu.sp.betalp1[k]~ dnorm(mu.co.betalp1, tau.co.betalp1)
    mu.sp.betalp2[k]~ dnorm(mu.co.betalp2, tau.co.betalp2)
    
    #Precision de lpsi y lp  
    tau.sp.lpsi[k] <- pow(sd.sp.lpsi[k],-2)   
    sd.sp.lpsi[k] ~ dunif(0,8)
    tau.sp.lp[k] <- pow(sd.sp.lp[k],-2)
    sd.sp.lp[k] ~ dunif(0,8)

    #Precision de los betas especie-especificos
    tau.sp.betalpsi1[k] <- pow(sd.sp.betalpsi1[k],-2)
    sd.sp.betalpsi1[k] ~ dunif(0,8)
    tau.sp.betalpsi2[k] <- pow(sd.sp.betalpsi2[k],-2)
    sd.sp.betalpsi2[k] ~ dunif(0,8)
    tau.sp.betalpsi3[k] <- pow(sd.sp.betalpsi3[k],-2)
    sd.sp.betalpsi3[k] ~ dunif(0,8)
    tau.sp.betalpsi4[k] <- pow(sd.sp.betalpsi4[k],-2)
    sd.sp.betalpsi4[k] ~ dunif(0,8)
    tau.sp.betalpsi5[k] <- pow(sd.sp.betalpsi5[k],-2)
    sd.sp.betalpsi5[k] ~ dunif(0,8)
    tau.sp.betalpsi6[k] <- pow(sd.sp.betalpsi6[k],-2)
    sd.sp.betalpsi6[k] ~ dunif(0,8)
    tau.sp.betalpsi7[k] <- pow(sd.sp.betalpsi7[k],-2)
    sd.sp.betalpsi7[k] ~ dunif(0,8)
    tau.sp.betalpsi8[k] <- pow(sd.sp.betalpsi8[k],-2)
    sd.sp.betalpsi8[k] ~ dunif(0,8)
    tau.sp.betalpsi9[k] <- pow(sd.sp.betalpsi9[k],-2)
    sd.sp.betalpsi9[k] ~ dunif(0,8)
    tau.sp.betalpsi10[k] <- pow(sd.sp.betalpsi10[k],-2)
    sd.sp.betalpsi10[k] ~ dunif(0,8)
    tau.sp.betalpsi11[k] <- pow(sd.sp.betalpsi11[k],-2)
    sd.sp.betalpsi11[k] ~ dunif(0,8)
    tau.sp.betalpsi12[k] <- pow(sd.sp.betalpsi12[k],-2)
    sd.sp.betalpsi12[k] ~ dunif(0,8)
    tau.sp.betalpsi13[k] <- pow(sd.sp.betalpsi13[k],-2)
    sd.sp.betalpsi13[k] ~ dunif(0,8)

    tau.sp.betalp1[k] <- pow(sd.sp.betalp1[k],-2)
    sd.sp.betalp1[k] ~ dunif(0,8)
    tau.sp.betalp2[k] <- pow(sd.sp.betalp2[k],-2)
    sd.sp.betalp2[k] ~ dunif(0,8)
    
    #Previas de las ordenadas
    lpsi[k] ~ dnorm(mu.sp.lpsi[k], tau.sp.lpsi[k]) #queda
    lp[k] ~ dnorm(mu.sp.lp[k], tau.sp.lp[k])
    
    #Previas de los betalpsi y betalp
    betalpsi1[k] ~ dnorm (mu.sp.betalpsi1[k], tau.sp.betalpsi1[k])
    betalpsi2[k] ~ dnorm (mu.sp.betalpsi2[k], tau.sp.betalpsi2[k])
    betalpsi3[k] ~ dnorm (mu.sp.betalpsi3[k], tau.sp.betalpsi3[k])
    betalpsi4[k] ~ dnorm (mu.sp.betalpsi4[k], tau.sp.betalpsi4[k])
    betalpsi5[k] ~ dnorm (mu.sp.betalpsi5[k], tau.sp.betalpsi5[k])
    betalpsi6[k] ~ dnorm (mu.sp.betalpsi6[k], tau.sp.betalpsi6[k])
    betalpsi7[k] ~ dnorm (mu.sp.betalpsi7[k], tau.sp.betalpsi7[k])
    betalpsi8[k] ~ dnorm (mu.sp.betalpsi8[k], tau.sp.betalpsi8[k])
    betalpsi9[k] ~ dnorm (mu.sp.betalpsi9[k], tau.sp.betalpsi9[k])
    betalpsi10[k] ~ dnorm (mu.sp.betalpsi10[k], tau.sp.betalpsi10[k])
    betalpsi11[k] ~ dnorm (mu.sp.betalpsi11[k], tau.sp.betalpsi11[k]) 
    betalpsi12[k] ~ dnorm (mu.sp.betalpsi12[k], tau.sp.betalpsi12[k])
    betalpsi13[k] ~ dnorm (mu.sp.betalpsi13[k], tau.sp.betalpsi13[k])
    betalp1[k] ~ dnorm(mu.sp.betalp1[k],tau.sp.betalp1[k])
    betalp2[k] ~ dnorm(mu.sp.betalp2[k],tau.sp.betalp2[k])
  } #nspec
  
# Modelo ecologico, ocurrencia en el sitio i 
#Likelihood

for (k in 1:nspec){
#loop sobre especies
    
    for (m in 1:nyears){ 
#loop sobre anios

        for (i in 1:nsite) {                                           
#loop sobre sitios

          for (l in 1:nseasons){
#loop sobre estacion

              logit(psi[i,l,k,m]) <- lpsi[k] + betalpsi1[k]*summer[i,m,l] + betalpsi2[k]*fall[i,m,l] + betalpsi3[k]*winter[i,m,l] + betalpsi4[k]*spring[i,m,l] + betalpsi5[k]*tipo[i,m] + betalpsi6[k]*zanja[i,m] + betalpsi7[k]*galpon[i,m] + betalpsi8[k]*tambo[i,m] + betalpsi9[k]*vegetacion[i,m] + betalpsi10[k]*comedero[i,m] + betalpsi11[k]*silo[i,m] + betalpsi12[k]*foza[i,m] + betalpsi13[k]*tipo[i,m]
       
              z[i,k,l,m] ~ dbern(psi[i,l,k,m])
  
# Modelos de observacion para sitio nsite = i, repeticion nrep = j, nspec=k, estacion nseasons = l, año nyears = m

            for (j in 1:nrep) {  
#loop sobre repeticiones nrep

               logit(p[i,j,k,l,m]) <-  lp[k] + betalp1[k]*temp[i,j,l,m] + betalp2[k]*noches[i,j]

               mu.p[i,j,k,l,m] <- p[i,j,k,l,m]*z[i,k,l,m]

               y[i,j,k,l,m] ~ dbern(mu.p[i,j,k,l,m])

               } #nrep
           } #nseasons
         } #nsite
      } #nyears
    } #nspec

} #model
",fill=TRUE)
sink()

# valores iniciales
zst <- array(1,c(nsite, nspec, nseasons, nyears))
inits <- function() list(z = zst)


# parametros a monitorear
# parámetros especie-especificos 1 
params1 <- c("lp", "betalp1", "betalp2",  "lpsi", "betalpsi1","betalpsi2", "betalpsi3","betalpsi4","betalpsi5","betalpsi6","betalpsi7")

params2 <- c("betalpsi8","betalpsi9","betalpsi10","betalpsi11","betalpsi12","betalpsi13")

params3 <- c("mu.sp.lpsi", "mu.sp.betalpsi1","mu.sp.betalpsi2", "mu.sp.betalpsi3", "mu.sp.betalpsi4","mu.sp.betalpsi5","mu.sp.betalpsi6","mu.sp.betalpsi7")

params4 <- c("mu.sp.betalpsi8","mu.sp.betalpsi9", "mu.sp.betalpsi10","mu.sp.betalpsi11","mu.sp.betalpsi12","mu.sp.lp","mu.sp.betalp1","mu.sp.betalp2")

#Parametros de la comunidad
params5 <- c( "mu.co.lpsi", "mu.co.betalpsi1", "mu.co.betalpsi2","mu.co.betalpsi3","mu.co.betalpsi4","mu.co.betalpsi5","mu.co.betalpsi6","mu.co.betalpsi7") 

params6 <- c("mu.co.betalpsi8", "mu.co.betalpsi9", "mu.co.betalpsi10", "mu.co.betalpsi11", "mu.co.betalpsi12", "mu.co.lp", "mu.co.betalp1", "mu.co.betalp2")


# ajustes de MCMC
ni <- 500000 #60000
nt <- 10 
nb <- 50000 #100
nc <- 3

beep(3)

out1 = jags(my.data, inits, params1, "mimodelo.jags", n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, parallel=TRUE);beep(3)

out2 = jags(my.data, inits, params2, "mimodelo.jags", n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, parallel=TRUE);beep(3)

out3 = jags(my.data, inits, params3, "mimodelo.jags", n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, parallel=TRUE);beep(3)

out4 = jags(my.data, inits, params4, "mimodelo.jags", n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, parallel=TRUE);beep(3)

out5 = jags(my.data, inits, params5, "mimodelo.jags", n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, parallel=TRUE); beep(3)

out6 = jags(my.data, inits, params6, "mimodelo.jags", n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, parallel=TRUE); beep(2)

#guardo resultados
save(out1, file='mi_modelo_1.rda') 
save(out2, file='mi_modelo_2.rda')
save(out3, file='mi_modelo_3.rda')
save(out4, file='mi_modelo_4.rda')
save(out5, file='mi_modelo_5.rda')
save(out6, file='mi_modelo_6.rda')


load('mi_modelo_1.rda')
load('mi_modelo_2.rda') 
load('mi_modelo_3.rda') 
load('mi_modelo_4.rda') 
load('mi_modelo_5.rda') 
load('mi_modelo_6.rda') 

#los outs que realmente me interesan
out1
out2
out3

#Valores de rhat
rhat1 <- t(abind(out1$Rhat$lp, out1$Rhat$lpsi,out1$Rhat$betalp1, out1$Rhat$betalp2, out1$Rhat$betalpsi1,out1$Rhat$betalpsi2,out1$Rhat$betalpsi3,out1$Rhat$betalpsi4,out1$Rhat$betalpsi5, out1$Rhat$betalpsi6,out1$Rhat$betalpsi7, along = 2)) 

dimnames(rhat1) <- list(parámetro=c("lp","lpsi","betalp1","betalp2","betalpsi1","betalpsi2","betalpsi3","betalpsi4","betalpsi5","betalpsi6","betalpsi7"), Especie= species.list)


#grafico la convegencia de las cadenas (hacer "enter" en el gráfico para que pasen)
x11(par(mfrow = c(3, 3)))
traceplot(out1)

x11(par(mfrow = c(3, 3)))
traceplot(out2)

x11(par(mfrow = c(3, 3)))
traceplot(out3)

x11(par(mfrow = c(3, 3)))
traceplot(out4)

x11(par(mfrow = c(3, 3)))
traceplot(out5)

x11(par(mfrow = c(3, 3)))
traceplot(out6)

#plot de las posteriores
x11()
densityplot(out1, parameters=c(out1$parameters), layout=NULL, ask=NULL)

x11()
densityplot(out2, parameters=c(out1$parameters), layout=NULL, ask=NULL)

x11()
densityplot(out3, parameters=c(out1$parameters), layout=NULL, ask=NULL)

x11()
densityplot(out4, parameters=c(out1$parameters), layout=NULL, ask=NULL)

x11()
densityplot(out5, parameters=c(out1$parameters), layout=NULL, ask=NULL)

x11()
densityplot(out6, parameters=c(out1$parameters), layout=NULL, ask=NULL)

#Whiskerplots
whiskerplot(out1, parameters=c(out1$parameters), quantiles=c(0.025,0.975), zeroline=TRUE)

#ploteos de posteriores que no solapan el cero
x11()
plotdens(out1, p="betalp2", legend = T, legendpos = "topleft", legendnames = c("sp_1","sp_2","sp_3","sp_4","sp_5","sp_6","sp_7"), minCI = 0.95, main = "Efecto de la noche de muestreo sobre la detectabilidad (betalp2)")
abline(v=0, col="red", lty=2, lwd=3)

#detección por noche
x11()
caterpillar(out1, p="betalp2", col = c(1:7), xlab = "especie", ylab="logit(psi)")
abline(h=0, col="red", lty=2, lwd=3)

#Ocupación por estación
x11()
par(mfrow=c(2,2))
caterpillar(out1, p="betalpsi1", col = c("1e8fffff","ffa400ff","006400ff","a020f0ff","ff0000ff","d2681eff","bebebeff"), xlab = "especie", ylab="logit(psi)", main = "Uso en verano")
abline(h=0, col="red", lty=2, lwd=3)

caterpillar(out1, p="betalpsi2", col = c(1:7), xlab = "especie", ylab="logit(psi)", main = "Uso en otoño")
abline(h=0, col="red", lty=2, lwd=3)

caterpillar(out1, p="betalpsi3", col = c(1:7), xlab = "especie", ylab="logit(psi)", main = "Uso en invierno")
abline(h=0, col="red", lty=2, lwd=3)

caterpillar(out1, p="betalpsi4", col = c(1:7), xlab = "especie", ylab="logit(psi)", main = "Uso en primavera")
abline(h=0, col="red", lty=2, lwd=3)


#Ocupación por ambiente
#
dev.off()
x11()
par(mfrow=c(2,3))
caterpillar(out1, p="betalpsi6", col = c("dodgerblue","orange","darkgreen","purple","red","chocolate","grey"), xlab = "", ylab="logit(simbol(psi))", main = "(a) Zanjas", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, xax=c(expression(italic("Rn"),italic("Aa"), italic("Mm"), italic("Rr"), italic("Cl"), italic("Da"), italic("Of")))) 
abline(h=0, col="red", lty=2, lwd=3)


caterpillar(out1, p="betalpsi7", col = c(1:7), xlab = "especie", ylab="logit(psi)", main = "Uso en ambiente galpón",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5,xax=c(expression(italic("R. norvegicus"),italic("A. azarae"), italic("M. musculus"),italic("R. rattus"), italic("C. laucha"), italic("D. albiventris"), italic("O. flavescens"))), xaxp=c(1,7,0.5))
abline(h=0, col="red", lty=2, lwd=3)

caterpillar(out2, p="betalpsi8", col = c(1:7), xlab = "especie", ylab="logit(psi)", main = "Uso en ambiente tambo", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5,xax=c(expression(italic("R. norvegicus"),italic("A. azarae"), italic("M. musculus"),italic("R. rattus"), italic("C. laucha"), italic("D. albiventris"), italic("O. flavescens"))))
abline(h=0, col="red", lty=2, lwd=3)

caterpillar(out2, p="betalpsi11", col = c(1:7), xlab = "especie", ylab="logit(psi)", main = "Uso en ambiente silo", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5,xax=c(expression(italic("R. norvegicus"),italic("A. azarae"), italic("M. musculus"),italic("R. rattus"), italic("C. laucha"), italic("D. albiventris"), italic("O. flavescens"))))
abline(h=0, col="red", lty=2, lwd=3)

caterpillar(out2, p="betalpsi12", col = c(1:7), xlab = "especie", ylab="logit(psi)", main = "Uso en ambiente fosa", cex.lab=1.5, cex.axis=1.5, cex.sub=1.5,xax=c(expression(italic("R. norvegicus"),italic("A. azarae"), italic("M. musculus"),italic("R. rattus"), italic("C. laucha"), italic("D. albiventris"), italic("O. flavescens"))))
abline(h=0, col="red", lty=2, lwd=3)

####Efectos "marginales" por especie en la detección

#Funcion para los CRI 95%

cri<-function(x) quantile(x,prob=c(0.025,0.975))

#Extraigo el número de muestras de las MCMC
nsamp<-out1$mcmc.info$n.samples

#Hago una secuencia entre los valores máximos y mínimos de la variable "noche"
o.noches <- seq(min(1),max(3), ,400)  

#Extraigo los parámetros y sus valores
tmp<-out1$sims.list

#veo la estructura (cada columna es una especie)
str(tmp)
str(tmp$lp)
str(tmp$betalp2)

#Abro ventana gráfica para poner los 7 gráficos de detección
dev.off()
x11()
par(mfrow=c(3, 3))
### Detección especie 1
#Extraigo la ordenada de la sp_1
lp_1<- tmp$lp[,1]

#Extraigo el beta de la sp_1
betalp_1<- tmp$betalp2[,1]

#Hago un reaglego para llenar con valores. Tiene tres dimensiones: 400 valores de filas, nsamp de muestras como columnas y 3 dimensiones mas (como hojas)
pred<- array(NA, dim=c(400,nsamp,3))

#lleno la dimension de nsamp con los valores
for(i in 1:nsamp){
  pred[,i,3] <- plogis(lp_1[i]+betalp_1[i]*o.noches)}

#Aplico la función del CRI 95% a los valores obtenidos
cri3<-apply(pred[,,3],1,cri)

#Ploteo
plot(o.noches, apply(pred[,,3],1,mean), ylab="Probabilidad de detección", xlab="Noche de muestreo",
     type="l",lty=1,lwd=3,col="dodgerblue",ylim=c(0,1),frame=F,xaxt="n",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5 )
axis(1, at=o.noches,labels=round(o.noches,1), cex=1.5, tck=0)
polygon(x= c(o.noches, rev(o.noches)), y= c(cri3[1,], rev(cri3[2,])), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.10), border = NA)


### Detección especie 2
#Extraigo la ordenada de la sp_2
lp_2<- tmp$lp[,2]

#Extraigo el beta de la sp_2
betalp_2<- tmp$betalp2[,2]

#Hago un reaglego para llenar con valores. Tiene tres dimensiones: 400 valores de filas, nsamp de muestras como columnas y 3 dimensiones mas (como hojas)
pred2<- array(NA, dim=c(400,nsamp,3))

#lleno la dimension de nsamp con los valores
for(i in 1:nsamp){
  pred2[,i,3] <- plogis(lp_2[i]+betalp_2[i]*o.noches)}

#Aplico la función del CRI 95% a los valores obtenidos
cri32<-apply(pred2[,,3],1,cri)

plot(o.noches, apply(pred2[,,3],1,mean), ylab="Probabilidad de detección", xlab="Noche de muestreo",
     type="l",lty=1,lwd=3,col="orange",ylim=c(0,1),frame=F,xaxt="n",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5 )
axis(1, at=o.noches,labels=round(o.noches,1), cex=1, tck=0)
polygon(x= c(o.noches, rev(o.noches)), y= c(cri32[1,], rev(cri32[2,])), 
        col =  adjustcolor("orange", alpha.f = 0.10), border = NA)

### Detección especie 3
#Extraigo la ordenada de la sp_3
lp_3<- tmp$lp[,3]

#Extraigo el beta de la sp_3
betalp_3<- tmp$betalp2[,3]

#Hago un reaglego para llenar con valores. Tiene tres dimensiones: 400 valores de filas, nsamp de muestras como columnas y 3 dimensiones mas (como hojas)
pred3<- array(NA, dim=c(400,nsamp,3))

#lleno la dimension de nsamp con los valores
for(i in 1:nsamp){
  pred3[,i,3] <- plogis(lp_3[i]+betalp_3[i]*o.noches)}

#Aplico la función del CRI 95% a los valores obtenidos
cri33<-apply(pred3[,,3],1,cri)

#ploteo
plot(o.noches, apply(pred3[,,3],1,mean), ylab="Probabilidad de detección", xlab="Noche de muestreo",
     type="l",lty=1,lwd=3,col="darkgreen",ylim=c(0,1),frame=F,xaxt="n",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5 )
axis(1, at=o.noches,labels=round(o.noches,1), cex=1, tck=0)
polygon(x= c(o.noches, rev(o.noches)), y= c(cri33[1,], rev(cri33[2,])), 
        col =  adjustcolor("darkgreen", alpha.f = 0.10), border = NA)

### Detección especie 4
#Extraigo la ordenada de la sp_4
lp_4<- tmp$lp[,4]

#Extraigo el beta de la sp_4
betalp_4<- tmp$betalp2[,4]

#Hago un reaglego para llenar con valores. Tiene tres dimensiones: 400 valores de filas, nsamp de muestras como columnas y 3 dimensiones mas (como hojas)
pred4<- array(NA, dim=c(400,nsamp,3))

#lleno la dimension de nsamp con los valores
for(i in 1:nsamp){
  pred4[,i,3] <- plogis(lp_4[i]+betalp_4[i]*o.noches)}

#Aplico la función del CRI 95% a los valores obtenidos
cri34<-apply(pred4[,,3],1,cri)

#Ploteo
plot(o.noches, apply(pred4[,,3],1,mean), ylab="Probabilidad de detección", xlab="Noche de muestreo",
     type="l",lty=1,lwd=3,col="purple",ylim=c(0,1),frame=F,xaxt="n",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5 )
axis(1, at=o.noches,labels=round(o.noches,1), cex=1, tck=0)
polygon(x= c(o.noches, rev(o.noches)), y= c(cri34[1,], rev(cri34[2,])), 
        col =  adjustcolor("purple", alpha.f = 0.10), border = NA)


### Detección especie 5
#Extraigo la ordenada de la sp_5
lp_5<- tmp$lp[,5]

#Extraigo el beta de la sp_4
betalp_5<- tmp$betalp2[,5]

#Hago un reaglego para llenar con valores. Tiene tres dimensiones: 400 valores de filas, nsamp de muestras como columnas y 3 dimensiones mas (como hojas)
pred5<- array(NA, dim=c(400,nsamp,3))

#lleno la dimension de nsamp con los valores
for(i in 1:nsamp){
  pred5[,i,3] <- plogis(lp_5[i]+betalp_5[i]*o.noches)}

#Aplico la función del CRI 95% a los valores obtenidos
cri35<-apply(pred5[,,3],1,cri)

#Ploteo
plot(o.noches, apply(pred5[,,3],1,mean), ylab="Probabilidad de detección", xlab="Noche de muestreo",
     type="l",lty=1,lwd=3,col="red",ylim=c(0,1),frame=F,xaxt="n",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5 )
axis(1, at=o.noches,labels=round(o.noches,1), cex=1, tck=0)
polygon(x= c(o.noches, rev(o.noches)), y= c(cri35[1,], rev(cri35[2,])), 
        col =  adjustcolor("red", alpha.f = 0.10), border = NA)

### Detección especie 6
#Extraigo la ordenada de la sp_6
lp_6<- tmp$lp[,6]

#Extraigo el beta de la sp_4
betalp_6<- tmp$betalp2[,6]

#Hago un reaglego para llenar con valores. Tiene tres dimensiones: 400 valores de filas, nsamp de muestras como columnas y 3 dimensiones mas (como hojas)
pred6<- array(NA, dim=c(400,nsamp,3))

#lleno la dimension de nsamp con los valores
for(i in 1:nsamp){
  pred6[,i,3] <- plogis(lp_6[i]+betalp_6[i]*o.noches)}

#Aplico la función del CRI 95% a los valores obtenidos
cri36<-apply(pred6[,,3],1,cri)

#Ploteo
plot(o.noches, apply(pred6[,,3],1,mean), ylab="Probabilidad de detección", xlab="Noche de muestreo",
     type="l",lty=1,lwd=3,col="chocolate",ylim=c(0,1),frame=F,xaxt="n",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5 )
axis(1, at=o.noches,labels=round(o.noches,1), cex=1, tck=0)
polygon(x= c(o.noches, rev(o.noches)), y= c(cri36[1,], rev(cri36[2,])), 
        col =  adjustcolor("chocolate", alpha.f = 0.10), border = NA)

### Detección especie 7
#Extraigo la ordenada de la sp_7
lp_7<- tmp$lp[,7]

#Extraigo el beta de la sp_4
betalp_7<- tmp$betalp2[,7]

#Hago un reaglego para llenar con valores. Tiene tres dimensiones: 400 valores de filas, nsamp de muestras como columnas y 3 dimensiones mas (como hojas)
pred7<- array(NA, dim=c(400,nsamp,3))

#lleno la dimension de nsamp con los valores
for(i in 1:nsamp){
  pred7[,i,3] <- plogis(lp_7[i]+betalp_7[i]*o.noches)}

#Aplico la función del CRI 95% a los valores obtenidos
cri37<-apply(pred7[,,3],1,cri)

#Ploteo
plot(o.noches, apply(pred7[,,3],1,mean), ylab="Probabilidad de detección", xlab="Noche de muestreo",
     type="l",lty=1,lwd=3,col="grey",ylim=c(0,1),frame=F,xaxt="n",cex.lab=1.5, cex.axis=1.5, cex.sub=1.5 )
axis(1, at=o.noches,labels=round(o.noches,1), cex=1, tck=0)
polygon(x= c(o.noches, rev(o.noches)), y= c(cri37[1,], rev(cri37[2,])), 
        col =  adjustcolor("grey", alpha.f = 0.10), border = NA)

# leyenda
legend(x=4, y=1.5, pch=15,cex=0.6,border=NA,bty="n",
       legend= substitute(paste(italic('Sp_1')))
       ,col=c("dodgerblue"), xpd=NA)
legend(x=4, y=1.30, pch=15,cex=0.6,border=NA,bty="n",
       legend=substitute(paste(italic('Sp_2')))
       ,col=c("orange"), xpd=NA)
legend(x=4, y=1.10, pch=15,cex=0.6,border=NA,bty="n",
       legend= substitute(paste(italic('Sp_3')))
       ,col=c("darkgreen"), xpd=NA)
legend(x=4, y=0.90, pch=15,cex=0.6,border=NA,bty="n",
       legend= substitute(paste(italic('Sp_4')))
       ,col=c("purple"), xpd=NA)
legend(x=4, y=0.70, pch=15,cex=0.6,border=NA,bty="n",
       legend= substitute(paste(italic('Sp_5')))
       ,col=c("red"), xpd=NA)
legend(x=4, y=0.50, pch=15,cex=0.6,border=NA,bty="n",
       legend= substitute(paste(italic('Sp_6')))
       ,col=c("chocolate"), xpd=NA)
legend(x=4, y=0.30, pch=15,cex=0.6,border=NA,bty="n",
       legend= substitute(paste(italic('Sp_7')))
       ,col=c("grey"), xpd=NA)

R.Version()
RStudio.Version()
