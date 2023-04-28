# remueve datos viejos
rm(list=ls(all=TRUE))

# indicar directorio de trabajo
setwd('C:\\Users\\andrea\\Documents\\GitHub\\Curso-Ocupacion23\\EjerciciosEjemplos\\Modulo5-OcupacionSimple')

#leer el csv
data<-read.csv('datos_MATRIZ.csv', sep = ",", header = T)
#explorar datos
head(data, 10)

str(data)

# remplazar NA por cero
data[is.na(data)] <- 0
head(data)

datasample<-read.csv('datOS_LONG_COMPLETOS.csv', sep = ";", header = T)
head(datasample, 10)

#cargo los paquetes que voy a usar
library(reshape)
library(plyr)

# extraigo el nombre de las especies
especies.names<-as.vector(names(data[6:18]))
especies.names

#armar el vector del nombre de las especies para cada uno de los 72 
#registros (36 sitios x 2 repeticiones) y lo ordeno alfabeticamente
sp.names<-sort(factor(rep(especies.names,72)))
class(sp.names)


# extraigo los datos de registros y los pongo en una sola fila de registros
sp.obs<-data[6:18]

# los ordeno alfabeticamente de modo q cada especie aparezca junta
sp.obs<-sp.obs[ , order(names(sp.obs))]
head(sp.obs)


# asi paso los datos de columna a un solo vector con "unlist"
especies.obs<-unlist(sp.obs,use.names = FALSE)
class(especies.obs)


#72 registros * 13 especies = 936 registros para cehquear q esta bien
nreg<-length(especies.obs)
nreg


# unir las especies con sus registros
sp<-data.frame(sp.names, especies.obs)


# colectar informacion de cada punto
data.info<-data[1:72,1:5]
head(data.info)


# repetirla para cada una de las 13 especcies
rep.data.info<-do.call("rbind", replicate(13, data.info, simplify = FALSE))


#armar la planilla longitudinal con los registros para cada especie x sitio
data.aves.long<-data.frame(rep.data.info, sp)
str(data.aves.long)



#exportar tabla
write.csv(data.aves.long, file="data.aves.long.csv")


## Formato longitudinal incompleto

# Otra opcion cuando se ingresan los datos, es ingresar por cada punto y repeticion, 
# las especies observadas.
data2<-read.csv("datos_LONG_INCOMPLETOS.csv") #,header = T, sep = ";")
head(data2, 10)


# Para analisis de ocupacion de una sola especie, podemos seleccionar un subset de especies
data.LEILO <- subset(data.aves.long,data.aves.long$sp.names == "LEILO")
head(data.LEILO)

## CONSTRUCCION DE LA MATRIZ Y
#  Como construir la matriz "y" para analizar con Occupancy

# Reordeno los datos con el paquete "reshape". Los datos de deteccion/no deteccion son reordenados en una
# matriz de N dimensiones,
# En el ejemplo que estabmos trabajando tenemos 4 dimensiones: El predio l es la 1er dimension, luego el
# punto i, la repeticion j y la especie k

data.melt=melt(data.aves.long,id.var=c("PREDIO", "PUNTO", "REP", "sp.names"), measure.var="especies.obs")
y1=cast(data.melt, PREDIO ~ PUNTO ~ REP ~ sp.names)

# Verifico las dimensiones
dim(y1)


#chequeo con datos originales
y1[,,,1]


# si hay algun NA le pone cero (esto no hay que hacerlo siempre)
y1[is.na(y1)] <- 0

# Si trabajo con occupancy paso todos los numeros a 0 o 1
y<-aaply(y1,1,function(x) {x[x>1]<-1;x})
dim(y)


# Nombres de cada dimension
lista.predios<-as.vector(names(y[,1,1,1]))
lista.especies<-as.vector(names(y[1,1,1,]))
lista.sitios<-as.vector(names(y[1,,1,1]))
lista.reps<-as.vector(names(y[1,1,,1]))

# Dimensiones de cada elemento
npredios<-length(lista.predios)
nspec<-length(lista.especies)
nsites<-length(lista.sitios)
nreps<-length(lista.reps)

#si quiero reordenar las dimensiones
#predio, punto, rep, sp
#sp, predio, punto, rep
y2<-aperm(y,c(4,1,2,3))
y<-y2




