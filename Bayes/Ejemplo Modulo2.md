---
title: "Ejemplos Modulo 2"
author: "Goijman, Serafini, Contreras"
date: '2023-03-14'
output: html_document
---

Cargamos la librería `jagsUI` para correr JAGS 
```{r setup, error=TRUE}
library(jagsUI)    #paquete JAGS
```
Ejemplo tamaño medio de zorzales (n=10)

Agrupar los datos en una lista para que pueda ser leida por JAGS
```{r, error=TRUE}
data <- list(size = c(7.9,8.1,11,10.6,9.2,8,9.8,10.1,10.9,9))
```

Modelo para pasarle a `jags`
es muy importante el nombre del archivo (también podríamos escribirlo en un archivo de texto)
```{r, error=TRUE, include=FALSE}
sink("zorzales.jags")
cat("

model
{
# previas no informativas
mean ~ dnorm (0, 1.0E-6)  # tamaño medio de los zorzales
varianza ~ dlnorm (0 ,1.0E-6)  # varianza tamaño zorzales

prec <- 1/varianza             # pasar de varianza a precision 

for (i in 1:10)            # para cada uno de los Zorzales
{
size[i] ~ dnorm (mean, prec)   # tamaño zorzal trazado de una distribución normal
}
}

",fill = TRUE)
sink()
```

Valores Iniciales para correr las cadenas de Markov (MCMC)
```{r, error=TRUE}
inits <- function() list(varianza=100, mean = 10)

```

Seteo para correr MCMC donde incicamos cuántos valores vamos a ir descartado en una secuencia (“thin”), qué largo tiene el “burn in” y cuántas cadenas queremos correr.
```{r, error=TRUE}
ni <- 10000  # número de iteraciones
nt <- 2      # tasa de thining
nb <- 1000   # iteraciones para el burn in
nc <- 3      # número de cadenas que corremos
```

Le decimos a `jags`que parámetros nos interesa monitorear
```{r, error=TRUE}
parameters <- c("mean", "varianza")
```
Para llamar a `JAGS` desde `R` usamos el paquete `jagsUI` que cargamos anteriormente y la función `jags` de ese paquete
```{r, error=TRUE, include=FALSE}
out <- jags(data, inits, parameters, "zorzales.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
```

Para ver los resultados le pedimos a `R` que muestre la salida de la función jags llamando a "out", que es el nombre que le dimos. 

```{r, error=TRUE}
out
```
Se indica el nombre del modelo que se usó, cuántas cadenas se simularon y cuántas iteraciones hubo. También tenemos el de "muestras totales" teniendo en cuenta el burn-in y el thinning. 
Luego podemos ver tabla con la lista de parámetros que le pedimos que registre (en este caso la media y la varianza) y la devianza. En la tabla aparece la media, desvío y cuantiles estimados a partir de las cadenas de Markov. 
También aparece una columna (overlap0) que nos dice si la posterior incluye al cero o no, y otra (f) que nos dice qué fracción de la posterior es del mismo signo que la media. 
Rhat estima si las cadenas convergieron a una distribución estable y n.eff estima el número efectivo de muestras de la posterior que surgen de las cadenas. 
Antes de sacar cualquier conclusión, tenemos que chequear que haya convergencia de las cadenas (Rhat ≤1.1). 

También podemos ver todo lo que tiene guardada la salida con 
```{r, error=TRUE}
str(out)
```

Podemos (y debemos) inspeccionar visualmente la salida del MCMC
```{r, error=TRUE}
plot(out)
```

Cuando tengamos muchos parámetros, podemos "llamar" a aquellos que nos interesen
```{r, error=TRUE}
traceplot(out, parameters=c("mean", "varianza")) 
```
