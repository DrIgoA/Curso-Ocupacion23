##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########         Ejemplo 2 - Módulo 2         #########
###########               Likelihood             #########
##########################################################
 
# En una muestra de n = 10 individuos se parasitan k = 4 individuos. 

n = 10  #individuos totales
k = 4   #individuos parasitados
 
# Calcula la probabilidad de obtener estos datos dado que p= 0, 0.1, …, 1. 
# Luego haz un gráfico de esta probabilidad como # función de p.

# genera valores de probabilidad desde 0 hasta 1
p <- seq(0,1,0.001)
p       # Podemos ver los datos generados

# Probabilidad Binomial 
dbinom(k,n,p)

# Graficando Probabilidad binomial de k dado p
plot(p, dbinom(k,n,p), type="l", main="Probabilidad binomial de k dado p")

# Actividad ----
# 1) Prueba diferentes valores de tamaño de muestra (n) y número de individuos
# parasitados (k) y observa cómo cambia el gráfico. ¿Cómo puede usar este gráfico
# para encontrar la “mejor suposición” de la prevalencia del parásito en toda la
# población con base en estos datos?

# 2) Intenta aumentar tanto el tamaño de la muestra (n) como el número de individuos
# parasitados (k) manteniendo las mismas proporciones (ej., k/n = 20/8, 30/12,
# 40/16,…). ¿Cómo cambia la forma de la curva cuando aumenta el tamaño de la
# muestra?  ¿Como interpretas esto?


# 3) Calcula también la probabilidad de los datos dados los diferentes valores de p 
# cuando los individuos parasitados aparecen en un orden particular (p. ej., primero
# obtienes 4 parasitados, luego 6 no parasitados) y grafica esto. ¿En qué se diferencian
# las dos curvas? ¿Esto tiene sentido? ¿Obtenemos alguna información sobre el parámetro
# sabiendo el orden en el que tomamos muestras de individuos parasitados y no parasitados? 

# Funcion de verosimilitud (likelihood)
L = function(p,k,n) p^k*(1-p)^(n-k)

# Funcion logaritmo de la verosimilitud (log-likelihood)
l = function(p,k,n) k*log(p) + (n-k)*log(1-p)

L(p,k,n)
# Graficando la funcion de verosimilitud
plot(p, L(p,k,n), type="l")

# Graficando la funcion del logaritmo de la verosimilitud (log lokelihood)
plot(p, l(p,k,n), type="l", main="Log likelihood")

# Las funciones de optimizacion en R encuentran el minimo, no el maximo. 
# asi que debemos crear nuevas funciones que nos den la likelihood y log-likelihoods negativas
# y luego minimizarlas

# Likelihood negativa:
mL = function(p,k,n) -p^k*(1-p)^(n-k)

# log-likelihood negativa:
ml = function(p,k,n) -(k*log(p) + (n-k)*log(1-p))

# Utilizar 'optimize' para encontar la soluci?n
optimize(mL, interval = c(0,1), k=k, n=n)
optimize(ml, interval = c(0,1), k=k, n=n)

k/n
