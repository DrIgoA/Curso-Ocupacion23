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

# Funcion de verosimilitud (likelihood)
L = function(p,k,n) p^k*(1-p)^(n-k)

# Funcion logaritmo de la verosimilitud (log-likelihood)
l = function(p,k,n) k*log(p) + (n-k)*log(1-p)

# Probabilidad desde 0 hasta 1
p = seq(0,1,0.001)

# Podemos ver los datos generados
p
dbinom(k,n,p)
L(p,K,N)

# Graficando Probabilidad binomial de k dado p
plot(p, dbinom(k,n,p), type="l", main="Probabilidad binomial de k dado p")

# Graficando la funcion de verosimilitud
windows()
plot(p, L(p,K,N), type="l")

# Graficando la funci?n del logaritmo de la verosimilitud (log lokelihood)
windows()
plot(p, l(p,K,N), type="l", main="Log likelihood")

# Las funciones de optimizaci?n en R encuentran el m?nimo, no el m?ximo. 
# asique debemos crear nuevas funciones que nos den la likelihood y log-likelihoods negativas
# y luego minimizarlas

# Likelihood negativa:
mL = function(p,k,n) -p^k*(1-p)^(n-k)

# log-likelihood negativa:
ml = function(p,k,n) -(k*log(p) + (n-k)*log(1-p))

# Utilizar 'optimize' para encontar la soluci?n
optimize(mL, interval = c(0,1), k=K, n=N)
optimize(ml, interval = c(0,1), k=K, n=N)

k/n
