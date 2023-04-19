##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########         Ejemplo 1 - Módulo 2         #########
###########             Probabilidades           #########
##########################################################

# Para las funciones de probabilidad estándares, en r hay funciones que permiten
# generar valores dentro de esa probabilidad (r), calcular probablidades (d),
# calcular probabilidades acumuladas (p) y determinar cuantiles (q).

# Cada función tiene nombres y argumentos específicos:
# Normal: rnorm(), dnorm(), pnorm(), qnorm() 
# Poisson: rpois(), dpois(), ppois(), qpois()
# Binomial: rbin(), dbin(), pbin(), qbin()

# DISTRIBUCIONES DISCRETAS----

# ---------#
#  Poisson #
# ---------#
# Simulemos conteos de individuos en 20 grillas, que poseen una media de
# 3 individuos por grilla 

y <- rpois(n = 20, lambda = 3)
hist(y, xlab = "Valores", ylab = "Frecuencia", main = "",
     bty =  "n",  las = 1)

# Para ver en detalle argumentos para  modificar los gráficos los puede
# explorar en ?par

# Para calcular la probabilidad de obtener 5 conteos de individuos por grilla utilizamos
dpois(5, lambda = 3) 

# Para calcular la probabilidad de obtener un valor menor o igual a 5,
ppois(5, lambda = 3)

# Actividad ----
# 1) ¿Cómo cambia la forma de la distribución de Poisson a medida que cambia λ?
# 2) Genere 10 valores con un valor de lambda de 1.4 
# 3) ¿Cuál es la probabilidad de obtener y = 0 para para λ = 1.3?
# 4) ¿Cuál es la probabilidad de y > 5 para λ=1.3?

# DISTRIBUCIONES CONTINUAS----
# ---------#
#  Normal  #
# ---------#

# 1) Genere 10 valores aleatorios con una distribución normal con media 5 y 
# desvío 0.3 
# 2) Calcule el intervalo del 80 y 95 % 
