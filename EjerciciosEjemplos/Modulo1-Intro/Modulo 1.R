##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########           Ejemplo Muestreo           #########
###########                                      #########
##########################################################
########              Ejemplos de:                  ######
########  Applied hierarchical modeling in ecology  ######
########  Modeling distribution, abundance and      ######
########  species richness using R and BUGS         ######
########  Volume 1: Prelude and Static models       ######
########      Marc K?ry & J. Andy Royle             ######
##########################################################

# La siguiente simulación permite experimentar con la relación entre patrones de 
# puntos, la ocurrencia y abundancia de organismos en un ?rea determinada en 
# relación a los sistemas de muestreos seleccionados para poder hacer el 
# revelamiento de la especie particular.

# Esta funci?n se encuentra en la librer?a AHMbook
# Antes de poder utilizarla, es necesario instalarla
 
install.packages("AHMbook")

library(AHMbook)
# La función "sim.fn" simula ubicaciones de plantas o animales en grillas de 
# celdas formando un cuadrante de un tamaño específico (quad.size) de acuerdo
# a una distribución de Poisson, donde los organismos están distribuidos al 
# azar dentro del cuadrante. El proceso se caracteriza por una constante "intensity",
# que es el promedio de animales o plantas (puntos) por unidad de ?rea. 

# Par?metros por default en la funci?n (corresponden a los valores
# que tomar?n los argumentos si no se especifican)
# quad.size = 10, cell.size = 1, intensity = 1

# EJEMPLO ----
set.seed(82)                    # N?mero entero para inicializar la generaci?n de valores aleatorios

tmp <-sim.fn(quad.size=16,       # tama?o total del cuadrante
             cell.size=2,        #tama?o de las celdas de la grilla
             intensity = 0.5)    # constante, valore promedio de animales por celda


# ACTIVIDAD ----
# Modifique la función de simulación para responder las siguientes preguntas
# cambiando el tamaño de la celda de muestreo (cell.size)

# 1 - ?Que observamos en relación a los patrones de abundancia y ocurrencia
# cuando se modifican los tamaños de la celda de muestreo?

# 2 - ?Qué sucede con los patrones de abundancia y ocurrencia cuando las especies
# son raras? ?Se mantiene siempre esta relación?
