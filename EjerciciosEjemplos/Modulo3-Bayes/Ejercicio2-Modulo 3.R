##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########      Ejercicio  jerarquico           #########
##########################################################
########              Ejemplo de:                   ######
########          Bayesian Methos for ecology       ######
########            Michael A. MacCarthy            ######
##########################################################

#---------------------------------------------------------------------------
# Estimación de media con distribucion Poisson
# ---------------------------------------------------------------------------

# Investigadores en NY midieron la densidad de arboles en 10 cuadrantes (400m2 c/U)
# en Van Cortlandt Park

y = c(6,0,1,2,1,7,1,5,2,0)

#-----------------------
# Actividad 
# ----------------------
# 1. Calcule la media de densidad de arboles. Construya dos modelos, el primero 
# de ellos calculando la media  considerando considerando cada cuadrante como
# independiente y en el segundo, considerando que todos los cuadrantes pertenecen 
# al mismo parque, y por lo tanto, la densidad en cada cuadrante pueden ser consideradas
# como realizaciones de una unica distribucion.
