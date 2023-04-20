##########################################################
##### CURSO Modelado y estimación de ocupación para  #####
#####  poblaciones y comunidades de especies bajo    #####
#####           enfoque Bayesiano.                   #####
#######      CCT Mendoza - ABRIL 2023                #####
##########################################################
###########           Ejercicio Media            #########
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

summary(lm(y1000 ~ 1))

# 2) Implemente el modelo con estadistica bayesiana, considerando
#    una previa no informativa con una distribcion uniforme = (0,5000)


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
# 6. Compara graficamente las distribuciones previa y poseterior.
#     