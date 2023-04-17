############################################################
###########          Ejemplo Muestreo            ###########
############################################################
########              Basado en                     ########
########      Kery, M., & Royle, J. A. (2016).      ########
########       Applied Hierarchical Modeling        ######## 
########              in Ecology                    ######## 
########                Volume 1                    ########
############################################################

# La siguiente simulación que permite experimentar con la relación entre patrones de 
# puntos y la ocurrencia y abundancia de organismos en un área determinada en 
# relación a los sistemas de muestreos seleccionados para poder hacer el 
# revelamiento de la especie particular.

# Esta función de encuentra dentro de la librería AHMbook
# Antes de poder utilizar esta librería es necesario instalarla

library(AHMbook)
# La función "sim.fn" simula ubicaciones de plantas o animales en grillas de 
# celdas formando un cuadrante de un tamaño específico (quad.size) de acuerdo
# a una distribución de Poisson, donde los organismos están distribuidos al 
# azar dentro del cuadrante. El proceso se caracteriza por una constante "intensity",
# que es el promedio de animales o plantas (puntos) por unidad de área. 

# Parámetros por default en la función
# quad.size=10, cell.size=1, intensity =1
# Si no se especifican alguno de los parámetros de esta función, se tomarán como
# valores los indicados por default
# quad.size = .


# EJEMPLO ----
set.seed(82)
tmp <-sim.fn(quad.size=16, cell.size=2, intensity = 0.5)

# Por los parámetros seleccionados, esperaríamos que haya 128 individuos en el 
# cuadrante completo. Sin embargo, debido a la aleatoriedad del proceso, hay sólo 
# 114.


# ACTIVIDAD ----
# Modifique la función de simulación para responder las siguientes preguntas
# cambiando el tamaño de la celda de muestreo (cell.size)

# 1 - ¿Que observamos en relación a los patrones de abundancia y ocurrencia
# cuando se modifican los tamaños de la celda de muestreo?

# 2 - ¿Qué sucede con los patrones de abundancia y ocurrencia cuando las especies
# son raras? ¿Se mantiene siempre esta relación?


simrep <- 100                 # Run 50 simulation reps
grain <- c(0.1,0.2,0.25,0.5,1,2) # values will be fed into 'cell.size' argument
int <- seq(0.1, 3,,6)         # values will be fed into 'lambda' argument
n.levels <- length(grain)     # number of factor levels in simulation
results <- array(NA, dim = c(n.levels, n.levels, 2, simrep)) # 4-D array !
for(i in 1:n.levels){         # Loop over levels of factor grain
  for(j in 1:n.levels){       # Loop over levels of factor intensity
    for(k in 1:simrep){
      cat("\nDim 1:",i, ", Dim 2:", j, ", Simrep", k)
      tmp <- sim.fn(cell.size = grain[i], intensity = int[j], show.plot = F)
      results[i,j,1:2,k] <- c(mean(tmp$N), tmp$psi)
    }
  }
}


# Plot these two prediction matrices (NOT IN BOOK)
par(mfrow = c(2, 2), mar = c(5,5,2,2), cex.lab = 1.5, cex.axis = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
# Plot mean abundance in sampled quadrats
z1 <- apply(results[,,1,], c(1,2), mean)   # mean abundance
image(x=grain, y=int, z=z1, col = mapPalette(100), axes = T, xlab = "Grain size (cell.size)", ylab = "Intensity of PPP")
contour(x=grain, y=int, z=z1, add = T, col = "blue", labcex = 1.5, lwd = 1.5)
# Plot mean occupancy in sampled quadrats
z2 <- apply(results[,,2,], c(1,2), mean)   # mean occupancy
image(x=grain, y=int, z=z2, col = mapPalette(100), axes = T, xlab = "Grain size (cell.size)", ylab = "Intensity of PPP")
contour(x=grain, y=int, z=z2, add = T, col = "blue", labcex = 1.5, lwd = 1.5)
# Plot relationship between occupancy and abundance for whole range of abundance
plot(results[,,1,], results[,,2,], xlab = "Mean abundance", ylab = "Occupancy", frame = FALSE)
lines(smooth.spline(results[,,2,] ~ results[,,1,], df = 4), lwd = 3, col = "blue")
#abline(0, 1, lwd = 3)
# ... and only for very small abundance
keep <- results[,,1,] < 0.25
plot(results[,,1,][keep], results[,,2,][keep], xlab = "Mean abundance", ylab = "Occupancy", frame = FALSE)
abline(0, 1, lwd = 3)
lines(smooth.spline(results[,,2,][keep] ~ results[,,1,][keep], df = 4), lwd = 3, col = "blue")


sim.fn(quad.size=10, cell.size=1, intensity = 0.5)
sim.fn(quad.size=10, cell.size=2, intensity = 0.5)
sim.fn(quad.size=10, cell.size=5, intensity = 0.5)
sim.fn(quad.size=10, cell.size=10, intensity = 0.5)

sim.fn(quad.size=10, cell.size=1, intensity = 0.1)
sim.fn(quad.size=10, cell.size=1, intensity = 1)
sim.fn(quad.size=10, cell.size=1, intensity = 5)
sim.fn(quad.size=10, cell.size=1, intensity = 10)
