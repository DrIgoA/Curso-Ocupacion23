########################################
### Analisis de los datos 
########################################

# Especificar modelo en lenguaje JAGS
sink("model-test.jags")
cat("
model {

# Priors
alpha.occ ~ dnorm(0, 0.001)
beta.occ ~ dunif(-10, 10)
alpha.p ~ dunif(-10, 10)

# Verosimilitud
for (i in 1:R) {
   # Modelo del Estado Verdadero para el estado verdadero parcialmente observado
   z[i] ~ dbern(psi[i])             # Ocupacion Real z en el sitio i
   logit(psi[i]) <- alpha.occ + beta.occ * X[i]

   for (j in 1:T) {
      # Modelo de Observacion para las observaciones actuales
      y[i,j] ~ dbern(p.eff[i,j])    # Deteccion-no deteccion en i  j
      p.eff[i,j] <- z[i] * p[i,j]
      logit(p[i,j]) <- alpha.p 
      } #j
   } #i

# Cantidades Derivadas
occ.fs <- sum(z[])       # Numero de sitios ocupados entre los estudiados
}
",fill = TRUE)
sink()

# Juntar datos
win.data <- list(y = sodata$y, X = sodata$X, R = nrow(sodata$y), T = ncol(sodata$y))

# Valores iniciales
zst <- apply(sodata$y, 1, max)   #Valores iniciales buenos son esenciales
inits <- function(){list(z = zst, alpha.occ = runif(1, -3, 3), beta.occ = runif(1, -3, 3), alpha.p = runif(1, -3, 3))}

# Parameteros monitoreados
params <- c("alpha.occ", "beta.occ", "alpha.p", "occ.fs")

# Seteo de MCMC 
ni <- 10000
nt <- 8
nb <- 5000
nc <- 3

# Llamar JAGS desde R (BRT 1 min)
out1 <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

naive.pred <- plogis(predict(glm(apply(sodata$y, 1, max) ~ X + I(X^2), family = binomial, data = sodata)))
lin.pred2 <- out1$mean$alpha.occ + out1$mean$beta.occ * sodata$X

plot(sodata$X, sodata$psi, ylim = c(0, 1), main = "Real (rojo), p<1 (negro), \n p=1 (azul)", ylab = "Probabilidad de Ocupacion", xlab = "Covariable", type = "l", lwd = 3, col = "red", las = 1, frame.plot = FALSE)
lines(sodata$X, naive.pred, ylim = c(0 ,1), type = "l", lty = 2, lwd = 3, col = "blue")
lines(sodata$X, plogis(lin.pred2), ylim = c(0, 1), type = "l", lty = 1, lwd = 2, col = "black")
 











########################################
### Analisis de los datos 
########################################

# Especificar modelo en lenguaje JAGS
sink("model-test.jags")
cat("
model {

# Priors
alpha.occ ~ dnorm(0, 0.001)
beta.occ ~ dunif(-10, 10)
alpha.p ~ dunif(-10, 10)

# Verosimilitud
for (i in 1:R) {
   # Modelo del Estado Verdadero para el estado verdadero parcialmente observado
   z[i] ~ dbern(psi[i])             # Ocupacion Real z en el sitio i
   logit(psi[i]) <- alpha.occ + beta.occ * X[i]

   for (j in 1:T) {
      # Modelo de Observacion para las observaciones actuales
      y[i,j] ~ dbern(p.eff[i,j])    # Deteccion-no deteccion en i  j
      p.eff[i,j] <- z[i] * p[i,j]
      logit(p[i,j]) <- alpha.p 
      } #j
   } #i

# Cantidades Derivadas
occ.fs <- sum(z[])       # Numero de sitios ocupados entre los estudiados
}
",fill = TRUE)
sink()

# Juntar datos
win.data <- list(y = sodata$y, X = sodata$X, R = nrow(sodata$y), T = ncol(sodata$y))

# Valores iniciales
zst <- apply(sodata$y, 1, max)   #Valores iniciales buenos son esenciales
inits <- function(){list(z = zst, alpha.occ = runif(1, -3, 3), beta.occ = runif(1, -3, 3), alpha.p = runif(1, -3, 3))}

# Parameteros monitoreados
params <- c("alpha.occ", "beta.occ", "alpha.p", "occ.fs")

# Seteo de MCMC 
ni <- 10000
nt <- 8
nb <- 5000
nc <- 3

# Llamar JAGS desde R (BRT 1 min)
out1 <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

naive.pred <- plogis(predict(glm(apply(sodata$y, 1, max) ~ X + I(X^2), family = binomial, data = sodata)))
lin.pred2 <- out1$mean$alpha.occ + out1$mean$beta.occ * sodata$X

plot(sodata$X, sodata$psi, ylim = c(0, 1), main = "Real (rojo), p<1 (negro), \n p=1 (azul)", ylab = "Probabilidad de Ocupacion", xlab = "Covariable", type = "l", lwd = 3, col = "red", las = 1, frame.plot = FALSE)
lines(sodata$X, naive.pred, ylim = c(0 ,1), type = "l", lty = 2, lwd = 3, col = "blue")
lines(sodata$X, plogis(lin.pred2), ylim = c(0, 1), type = "l", lty = 1, lwd = 2, col = "black")
 