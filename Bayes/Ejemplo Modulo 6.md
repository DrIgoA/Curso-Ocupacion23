---
# Guía de un modelo de comunidades con DA en jags 

```
model {
```
### Previa de la meta comunidad
```
omega ~ dunif(0,1)
```
### Previas de los efectos especie-específicos en la ocupación y detección
#### k va de 1 hasta M (M es el número total de especies observadas y "potenciales"
#### estos hiperparámetros describen a la comunidad donde hay una ordenada (lpsi y lp) de ocupación y detección para cada especie k, y los efectos de las covariables también son modelados para cada especie. Cada uno de esos hiperparámetros, tienen a su vez una media y una precisión
```
for(k in 1:M){
  lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)    
  betalpsi1[k] ~ dnorm(mu.betalpsi1, tau.betalpsi1)
  betalpsi2[k] ~ dnorm(mu.betalpsi2, tau.betalpsi2)
  betalpsi3[k] ~ dnorm(mu.betalpsi3, tau.betalpsi3)
  lp[k] ~ dnorm(mu.lp, tau.lp)
  betalp1[k] ~ dnorm(mu.betalp1, tau.betalp1)
  betalp2[k] ~ dnorm(mu.betalp2, tau.betalp2)
  betalp3[k] ~ dnorm(mu.betalp3, tau.betalp3)
}

# Hyperpriors
# For the model of occupancy
mu.lpsi ~ dnorm(0,0.01)
tau.lpsi <- pow(sd.lpsi, -2)
sd.lpsi ~ dunif(0,8)   # as always, bounds of uniform chosen by trial and error
mu.betalpsi1 ~ dnorm(0,0.1)
tau.betalpsi1 <- pow(sd.betalpsi1, -2)
sd.betalpsi1 ~ dunif(0, 4)
mu.betalpsi2 ~ dnorm(0,0.1)
tau.betalpsi2 <- pow(sd.betalpsi2, -2)
sd.betalpsi2 ~ dunif(0,2)
mu.betalpsi3 ~ dnorm(0,0.1)
tau.betalpsi3 <- pow(sd.betalpsi3, -2)
sd.betalpsi3 ~ dunif(0,2)

# For the model of detection
mu.lp ~ dnorm(0,0.1)
tau.lp <- pow(sd.lp, -2)
sd.lp ~ dunif(0, 2)
mu.betalp1 ~ dnorm(0,0.1)
tau.betalp1 <- pow(sd.betalp1, -2)
sd.betalp1 ~ dunif(0,1)
mu.betalp2 ~ dnorm(0,0.1)
tau.betalp2 <- pow(sd.betalp2, -2)
sd.betalp2 ~ dunif(0,1)
mu.betalp3 ~ dnorm(0,0.1)
tau.betalp3 <- pow(sd.betalp3, -2)
sd.betalp3 ~ dunif(0,1)

# Superpopulation process: Ntotal species sampled out of M available
for(k in 1:M){
   w[k] ~ dbern(omega)
}

# Ecological model for true occurrence (process model)
for(k in 1:M){
  for (i in 1:nsite) {
    logit(psi[i,k]) <- lpsi[k] + betalpsi1[k] * ele[i] + 
      betalpsi2[k] * pow(ele[i],2) + betalpsi3[k] * forest[i]
    mu.psi[i,k] <- w[k] * psi[i,k]
    z[i,k] ~ dbern(mu.psi[i,k])
  }
}

# Observation model for replicated detection/nondetection observations
for(k in 1:M){
  for (i in 1:nsite){
    for(j in 1:nrep){
      logit(p[i,j,k]) <- lp[k] + betalp1[k] * DAT[i,j] + 
        betalp2[k] * pow(DAT[i,j],2) + betalp3[k] * DUR[i,j]
      mu.p[i,j,k] <- z[i,k] * p[i,j,k]
      Y[i,j,k] ~ dbern(mu.p[i,j,k])
    }
  }
}

# Derived quantities
#for(k in 1:M){
#   Nocc.fs[k] <- sum(z[,k])       # Number of occupied sites among the 267
#}
for (i in 1:nsite){
   Nsite[i] <- sum(z[i,])          # Number of occurring species at each site
}
n0 <- sum(w[(nspec+1):(nspec+nz)]) # Number of unseen species
Ntotal <- sum(w[])                 # Total metacommunity size

# Vectors to save (S for save; discard posterior samples for 
# all minus 1 of the potential species to save disk space)
# we do this for nz = 250 (i.e., M = 395)
lpsiS[1:(nspec+1)] <- lpsi[1:(nspec+1)]
betalpsi1S[1:(nspec+1)] <- betalpsi1[1:(nspec+1)]
betalpsi2S[1:(nspec+1)] <- betalpsi2[1:(nspec+1)]
betalpsi3S[1:(nspec+1)] <- betalpsi3[1:(nspec+1)]
lpS[1:(nspec+1)] <- lp[1:(nspec+1)]
betalp1S[1:(nspec+1)] <- betalp1[1:(nspec+1)]
betalp2S[1:(nspec+1)] <- betalp2[1:(nspec+1)]
betalp3S[1:(nspec+1)] <- betalp3[1:(nspec+1)]
}
```