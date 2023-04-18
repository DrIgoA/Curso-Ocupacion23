---
# Guía de un modelo de comunidades con DA en jags 

```
model {
```
Previa de la meta comunidad, es una probabilidad y por ello va de 0 a 1
```
omega ~ dunif(0,1)
```
### Previas de los efectos especie-específicos en la ocupación y detección
loop sobre el subíndice **k** que va de 1 hasta **M** (M es el número total de especies observadas N y "potenciales" nz)
Estos hiperparámetros describen a la comunidad, donde hay una ordenada (lpsi y lp) de ocupación y detección para cada especie k, con efectos de las covariables también para cada especie. Cada uno de esos hiperparámetros, tienen a su vez una media "mu" y una precisión "tau"
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
```
### Hiperprevias de los parámetros de la comunidad. 
Puede haber muchas opciones (normal, uniforme). Eso es bastante subjetivo y muchas veces hay que probar distintas opciones para asegurarse que haya convergencia del modelo. A veces saltan errores y simplemente hayq ue cambiar de distribución o los límites de las mismas.
```
# Modelo de ocupación
mu.lpsi ~ dnorm(0,0.01)
tau.lpsi <- pow(sd.lpsi, -2)
sd.lpsi ~ dunif(0,8)          # como siempre, los límites de la uniforme son elegidos por "prueba y error"
mu.betalpsi1 ~ dnorm(0,0.1)
tau.betalpsi1 <- pow(sd.betalpsi1, -2)
sd.betalpsi1 ~ dunif(0, 4)
mu.betalpsi2 ~ dnorm(0,0.1)
tau.betalpsi2 <- pow(sd.betalpsi2, -2)
sd.betalpsi2 ~ dunif(0,2)
mu.betalpsi3 ~ dnorm(0,0.1)
tau.betalpsi3 <- pow(sd.betalpsi3, -2)
sd.betalpsi3 ~ dunif(0,2)

# Modelo de detección
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

# Proceso se superpoblación: N especies totales muestreadas de M disponibles
for(k in 1:M){
   w[k] ~ dbern(omega)
}
```
### Modelo Ecológico para la ocurrencia verdadera (modelo de proceso). 
Se modela la ocupación **psi** por sitio **i** que va desde 1 hasta número de sitios (nsite), y por especie. Esto se hace dentro de los loops por especie y por sitio.
En este caso las covariables tienen valores por sitio **i** y los parámetros para cada especie **k**
```
for(k in 1:M){
  for (i in 1:nsite) {
    logit(psi[i,k]) <- lpsi[k] + betalpsi1[k] * ele[i] + betalpsi2[k] * pow(ele[i],2) + betalpsi3[k] * forest[i]
    mu.psi[i,k] <- w[k] * psi[i,k]
    z[i,k] ~ dbern(mu.psi[i,k])
  }
}
```
### Modelo de observación para la observaciones replicadas de detección/no detección 
Se modela la detección **p** por sitio **i** que va desde 1 hasta número de sitios (nsite), por especie, y por repetición **j** que va desde 1 al número de ocasiones (nrep), que también pueden modelarse como desbalanceadas. Esto se hace dentro de los loops por especie, por sitio, y por repetición.
En este caso las covariables tienen valores por sitio **i** y repetición **j** y los parámetros para cada especie **k**
La probabilidad de detección es condicional a la ocurrencia verdadera z
```
for(k in 1:M){
  for (i in 1:nsite){
    for(j in 1:nrep){
      logit(p[i,j,k]) <- lp[k] + betalp1[k] * DAT[i,j] + betalp2[k] * pow(DAT[i,j],2) + betalp3[k] * DUR[i,j]
      mu.p[i,j,k] <- z[i,k] * p[i,j,k]
      Y[i,j,k] ~ dbern(mu.p[i,j,k])
    }
  }
}

````
Cantidades derivadas. Podemos estimar el número total de sitios ocupados por cada especie, y el número de especies por sitio (riqueza)
```
for(k in 1:M){
   Nocc.fs[k] <- sum(z[,k])       
}
for (i in 1:nsite){
   Nsite[i] <- sum(z[i,])         
}
n0 <- sum(w[(nspec+1):(nspec+nz)])  #Number of unseen species
Ntotal <- sum(w[])                  #Total metacommunity size
```
