# 4.3.2. Analysis of real data set
# Read in the tit data and have a look at them

# indicar directorio de trabajo
setwd("C:\\Users\\andrea\\Documents\\GitHub\\Curso-Ocupacion23\\Bayes")

tits <- read.table("tits.txt", header = TRUE)
str(tits)
C <- as.matrix(tits[5:13])
obs <- as.matrix(tits[14:22])
first <- as.matrix(tits[23:31])

matplot(1999:2007, t(C), type = "l", lty = 1, lwd = 2, main = "", las = 1, ylab = "Territory counts", xlab = "Year", ylim = c(0, 80), frame = FALSE)

table(obs)
length(table(obs))

apply(first, 2, sum, na.rm = TRUE)

a <- as.numeric(levels(factor(obs)))     # All the levels, numeric
newobs <- obs                            # Gets ObsID from 1:271
for (j in 1:length(a)){newobs[which(obs==a[j])] <- j }
table(newobs)

newobs[is.na(newobs)] <- 272
table(newobs)
first[is.na(first)] <- 0
table(first)

#  (a) Null or intercept-only model
# Specify model in BUGS language
sink("GLM0.jags")
cat("
model {

# Prior
alpha ~ dnorm(0, 0.01)    # log(mean count)

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

# Initial values
inits <- function() list(alpha = runif(1, -10, 10))

# Parameters monitored
params <- c("alpha")

# MCMC settings
ni <- 1200
nt <- 2
nb <- 200
nc <- 3

# Call JAGS from R (BRT < 1 min)
out0 <- jags(win.data, inits, params, "GLM0.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# Summarize posteriors
print(out0, dig = 3)


#  (b) Fixed site effects
# Specify model in BUGS language
sink("GLM1.jags")
cat("
model {

# Priors
for (j in 1:nsite){
   alpha[j] ~ dnorm(0, 0.01)     # Site effects
   }

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha[j]
      }  #j
   }  #i
} 
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

# Initial values (not required for all)
inits <- function() list(alpha = runif(235, -1, 1))

# Parameters monitored
params <- c("alpha")

# MCMC settings
ni <- 1200
nt <- 2
nb <- 200
nc <- 3

# Call JAGS from R (BRT < 1 min)
out1 <- jags(win.data, inits, params, "GLM1.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# Summarize posteriors
print(out1, dig = 2)


#  (c) Fixed site and fixed year effects
# Specify model in BUGS language
sink("GLM2.jags")
cat("
model {

# Priors
for (j in 1:nsite){           # site effects
   alpha[j] ~ dnorm(0, 0.01)
   }

for (i in 2:nyear){           # nyear-1 year effects
   eps[i] ~ dnorm(0, 0.01)
   }
eps[1] <- 0                   # Aliased

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha[j] + eps[i]
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

# Initial values
inits <- function() list(alpha = runif(235, -1, 1), eps = c(NA, runif(8, -1, 1)))

# Parameters monitored
params <- c("alpha", "eps")

# MCMC settings
ni <- 1200
nt <- 2
nb <- 200
nc <- 3

# Call JAGS from R (BRT < 1 min)
out2 <- jags(win.data, inits, params, "GLM2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out2, dig = 2)


#  (d) Random site effects (no year effects)
# Specify model in BUGS language
sink("GLMM1.jags")
cat("
model {

# Priors
for (j in 1:nsite){
   alpha[j] ~ dnorm(mu.alpha, tau.alpha)   # Random site effects
   }
mu.alpha ~ dnorm(0, 0.01)
tau.alpha <- 1/ (sd.alpha * sd.alpha)
sd.alpha ~ dunif(0, 5)

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha[j]
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

# Initial values
inits <- function() list(mu.alpha = runif(1, 2, 3))

# Parameters monitored
params <- c("alpha", "mu.alpha", "sd.alpha")

# MCMC settings
ni <- 1200
nt <- 2
nb <- 200
nc <- 3

# Call JAGS from R (BRT < 1 min)
out3 <- jags(win.data, inits, params, "GLMM1.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out3, dig = 2)


#  (e) Random site and random year effects
# Specify model in BUGS language
sink("GLMM2.jags")
cat("
model {

# Priors
mu ~ dnorm(0, 0.01)                  # Grand mean

for (j in 1:nsite){
   alpha[j] ~ dnorm(0, tau.alpha)    # Random site effects
   }
tau.alpha <- 1/ (sd.alpha * sd.alpha)
sd.alpha ~ dunif(0, 5)

for (i in 1:nyear){
   eps[i] ~ dnorm(0, tau.eps)        # Random year effects
   }
tau.eps <- 1/ (sd.eps * sd.eps)
sd.eps ~ dunif(0, 3)

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
       C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + alpha[j] + eps[i]
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C))

# Initial values (not required for all)
inits <- function() list(mu = runif(1, 0, 4), alpha = runif(235, -2, 2), eps = runif(9, -1, 1))

# Parameters monitored
params <- c("mu", "alpha", "eps", "sd.alpha", "sd.eps")

# MCMC settings
ni <- 6000
nt <- 5
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
out4 <- jags(win.data, inits, params, "GLMM2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out4, dig = 2)


# (f) Random site and random year effects and first-year fixed observer effect
# Specify model in BUGS language
sink("GLMM3.jags")
cat("
model {

# Priors
mu ~ dnorm(0, 0.01)                 # Overall mean
beta2 ~ dnorm(0, 0.01)              # First-year observer effect

for (j in 1:nsite){
   alpha[j] ~ dnorm(0, tau.alpha)   # Random site effects
   }
tau.alpha <- 1/ (sd.alpha * sd.alpha)
sd.alpha ~ dunif(0, 5)

for (i in 1:nyear){
   eps[i] ~ dnorm(0, tau.eps)      # Random year effects
   }
tau.eps <- 1/ (sd.eps * sd.eps)
sd.eps ~ dunif(0, 5)

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + beta2 * first[i,j] + alpha[j] + eps[i]
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C), first = t(first))

# Initial values
inits <- function() list(mu = runif(1, 0, 4), beta2 = runif(1, -1, 1), alpha = runif(235, -2, 2), eps = runif(9, -1, 1))

# Parameters monitored
params <- c("mu", "beta2", "alpha", "eps", "sd.alpha", "sd.eps")

# MCMC settings
ni <- 6000
nt <- 5
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 3 min)
out5 <- jags(win.data, inits, params, "GLMM3.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out5, dig = 2)


#  (g) Random site and random year effects, first-year fixed observer effect and overall linear time trend
# Specify model in BUGS language
sink("GLMM4.jags")
cat("
model {

# Priors
mu ~ dnorm(0, 0.01)                  # Overall intercept
beta1 ~ dnorm(0, 0.01)               # Overall trend 
beta2 ~ dnorm(0, 0.01)               # First-year observer effect

for (j in 1:nsite){
   alpha[j] ~ dnorm(0, tau.alpha)    # Random site effects
   }
tau.alpha <- 1/ (sd.alpha * sd.alpha)
sd.alpha ~ dunif(0, 5)

for (i in 1:nyear){
   eps[i] ~ dnorm(0, tau.eps)        # Random year effects
   }
tau.eps <- 1/ (sd.eps * sd.eps)
sd.eps ~ dunif(0, 3)

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + beta1 * year[i] + beta2 * first[i,j] + alpha[j] + eps[i]
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C), first = t(first), year = ((1:9)-5) / 4)

# Initial values
inits <- function() list(mu = runif(1, 0, 4), beta1 = runif(1, -1, 1), beta2 = runif(1, -1, 1), alpha = runif(235, -2, 2), eps = runif(9, -1, 1))

# Parameters monitored
params <- c("mu", "beta1", "beta2", "alpha", "eps", "sd.alpha", "sd.eps")

# MCMC settings
ni <- 12000
nt <- 6
nb <- 6000
nc <- 3

# Call JAGS from R (BRT 7 min)
out6 <- jags(win.data, inits, params, "GLMM4.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out6, dig = 2)


# (h) The full model 
# Specify model in BUGS language
sink("GLMM5.jags")
cat("
model {

# Priors
mu ~ dnorm(0, 0.01)                  # Overall intercept
beta1 ~ dnorm(0, 0.01)               # Overall trend 
beta2 ~ dnorm(0, 0.01)               # First-year observer effect

for (j in 1:nsite){
   alpha[j] ~ dnorm(0, tau.alpha)    # Random site effects
   }
tau.alpha <- 1/ (sd.alpha * sd.alpha)
sd.alpha ~ dunif(0, 3)

for (i in 1:nyear){
   eps[i] ~ dnorm(0, tau.eps)        # Random year effects
   }
tau.eps <- 1/ (sd.eps * sd.eps)
sd.eps ~ dunif(0, 1)

for (k in 1:nobs){
   gamma[k] ~ dnorm(0, tau.gamma)   # Random observer effects
   }
tau.gamma <- 1/ (sd.gamma * sd.gamma)
sd.gamma ~ dunif(0, 1)

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- mu + beta1 * year[i] + beta2 * first[i,j] + alpha[j] + gamma[newobs[i,j]] + eps[i]
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(C), nsite = nrow(C), nyear = ncol(C), nobs = 272, newobs = t(newobs), first = t(first), year = ((1:9)-5) / 4)

# Initial values
inits <- function() list(mu = runif(1, 0, 4), beta1 = runif(1, -1, 1), beta2 = runif(1, -1, 1), alpha = runif(235, -1, 1), gamma = runif(272, -1, 1), eps = runif(9, -1, 1))

# Parameters monitored
params <- c("mu", "beta1", "beta2", "alpha", "gamma", "eps", "sd.alpha", "sd.gamma", "sd.eps")

# MCMC settings
ni <- 12000
nt <- 6
nb <- 6000
nc <- 3

# Call JAGS from R (BRT 11 min)
out7 <- jags(win.data, inits, params, "GLMM5.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out7, dig = 2)
