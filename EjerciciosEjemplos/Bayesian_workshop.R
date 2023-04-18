##########################################################################
#
# Bayesian Population Analysis using WinBUGS - 
# a hierarchical perspective
#
# Marc Kéry & Michael Schaub
#
#
# WORKSHOP JULY 8-12 2013 , Patuxent, Laurel MD
#
#
##########################################################################


setwd("C:/Users/Andrea/Documents/Modeling/Kery_Schaub_2012")

# 3.3. Poisson GLM in R and WinBUGS for modeling times series of counts
# 3.3.1. Generation and analysis of simulated data
data.fn <- function(n = 40, alpha = 3.5576, beta1 = -0.0912, beta2 = 0.0091, beta3 = -0.00014){
# n: Number of years
# alpha, beta1, beta2, beta3: coefficients of a 
#    cubic polynomial of count on year

# Generate values of time covariate
year <- 1:n

# Signal: Build up systematic part of the GLM
log.expected.count <- alpha + beta1 * year + beta2 * year^2 + beta3 * year^3
expected.count <- exp(log.expected.count)

# Noise: generate random part of the GLM: Poisson noise around expected counts
C <- rpois(n = n, lambda = expected.count)

# Plot simulated data
plot(year, C, type = "b", lwd = 2, col = "black", main = "", las = 1, ylab = "Population size", xlab = "Year", cex.lab = 1.2, cex.axis = 1.2)
lines(year, expected.count, type = "l", lwd = 3, col = "red")

return(list(n = n, alpha = alpha, beta1 = beta1, beta2 = beta2, beta3 = beta3, year = year, expected.count = expected.count, C = C))
}

data <- data.fn()

fm <- glm(C ~ year + I(year^2) + I(year^3), family = poisson, data = data)
summary(fm)


#######################
###NOW DO IT IN WINBUGS

library(R2WinBUGS)

##not necessary when you use the defualt english speakin computer
bugs.dir <-"C:/Program Files/WinBUGS14"  


# Specify model in BUGS language
sink("GLM_Poisson.txt")
cat("
model {

# Priors
alpha ~ dunif(-2000, 2000)
beta1 ~ dunif(-10, 10)
beta2 ~ dunif(-10, 10)
beta3 ~ dunif(-10, 10)

# Likelihood: Note key components of a GLM on one line each
for (i in 1:n){
   C[i] ~ dpois(lambda[i])          # 1. Distribution for random part
   log(lambda[i]) <- log.lambda[i]  # 2. Link function
   log.lambda[i] <- alpha + beta1 * year[i] + beta2 * pow(year[i],2) + beta3 * pow(year[i],3)                      # 3. Linear predictor
   } #i
}
",fill = TRUE)
sink()

# Bundle data
mean.year <- mean(data$year)             # Mean of year covariate
sd.year <- sd(data$year)                 # SD of year covariate
win.data <- list(C = data$C, n = length(data$C), year = (data$year - mean.year) / sd.year)

# Initial values
inits <- function() list(alpha = runif(1, -2, 2), beta1 = runif(1, -3, 3))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "lambda")

# MCMC settings
ni <- 10000
nt <- 2
nb <- 1000
nc <- 3

# Call WinBUGS from R
out <- bugs(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.txt", n.chains = nc, 
n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())


print(out, dig = 3)

plot(1:40, data$C, type = "b", lwd = 2, col = "black", main = "", las = 1, ylab = "Population size", xlab = "Year")
R.predictions <- predict(glm(C ~ year + I(year^2) + I(year^3), family = poisson, data = data), type = "response")
lines(1:40, R.predictions, type = "l", lwd = 3, col = "green")
WinBUGS.predictions <- out$mean$lambda
lines(1:40, WinBUGS.predictions, type = "l", lwd = 3, col = "blue", lty = 2)
cbind(R.predictions, WinBUGS.predictions)

############ JAGS

library(R2jags)

# Specify model in BUGS language
sink("GLM_Poisson.jags")
cat("
model {

# Priors
alpha ~ dunif(-20, 20)
beta1 ~ dunif(-10, 10)
beta2 ~ dunif(-10, 10)
beta3 ~ dunif(-10, 10)

# Likelihood: Note key components of a GLM on one line each
for (i in 1:n){
   C[i] ~ dpois(lambda[i])          # 1. Distribution for random part
   log(lambda[i]) <- log.lambda[i]  # 2. Link function
   log.lambda[i] <- alpha + beta1 * year[i] + beta2 * pow(year[i],2) + beta3 * pow(year[i],3)                      # 3. Linear predictor
   } #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = data$C, n = length(data$C), year = data$year)

# Initial values
inits <- function() list(alpha = runif(1, -2, 2), beta1 = runif(1, -3, 3))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "lambda")

# MCMC settings
ni <- 2000
nt <- 2
nb <- 1000
nc <- 3

# Call JAGS from R
out <- jags(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Bundle data
mean.year <- mean(data$year)             # Mean of year covariate
sd.year <- sd(data$year)                 # SD of year covariate
win.data <- list(C = data$C, n = length(data$C), year = (data$year - mean.year) / sd.year)

# Call JAGS from R (BRT < 1 min)
out <- jags(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(out, dig = 3)

# New MCMC settings with essentially no burnin
ni <- 100
nt <- 1
nb <- 1

# Call JAGS from R (BRT < 1 min)
tmp <- jags(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

plot(1:40, data$C, type = "b", lwd = 2, col = "black", main = "", las = 1, ylab = "Population size", xlab = "Year")
R.predictions <- predict(glm(C ~ year + I(year^2) + I(year^3), family = poisson, data = data), type = "response")
lines(1:40, R.predictions, type = "l", lwd = 3, col = "green")
JAGS.predictions <- out$BUGSoutput$mean$lambda
lines(1:40, JAGS.predictions, type = "l", lwd = 3, col = "blue", lty = 2)
cbind(R.predictions, JAGS.predictions)

# An option in JAGS to see the traceplots of all parameters is the following:
traceplot(tmp)


#################################################################
#################################################################


# 3.5. Binomial GLM for modeling bounded counts or proportions
# 3.5.1. Generation and analysis of simulated data
data.fn <- function(nyears = 40, alpha = 0, beta1 = -0.1, beta2 = -0.9){
# nyears: Number of years
# alpha, beta1, beta2: coefficients

# Generate untransformed and transformed values of time covariate
year <- 1:nyears
YR <- (year-round(nyears/2)) / (nyears / 2)

# Generate values of binomial totals (N)
N <- round(runif(nyears, min = 20, max = 100))

# Signal: build up systematic part of the GLM
exp.p <- plogis(alpha + beta1 * YR + beta2 * (YR^2))

# Noise: generate random part of the GLM: Binomial noise around expected counts (which is N)
C <- rbinom(n = nyears, size = N, prob = exp.p)

# Plot simulated data
plot(year, C/N, type = "b", lwd = 2, col = "black", main = "", las = 1, ylab = "Proportion successful pairs", xlab = "Year", ylim = c(0, 50))
points(year, exp.p, type = "l", lwd = 3, col = "red")

return(list(nyears = nyears, alpha = alpha, beta1 = beta1, beta2 = beta2, year = year, YR = YR, exp.p = exp.p, C = C, N = N))
}

data <- data.fn(nyears = 40, alpha = 1, beta1 = -0.03, beta2 = -0.9)


###################### Binomial model

# Specify model in BUGS language
sink("GLM_Binomial.txt")
cat("
model {

# Priors
alpha ~ dnorm(0, 0.001)
beta1 ~ dnorm(0, 0.001)
beta2 ~ dnorm(0, 0.001)

# Likelihood
for (i in 1:nyears){
   C[i] ~ dbin(p[i], N[i])          # 1. Distribution for random part
   logit(p[i]) <- alpha + beta1 * year[i] + beta2 * pow(year[i],2) # link function and linear predictor
   }
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = data$C, N = data$N, nyears = length(data$C), year = data$YR)

# Initial values
inits <- function() list(alpha = runif(1, -1, 1), beta1 = runif(1, -1, 1), beta2 = runif(1, -1, 1))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "p")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call WinBUGS from R (BRT < 1 min)
out <- bugs(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Binomial.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Plot predictions
WinBUGS.predictions <- out$mean$p
lines(1:length(data$C), WinBUGS.predictions, type = "l", lwd = 3, col = "blue", lty = 2)


###################### Poisson model exercise - transform the binomial into poisson

# Specify model in BUGS language
sink("GLM_Poisson.txt")
cat("
model {

# Priors
alpha ~ dnorm(0, 0.001)
beta1 ~ dnorm(0, 0.001)
beta2 ~ dnorm(0, 0.001)

# Likelihood
for (i in 1:nyears){           # now we have to ignore N
   C[i] ~ dpois(p[i])          # 1. Distribution for random part
   log(p[i]) <- alpha + beta1 * year[i] + beta2 * pow(year[i],2) # link function and linear predictor
   }
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = data$C, nyears = length(data$C), year = data$YR)

# Initial values
inits <- function() list(alpha = runif(1, -1, 1), beta1 = runif(1, -1, 1), beta2 = runif(1, -1, 1))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "p")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call WinBUGS from R (BRT < 1 min)
out <- bugs(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Plot predictions
WinBUGS.predictions <- out$mean$p
lines(1:length(data$C), WinBUGS.predictions, type = "l", lwd = 3, col = "green", lty = 2)


###################### Normal model

# Specify model in BUGS language
sink("GLM_Normal.txt")
cat("
model {

# Priors
alpha ~ dnorm(0, 0.001)
beta1 ~ dnorm(0, 0.001)
beta2 ~ dnorm(0, 0.001)
tau<-1/(sd*sd)
sd~dunif(0,100)

# Likelihood
for (i in 1:nyears){
   C[i] ~ dnorm(p[i], tau)          # 1. Distribution for random part
   p[i] <- alpha + beta1 * year[i] + beta2 * pow(year[i],2) # link function and linear predictor
   }
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = data$C, nyears = length(data$C), year = data$YR)

# Initial values
inits <- function() list(alpha = runif(1, -1, 1), beta1 = runif(1, -1, 1), beta2 = runif(1, -1, 1),tau = runif(1, -1, 1))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "p")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call WinBUGS from R (BRT < 1 min)
out <- bugs(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Normal.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Plot predictions
WinBUGS.predictions <- out$mean$p
lines(1:length(data$C), WinBUGS.predictions, type = "l", lwd = 3, col = "red", lty = 2)

str(out)

##str(outSsims.array)????


############################################################################
#
# 4. Introduction to random effects: Conventional Poisson GLMM for count data
#
##############################################################################


# 4.1. Introduction
# 4.1.1. An example
# Define and plot data
mass <- c(25, 14, 68, 79, 64, 139, 49, 119, 111)
pop <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3))
length <- c(1, 14, 22, 2, 9, 20, 2, 13, 22)
plot(length, mass, col = c(rep("red", 3), rep("blue", 3), rep("green", 3)), 
     xlim = c(-1, 25), ylim = c(0, 140), cex = 1.5, lwd = 2, frame.plot = FALSE, 
     las = 1, pch = 16, xlab = "Length", ylab = "Mass")






# Fit fixed-effects model, print regression parameter estimates and plot regression lines
summary(lm <- lm(mass ~ pop-1 + length))
abline(lm$coef[1], lm$coef[4], col = "red", lwd = 3, lty = 2)
abline(lm$coef[2], lm$coef[4], col = "blue", lwd = 3, lty = 2)
abline(lm$coef[3], lm$coef[4], col = "green", lwd = 3, lty = 2)

library(lme4)

# Fit mixed model, print random effects and plot regression lines
summary(lmm <- lmer(mass ~ length + (1|pop)))  #pop is the random effect 
							     #and is itercept
ranef(lmm)
abline((lmm@fixef[1]+ranef(lmm)$pop)[1,], lmm@fixef[2], col = "red", lwd = 3)
abline((lmm@fixef[1]+ranef(lmm)$pop)[2,], lmm@fixef[2], col = "blue", lwd = 3)
abline((lmm@fixef[1]+ranef(lmm)$pop)[3,], lmm@fixef[2], col = "green", lwd = 3)



# 4.3. Mixed models with random effects for variability among groups (site and year effects)
# 4.3.1. Generation and analysis of simulated data
data.fn <- function(nsite = 5, nyear = 20, alpha = 4.18456, beta1 = 1.90672, beta2 = 0.10852, beta3 = -1.17121, sd.site = 0.5, sd.year = 0.2){
   # nsite: Number of populations
   # nyear: Number of years
   # alpha, beta1, beta2, beta3: cubic polynomial coefficients of year
   # sd.site: standard deviation of the normal distribution assumed for the population intercepts alpha
   # sd.year: standard deviation of the normal distribution assumed for the year effects
   # We standardize the year covariate so that it runs from about 1 to 1

   # Generate data structure to hold counts and log(lambda)
   C <- log.expected.count <- array(NA, dim = c(nyear, nsite))

   # Generate covariate values
   year <- 1:nyear
   yr <- (year-20)/20	# Standardize
   site <- 1:nsite

   # Draw two sets of random effects from their respective distribution
   alpha.site <- rnorm(n = nsite, mean = alpha, sd = sd.site)
   eps.year <- rnorm(n = nyear, mean = 0, sd = sd.year)

   # Loop over populations
   for (j in 1:nsite){
      # Signal (plus first level of noise): build up systematic part of the GLM including random site and year effects
      log.expected.count[,j] <- alpha.site[j] + beta1 * yr + beta2 * yr^2 + beta3 * yr^3 + eps.year
      expected.count <- exp(log.expected.count[,j])

      # Second level of noise: generate random part of the GLM: Poisson noise around expected counts
      C[,j] <- rpois(n = nyear, lambda = expected.count)
      }

   # Plot simulated data
   matplot(year, C, type = "l", lty = 1, lwd = 2, main = "", las = 1, ylab = "Population size", xlab = "Year")

   return(list(nsite = nsite, nyear = nyear, alpha.site = alpha.site, beta1 = beta1, beta2 = beta2, beta3 = beta3, year = year, sd.site = sd.site, sd.year = sd.year, expected.count = expected.count, C = C))
   }

data <- data.fn(nsite = 5, nyear = 20, sd.site = 0.3, sd.year = 0.2)

# Specify model in BUGS language
sink("GLMM_Poisson.txt")
cat("
model {

# Priors
for (j in 1:nsite){
   alpha[j] ~ dnorm(mu, tau.alpha)		# 4. Random site effects
   }
mu ~ dnorm(0, 0.01)				# Hyperparameter 1
tau.alpha <- 1 / (sd.alpha*sd.alpha)	        # Hyperparameter 2
sd.alpha ~ dunif(0, 2)
for (p in 1:3){
   beta[p] ~ dnorm(0, 0.01)
   }

tau.year <- 1 / (sd.year*sd.year)
sd.year ~ dunif(0, 1)				# Hyperparameter 3

# Likelihood
for (i in 1:nyear){
   eps[i] ~ dnorm(0, tau.year)                # 4. Random year effects
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])             # 1. Distribution for random part
      lambda[i,j] <- exp(log.lambda[i,j])     # 2. Link function
      log.lambda[i,j] <- alpha[j] + beta[1] * year[i] + beta[2] * pow(year[i],2) + beta[3] * pow(year[i],3) + eps[i]    # 3. Linear predictor including random site and random year effects
      }  #j
   }  #i
}
",fill = TRUE)
sink()


# Bundle data
win.data <- list(C = data$C, nsite = ncol(data$C), nyear = nrow(data$C), year = (data$year-20) / 20) # Note year standardized

# Initial values
inits <- function() list(mu = runif(1, 0, 2), alpha = runif(data$nsite, -1, 1), beta = runif(3, -1, 1), sd.alpha = runif(1, 0, 0.1), sd.year = runif(1, 0, 0.1))

# Parameters monitored (may want to add "lambda")
params <- c("mu", "alpha", "beta", "sd.alpha", "sd.year")

# MCMC settings (may have to adapt)
ni <- 50000
nt <- 5
nb <- 20000
nc <- 3


# Call WinBUGS from R (BRT 98 min)
out <- bugs(win.data, inits, params, "GLMM_Poisson.txt", n.chains = nc, 
n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)


########################################################################
#
# 5. State-space models
#
###########################################################################

# 5.2. A simple model
n.years <- 25           # Number of years
N1 <- 30                # Initial population size
mean.lambda <- 1.02     # Mean annual population growth rate
sigma2.lambda <- 0.02   # Process (temporal) variation of the growth rate
sigma2.y <- 20          # Variance of the observation error

y <- N <- numeric(n.years)
N[1] <- N1
lambda <- rnorm(n.years-1, mean.lambda, sqrt(sigma2.lambda))  #growth rates simulated
for (t in 1:(n.years-1)){
   N[t+1] <- N[t] * lambda[t]
   }

for (t in 1:n.years){                               #observed data conditonal on the true pop size
   y[t] <- rnorm(1, N[t], sqrt(sigma2.y))
   }

cbind(N,y)

# Specify model in BUGS language
sink("ssm.bug")
cat("
model { 
# Priors and constraints
N.est[1] ~ dunif(0, 500)            # Prior for initial population size
mean.lambda ~ dunif(0, 10)          # Prior for mean growth rate
sigma.proc ~ dunif(0, 10)           # Prior for sd of state process
sigma2.proc <- pow(sigma.proc, 2)
tau.proc <- pow(sigma.proc, -2)
sigma.obs ~ dunif(0, 100)           # Prior for sd of observation process standard deviation
sigma2.obs <- pow(sigma.obs, 2)
tau.obs <- pow(sigma.obs, -2)

# Likelihood
# State process   change in pop from one year to the next
for (t in 1:(T-1)){
   lambda[t] ~ dnorm(mean.lambda, tau.proc)   #lambda comes from normal dist
   N.est[t+1] <- N.est[t] * lambda[t] 
   }
# Observation process
for (t in 1:T) {
   y[t] ~ dnorm(N.est[t], tau.obs)    #the observation depends on the true pop N
   }
}
",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(y = y, T = n.years)

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 5), mean.lambda = runif(1, 0.1, 2), sigma.obs = runif(1, 0, 10), N.est = c(runif(1, 20, 40), rep(NA, (n.years-1))))} 

# Parameters monitored
parameters <- c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT <1 min)
ssm <- bugs(bugs.data, inits, parameters, "ssm.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Define function to draw a graph to summarize results
graph.ssm <- function(ssm, N, y){
   fitted <- lower <- upper <- numeric()
   n.years <- length(y)
   for (i in 1:n.years){
      fitted[i] <- mean(ssm$sims.list$N.est[,i])
      lower[i] <- quantile(ssm$sims.list$N.est[,i], 0.025)
      upper[i] <- quantile(ssm$sims.list$N.est[,i], 0.975)}
   m1 <- min(c(y, fitted, N, lower))
   m2 <- max(c(y, fitted, N, upper))
   par(mar = c(4.5, 4, 1, 1), cex = 1.2)
   plot(0, 0, ylim = c(m1, m2), xlim = c(0.5, n.years), ylab = "Population size", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2, frame = FALSE, axes = FALSE)
   axis(2, las = 1)
   axis(1, at = seq(0, n.years, 5), labels = seq(0, n.years, 5))
   axis(1, at = 0:n.years, labels = rep("", n.years + 1), tcl = -0.25)
   polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
   points(N, type = "l", col = "red", lwd = 2)
   points(y, type = "l", col = "black", lwd = 2)
   points(fitted, type = "l", col = "blue", lwd = 2)
   legend(x = 1, y = m2, legend = c("True", "Observed", "Estimated"), lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("red", "black", "blue"), bty = "n", cex = 1)
}

# Execute function: Produce figure 
graph.ssm(ssm, N, y)

# 5.3. Systematic bias in the observation process
n.years <- 25  # Number of years
N <- rep(50, n.years) 

p <- 0.7
y <- numeric(n.years)
for (t in 1:n.years){
   y[t] <- rbinom(1, N[t], p)
   }
y

# Bundle data
bugs.data <- list(y = y, T = n.years)

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 5), mean.lambda = runif(1, 0.1, 2), sigma.obs = runif(1, 0, 10), N.est = c(runif(1, 30, 60), rep(NA, (n.years-1))))}

# Parameters monitored
parameters <- c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT <1 min)
ssm <- bugs(bugs.data, inits, parameters, "ssm.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(ssm, digits = 3)

# Produce figure
graph.ssm(ssm, N, y)


# 5.3. Systematic bias in the observation process
n.years <- 25  # Number of years
N <- rep(50, n.years)   #we make the pop constant

p <- 0.7    #constant probability of detection
y <- numeric(n.years)
for (t in 1:n.years){
   y[t] <- rbinom(1, N[t], p)
   }
y

# Bundle data
bugs.data <- list(y = y, T = n.years)

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 5), mean.lambda = runif(1, 0.1, 2), sigma.obs = runif(1, 0, 10), N.est = c(runif(1, 30, 60), rep(NA, (n.years-1))))}

# Parameters monitored
parameters <- c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT <1 min)
ssm <- bugs(bugs.data, inits, parameters, "ssm.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(ssm, digits = 3)

# Produce figure
graph.ssm(ssm, N, y)

 n.years <- 25  # Number of years
N <- rep(50, n.years)

lp <- -0.5 + 0.1*(1:n.years)  # Increasing trend of logit p
p <- plogis(lp)
y <- numeric(n.years)
for (t in 1:n.years){
   y[t] <- rbinom(1, N[t], p[t])
   }

# Bundle data
bugs.data <- list(y = y, T = n.years)

# Call WinBUGS from R (BRT <1 min)
ssm <- bugs(bugs.data, inits, parameters, "ssm.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Produce figure
graph.ssm(ssm, N, y)
points(N*p, col = "black", type = "l", lwd = 2, lty = 2)
legend(x = 1, y = 45.5, legend = "Np", lwd = 2, col = "black", lty = 2, bty = "n")


# 5.4. Real example: House martin population counts in the village of Magden
# Specify model in BUGS language
sink("ssm.bug")
cat("
model {
# Priors and constraints
logN.est[1] ~ dnorm(5.6, 0.01)       # Prior for initial population size
mean.r ~ dnorm(1, 0.001)             # Prior for mean growth rate
sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
sigma2.proc <- pow(sigma.proc, 2)
tau.proc <- pow(sigma.proc, -2)
sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
sigma2.obs <- pow(sigma.obs, 2)
tau.obs <- pow(sigma.obs, -2)

# Likelihood
# State process
for (t in 1:(T-1)){
   r[t] ~ dnorm(mean.r, tau.proc)
   logN.est[t+1] <- logN.est[t] + r[t]
   }
# Observation process
for (t in 1:T) {
   y[t] ~ dnorm(logN.est[t], tau.obs)
   }

# Population sizes on real scale
for (t in 1:T) {
   N.est[t] <- exp(logN.est[t])
   }
}
",fill = TRUE)
sink()

# House martin population data from Magden
pyears <- 6 # Number of future years with predictions
hm <- c(271, 261, 309, 318, 231, 216, 208, 226, 195, 226, 233, 209, 226, 192, 191, 225, 245, 205, 191, 174, rep(NA, pyears))
year <- 1990:(2009 + pyears)

# Bundle data
bugs.data <- list(y = log(hm), T = length(year))

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 1), mean.r = rnorm(1), sigma.obs = runif(1, 0, 1), logN.est = c(rnorm(1, 5.6, 0.1), rep(NA, (length(year)-1))))}

# Parameters monitored
parameters <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 200000
nt <- 6
nb <- 100000
nc <- 3

# Call WinBUGS from R (BRT 3 min)
hm.ssm <- bugs(bugs.data, inits, parameters, "ssm.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(hm.ssm, digits = 3)

# Draw figure
fitted <- lower <- upper <- numeric()
year <- 1990:2015
n.years <- length(hm)
for (i in 1:n.years){
   fitted[i] <- mean(hm.ssm$sims.list$N.est[,i])
   lower[i] <- quantile(hm.ssm$sims.list$N.est[,i], 0.025)
   upper[i] <- quantile(hm.ssm$sims.list$N.est[,i], 0.975)}
m1 <- min(c(fitted, hm, lower), na.rm = TRUE)
m2 <- max(c(fitted, hm, upper), na.rm = TRUE)
par(mar = c(4.5, 4, 1, 1))
plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years), ylab = "Population size", xlab = "Year", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
axis(2, las = 1)
axis(1, at = 1:n.years, labels = year)
polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
points(hm, type = "l", col = "black", lwd = 2)
points(fitted, type = "l", col = "blue", lwd = 2)
legend(x = 1, y = 150, legend = c("Counts", "Estimates"), lty = c(1, 1), lwd = c(2, 2), col = c("black", "blue"), bty = "n", cex = 1)

# Probability of N(2015) < N(2009)
mean(hm.ssm$sims.list$N.est[,26] < hm.ssm$mean$N.est[20])

#sims.list are all the random draws
hm.ssm$sims.list$N.est[,26] 


# Exercise 2, using Real example: House martin population counts in the village of Magden
# Specify model in BUGS language

sink("ssm.bug")
cat("
model {
# Priors and constraints
logN.est[1] ~ dnorm(5.6, 0.01)       # Prior for initial population size
mean.r ~ dnorm(1, 0.001)             # Prior for mean growth rate
sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
sigma2.proc <- pow(sigma.proc, 2)
tau.proc <- pow(sigma.proc, -2)

######

sigma ~ dunif(0, 1)             # Prior for sd of state process
sigma2 <- pow(sigma, 2)
tau<-pow(sigma, -2)

#####
sigma.obs1 ~ dunif(0, 1)              # Prior for sd of observation process
sigma2.obs1 <- pow(sigma.obs1, 2)
tau.obs1 <- pow(sigma.obs1, -2)

sigma.obs2 ~ dunif(0, 1)              # Prior for sd of observation process
sigma2.obs2 <- pow(sigma.obs2, 2)
tau.obs2 <- pow(sigma.obs2, -2)

# Likelihood
# State process
for (t in 1:(T-1)){
   r[t] ~ dnorm(mean.r, tau.proc)
   logN.est[t+1] <- logN.est[t] + r[t]
   }
# Observation process
for (t in 1:9) {
   y[t] ~ dnorm(logN.est[t], tau.obs1)
   }
   
for (t in 10:T) {
  y[t] ~ dnorm(logN.est[t], tau.obs2)
   }
# Population sizes on real scale
for (t in 1:T) {
   N.est[t] <- exp(logN.est[t])
   }
}
",fill = TRUE)
sink()

# House martin population data from Magden
pyears <- 6 # Number of future years with predictions
hm <- c(271, 261, 309, 318, 231, 216, 208, 226, 195, 226, 233, 209, 226, 192, 191, 225, 245, 205, 191, 174, rep(NA, pyears))
year <- 1990:(2009 + pyears)

# Bundle data
bugs.data <- list(y = log(hm), T = length(year))

# Initial values
inits <- function(){list(sigma = runif(1, 0, 1),sigma.year1 = runif(1, 0, 1),sigma.year2 = runif(1, 0, 1),sigma.proc = runif(1, 0, 1), mean.r = rnorm(1), sigma.obs = runif(1, 0, 1), logN.est = c(rnorm(1, 5.6, 0.1), rep(NA, (length(year)-1))))}

# Parameters monitored
parameters <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N.est","sigma","sigma.year1")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT 3 min)
hm.ssm <- bugs(bugs.data, inits, parameters, "ssm.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(hm.ssm, digits = 3)


#######################################################################
#
# 6. Estimation of the size of a closed population
# 
#######################################################################


# 6.3. Analysis of a real data set: model Mtbh for species richness estimation
# Read in data and look at them
p610 <- read.table("p610.txt", header = TRUE)
y <- p610[,5:9]                           # Grab counts
y[y > 1] <- 1                             # Counts to det-nondetections
C <- sum(apply(y, 1, max)) ; print(C)     # Number of observed species
table(apply(y, 1, sum))                   # Capture-frequencies

# Specify model in BUGS language
sink("M_tbh.txt")
cat("
model {

# Priors
omega ~ dunif(0, 1)
for  (j in 1:T){
   alpha[j] <- log(mean.p[j] / (1-mean.p[j])) # Define logit 
   mean.p[j] ~ dunif(0, 1) 	# Detection intercepts
   }
gamma ~ dnorm(0, 0.01)
tau <- 1 / (sd * sd)
sd ~ dunif(0, 3)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   eps[i] ~ dnorm(0, tau)I(-16, 16)

   # First occasion: no term for recapture (gamma)
   y[i,1] ~ dbern(p.eff[i,1])
   p.eff[i,1] <- z[i] * p[i,1]
   p[i,1] <- 1 / (1 + exp(-lp[i,1]))
   lp[i,1] <- alpha[1] + eps[i]

   # All subsequent occasions: includes recapture term (gamma)
   for (j in 2:T){
      y[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp[i,j]))   
      lp[i,j] <- alpha[j] + eps[i] + gamma * y[i,(j-1)]
      } #j
   } #i

# Derived quantities
N <- sum(z[])
} 
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = as.matrix(y), M = nrow(y), T = ncol(y))

# Initial values
inits <- function() list(z = rep(1, nrow(y)), sd = runif(1, 0.1, 0.9))

# Parameters monitored
params <- c("N", "mean.p", "gamma", "sd", "omega")

# MCMC settings
ni <- 50000
nt <- 4
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT 24 min)
out <- bugs(win.data, inits, params, "M_tbh.txt", n.chains = nc, 
n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors and plot posterior for N
print(out, dig = 3)
par(mfrow = c(1,2))
hist(out$sims.list$N, breaks = 35, col = "gray", main = "", xlab = "Community size", las = 1, xlim = c(30, 100), freq = FALSE)
abline(v = C, col = "black", lwd = 3)

# Define model
sink("M0.txt")
cat("
model {

# Priors
omega ~ dunif(0, 1)
p ~ dunif(0, 1)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   for (j in 1:T){
      y[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p
      } #j
   } #i

# Derived quantities
N <- sum(z[])
} # end model
",fill = TRUE)
sink()

# Initial values
inits <- function() list(z = rep(1, nrow(y)))

# Define parameters to be monitored
params <- c("N", "p", "omega")

# MCMC settings
ni <- 50000
nt <- 4
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT 1 min)
out0 <- bugs(win.data, inits, params, "M0.txt", n.chains = nc, 
n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = bugs.dir, working.directory = getwd())

# Inspect output
print(out0, dig = 3)


####################################################################
#
# 7. Estimation of survival probabilities using capture-recapture data
#
#####################################################################

# 7.3. Models with constant parameters
# Define parameter values
n.occasions <- 6                   # Number of capture occasions
marked <- rep(50, n.occasions-1)   # Annual number of newly marked individuals
phi <- rep(0.65, n.occasions-1)
p <- rep(0.4, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked){
   n.occasions <- dim(PHI)[2] + 1
   CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
   # Define a vector with the occasion of marking
   mark.occ <- rep(1:length(marked), marked[1:length(marked)])
   # Fill the CH matrix
   for (i in 1:sum(marked)){
      CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
      if (mark.occ[i]==n.occasions) next
      for (t in (mark.occ[i]+1):n.occasions){
         # Bernoulli trial: does individual survive occasion?
         sur <- rbinom(1, 1, PHI[i,t-1])
         if (sur==0) break		# If dead, move to next individual 
         # Bernoulli trial: is individual recaptured? 
         rp <- rbinom(1, 1, P[i,t-1])
         if (rp==1) CH[i,t] <- 1
         } #t
      } #i
   return(CH)
   }

# Execute function
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-c-c.bug")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- mean.phi
      p[i,t] <- mean.p
      } #t
   } #i

mean.phi ~ dunif(0, 1)         # Prior for mean survival
mean.p ~ dunif(0, 1)           # Prior for mean recapture

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2])

# Function to create a matrix of initial values for latent state z
ch.init <- function(ch, f){
   for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
   return(ch)
   }

# Initial values
inits <- function(){list(z = ch.init(CH, f), mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1))}

# Parameters monitored
parameters <- c("mean.phi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 1 min)
cjs.c.c <- bugs(bugs.data, inits, parameters, "cjs-c-c.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(cjs.c.c, digits = 3)


###CODE FOR EXERCISE 1

# Specify model in BUGS language
sink("cjs-c-c.bug")
cat("
model {

mean.phi ~ dunif(0, 1)         # Prior for mean survival
mean.p ~ dunif(0, 1)           # Prior for mean recapture

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- mean.phi* z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- mean.p * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2])

# Function to create a matrix of initial values for latent state z
ch.init <- function(ch, f){
   for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
   return(ch)
   }

# Initial values
inits <- function(){list(z = ch.init(CH, f), mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1))}

# Parameters monitored
parameters <- c("mean.phi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 1 min)
cjs.c.c <- bugs(bugs.data, inits, parameters, "cjs-c-c.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(cjs.c.c, digits = 3)


known.state.cjs <- function(ch){
   state <- ch
   for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
      }
   state[state==0] <- NA
   return(state)
   }




# Function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch,f){
   for (i in 1:dim(ch)[1]){
      if (sum(ch[i,])==1) next
      n2 <- max(which(ch[i,]==1))
      ch[i,f[i]:n2] <- NA
      }
   for (i in 1:dim(ch)[1]){
   ch[i,1:f[i]] <- NA
   }
   return(ch)
   }





# 7.7. Models with age effects
# Define parameter values
n.occasions <- 10                   # Number of capture occasions
marked.j <- rep(200, n.occasions-1) # Annual number of newly marked juveniles
marked.a <- rep(30, n.occasions-1)  # Annual number of newly marked adults
phi.juv <- 0.3                      # Juvenile annual survival
phi.ad <- 0.65                      # Adult annual survival
p <- rep(0.5, n.occasions-1)        # Recapture
phi.j <- c(phi.juv, rep(phi.ad, n.occasions-2))
phi.a <- rep(phi.ad, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI.J <- matrix(0, ncol = n.occasions-1, nrow = sum(marked.j))
for (i in 1:length(marked.j)){
   PHI.J[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:(n.occasions-1)] <- matrix(rep(phi.j[1:(n.occasions-i)],marked.j[i]), ncol = n.occasions-i, byrow = TRUE)
   }
P.J <- matrix(rep(p, sum(marked.j)), ncol = n.occasions-1, nrow = sum(marked.j), byrow = TRUE)
PHI.A <- matrix(rep(phi.a, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)
P.A <- matrix(rep(p, sum(marked.a)), ncol = n.occasions-1, nrow = sum(marked.a), byrow = TRUE)

# Apply simulation function
CH.J <- simul.cjs(PHI.J, P.J, marked.j)
CH.A <- simul.cjs(PHI.A, P.A, marked.a) 

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f.j <- apply(CH.J, 1, get.first)
f.a <- apply(CH.A, 1, get.first)

# Create matrices X indicating age classes
x.j <- matrix(NA, ncol = dim(CH.J)[2]-1, nrow = dim(CH.J)[1])
x.a <- matrix(NA, ncol = dim(CH.A)[2]-1, nrow = dim(CH.A)[1])
for (i in 1:dim(CH.J)[1]){
   for (t in f.j[i]:(dim(CH.J)[2]-1)){
      x.j[i,t] <- 2
      x.j[i,f.j[i]] <- 1   
      } #t
   } #i
for (i in 1:dim(CH.A)[1]){
   for (t in f.a[i]:(dim(CH.A)[2]-1)){
      x.a[i,t] <- 2
      } #t
   } #i

CH <- rbind(CH.J, CH.A)
f <- c(f.j, f.a)
x <- rbind(x.j, x.a)

# Specify model in BUGS language
sink("cjs-age.bug")
cat("
model {
# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- beta[x[i,t]]
      p[i,t] <- mean.p
      } #t
   } #i
for (u in 1:2){
   beta[u] ~ dunif(0, 1)              # Priors for age-specific survival
   }
mean.p ~ dunif(0, 1)                  # Prior for mean recapture
# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), x = x)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), beta = runif(2, 0, 1), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("beta", "mean.p")

# MCMC settings
ni <- 2000
nt <- 3
nb <- 1000
nc <- 3

# Call WinBUGS from R (BRT 3 min)
cjs.age <- bugs(bugs.data, inits, parameters, "cjs-age.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

print(cjs.age, digits = 3)

# Create matrix X indicating age classes
x <- matrix(NA, ncol = dim(CH)[2]-1, nrow = dim(CH)[1])
for (i in 1:dim(CH)[1]){
   for (t in f[i]:(dim(CH)[2]-1)){
      x[i,t] <- t-f[i]+1
      } #t 
   } #i

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- mu + beta*x[i,t]
      p[i,t] <- mean.p
      } #t
   } #i
mu ~ dnorm(0, 0.01)             # Prior for mean of logit survival
beta ~ dnorm(0, 0.01)           # Prior for slope parameter
for (i in 1:(n.occasions-1)){
   phi.age[i] <- 1 / (1+exp(-mu –beta*i))   # Logit back-transformation 
   }
mean.p ~ dunif(0, 1)                # Prior for mean recapture



########### EXERCISE 5  (copy code from website, this one is not finished)
# Create matrix X indicating age classes
x <- matrix(NA, ncol = dim(CH)[2]-1, nrow = dim(CH)[1])
for (i in 1:dim(CH)[1]){
   for (t in f[i]:(dim(CH)[2]-1)){
      x[i,t] <- t-f[i]+1
      } #t 
   } #i

# Specify model in BUGS language
sink("cjs-age.bug")
cat("
model {
# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- beta[x[i,t]]  #add logit, first beta separate
      p[i,t] <- mean.p
      } #t
   } #i
beta[1]~dnorm(0,0.001)
for (u in 2:n.ocassions-1){
   beta[u] ~ mu + gamma*(u-1)            #linear model for age
   }
#######################
mu~dnorm()
gamma~
#######################
mean.p ~ dunif(0, 1)                  # Prior for mean recapture
# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), x = x)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), beta = runif(2, 0, 1), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("beta", "mean.p")

# MCMC settings
ni <- 2000
nt <- 3
nb <- 1000
nc <- 3

# Call WinBUGS from R (BRT 3 min)
cjs.age <- bugs(bugs.data, inits, parameters, "cjs-age.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

print(cjs.age, digits = 3)


############################################################
#
# 9. Multistate capture-recapture models
#
##############################################################

# 9.2. Estimation of movement between two sites
# 9.2.1. Model description
# 9.2.2. Generation of simulated data
# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phiA <- 0.8
phiB <- 0.7
psiAB <- 0.3
psiBA <- 0.5
pA <- 0.7
pB <- 0.4
n.occasions <- 6
n.states <- 3
n.obs <- 3
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions)  
marked[,2] <- rep(60, n.occasions)
marked[,3] <- rep(0, n.occasions)

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
   # Dimension 1: state of departure
   # Dimension 2: state of arrival
   # Dimension 3: individual
   # Dimension 4: time

# 1. State process matrix
totrel <- sum(marked)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.STATE[,,i,t] <- matrix(c(
      phiA*(1-psiAB), phiA*psiAB,     1-phiA,
      phiB*psiBA,     phiB*(1-psiBA), 1-phiB,
      0,              0,              1       ), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.OBS[,,i,t] <- matrix(c(
      pA, 0,  1-pA,
      0,  pB, 1-pB,
      0,  0,  1       ), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# Define function to simulate multistate capture-recapture data
simul.ms <- function(PSI.STATE, PSI.OBS, marked, unobservable = NA){
   # Unobservable: number of state that is unobservable
   n.occasions <- dim(PSI.STATE)[4] + 1
   CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
   # Define a vector with the occasion of marking
   mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked))
   g <- colSums(marked)
   for (s in 1:dim(PSI.STATE)[1]){
      if (g[s]==0) next  # To avoid error message if nothing to replace
      mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <-
      rep(1:n.occasions, marked[1:n.occasions,s])
      } #s
   for (i in 1:sum(marked)){
      for (s in 1:dim(PSI.STATE)[1]){
         if (mark.occ[i,s]==0) next
         first <- mark.occ[i,s]
         CH[i,first] <- s
         CH.TRUE[i,first] <- s
         } #s
      for (t in (first+1):n.occasions){
         # Multinomial trials for state transitions
         if (first==n.occasions) next
         state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
         CH.TRUE[i,t] <- state
         # Multinomial trials for observation process
         event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
         CH[i,t] <- event
         } #t
      } #i
   # Replace the NA and the highest state number (dead) in the file by 0
   CH[is.na(CH)] <- 0
   CH[CH==dim(PSI.STATE)[1]] <- 0
   CH[CH==unobservable] <- 0
   id <- numeric(0)
   for (i in 1:dim(CH)[1]){
      z <- min(which(CH[i,]!=0))
      ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id))
      }
   return(list(CH=CH[-id,], CH.TRUE=CH.TRUE[-id,]))
   # CH: capture histories to be used
   # CH.TRUE: capture histories with perfect observation
   }

# Execute function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in A, 2 = seen alive in B, 3 = not seen
rCH <- CH          # Recoded CH
rCH[rCH==0] <- 3



# 9.2.3. Analysis of the model
# Specify model in BUGS language
sink("ms.bug")
cat("
model {

# -------------------------------------------------
# Parameters:
# phiA: survival probability at site A
# phiB: survival probability at site B
# psiAB: movement probability from site A to site B
# psiBA: movement probability from site B to site A
# pA: recapture probability at site A
# pB: recapture probability at site B
# -------------------------------------------------
# States (S):
# 1 alive at A
# 2 alive at B
# 3 dead
# Observations (O):  
# 1 seen at A 
# 2 seen at B
# 3 not seen
# -------------------------------------------------

# Priors and constraints
for (t in 1:(n.occasions-1)){
   phiA[t] <- mean.phi[1]
   phiB[t] <- mean.phi[2]
   psiAB[t] <- mean.psi[1]
   psiBA[t] <- mean.psi[2]
   pA[t] <- mean.p[1]
   pB[t] <- mean.p[2]
   }
for (u in 1:2){
   mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
   mean.psi[u] ~ dunif(0, 1)    # Priors for mean transitions
   mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
   }

# Define state-transition and observation matrices
for (i in 1:nind){  
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phiA[t] * (1-psiAB[t])
      ps[1,i,t,2] <- phiA[t] * psiAB[t]
      ps[1,i,t,3] <- 1-phiA[t]
      ps[2,i,t,1] <- phiB[t] * psiBA[t]
      ps[2,i,t,2] <- phiB[t] * (1-psiBA[t])
      ps[2,i,t,3] <- 1-phiB[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      
      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- pA[t]
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 1-pA[t]
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- pB[t]
      po[2,i,t,3] <- 1-pB[t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 1
      } #t
   } #i

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
   } #i
}
",fill = TRUE)
sink()

# Function to create known latent states z
known.state.ms <- function(ms, notseen){
   # notseen: label for ‘not seen’
   state <- ms
   state[state==notseen] <- NA
   for (i in 1:dim(ms)[1]){
      m <- min(which(!is.na(state[i,])))
      state[i,m] <- NA
      }
   return(state)
   }

# Function to create initial values for unknown z
ms.init.z <- function(ch, f){
   for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
   states <- max(ch, na.rm = TRUE)
   known.states <- 1:(states-1)
   v <- which(ch==states)
   ch[-v] <- NA
   ch[v] <- sample(known.states, length(v), replace = TRUE)
   return(ch)
   }

# Bundle data
bugs.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 3))

# Initial values
inits <- function(){list(mean.phi = runif(2, 0, 1), mean.psi = runif(2, 0, 1), mean.p = runif(2, 0, 1), z = ms.init.z(rCH, f))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.psi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 2000
nc <- 3

# Call WinBUGS from R (BRT 8 min)
ms <- bugs(bugs.data, inits, parameters, "ms.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
print(ms, digits = 3)




############################################################################
#
# 11. Integrated population models
# 
##############################################################################

# 11.3. Example of a simple IPM (counts, capture-recapture, reproduction)
# 11.3.1. Load data
# Population counts (from years 1 to 10)
y <- c(45, 48, 44, 59, 62, 62, 55, 51, 46, 42)

# Capture-recapture data (in m-array format, from years 1 to 10)
m <- matrix(c(11,  0,  0,  0,  0,  0,  0,  0,  0,  70,
               0, 12,  0,  1,  0,  0,  0,  0,  0,  52,
               0,  0, 15,  5,  1,  0,  0,  0,  0,  42,
               0,  0,  0,  8,  3,  0,  0,  0,  0,  51,
               0,  0,  0,  0,  4,  3,  0,  0,  0,  61,
               0,  0,  0,  0,  0, 12,  2,  3,  0,  66,
               0,  0,  0,  0,  0,  0, 16,  5,  0,  44,
               0,  0,  0,  0,  0,  0,  0, 12,  0,  46,
               0,  0,  0,  0,  0,  0,  0,  0, 11,  71,
              10,  2,  0,  0,  0,  0,  0,  0,  0,  13,
               0,  7,  0,  1,  0,  0,  0,  0,  0,  27,
               0,  0, 13,  2,  1,  1,  0,  0,  0,  14,
               0,  0,  0, 12,  2,  0,  0,  0,  0,  20,
               0,  0,  0,  0, 10,  2,  0,  0,  0,  21,
               0,  0,  0,  0,  0, 11,  2,  1,  1,  14,
               0,  0,  0,  0,  0,  0, 12,  0,  0,  18,
               0,  0,  0,  0,  0,  0,  0, 11,  1,  21,
               0,  0,  0,  0,  0,  0,  0,  0, 10,  26), ncol = 10, byrow = TRUE)

# Productivity data (from years 1 to 9)
J <- c(64, 132,  86, 154, 156, 134, 116, 106, 110)
R <- c(21, 28, 26, 38, 35, 33, 31, 30, 33) 


# 11.3.2. Analysis of the model
# Specify model in BUGS language
sink("ipm.bug")
cat("
model {
#-------------------------------------------------
#  Integrated population model
#  - Age structured model with 2 age classes: 
#		1-year old and adult (at least 2 years old)
#  - Age at first breeding = 1 year
#  - Prebreeding census, female-based
#  - All vital rates assumed to be constant
#-------------------------------------------------

#-------------------------------------------------
# 1. Define the priors for the parameters
#-------------------------------------------------
# Observation error
tauy <- pow(sigma.y, -2)
sigma.y ~ dunif(0, 50)
sigma2.y <- pow(sigma.y, 2)

# Initial population sizes
N1[1] ~ dnorm(100, 0.0001)I(0,)     # 1-year
Nad[1] ~ dnorm(100, 0.0001)I(0,)    # Adults

# Survival and recapture probabilities, as well as productivity
for (t in 1:(nyears-1)){
   sjuv[t] <- mean.sjuv
   sad[t] <- mean.sad
   p[t] <- mean.p
   f[t] <- mean.fec
   }

mean.sjuv ~ dunif(0, 1)
mean.sad ~ dunif(0, 1)
mean.p ~ dunif(0, 1)
mean.fec ~ dunif(0, 20)

#-------------------------------------------------
# 2. Derived parameters
#-------------------------------------------------
# Population growth rate
for (t in 1:(nyears-1)){
   lambda[t] <- Ntot[t+1] / Ntot[t]
   }

#-------------------------------------------------
# 3. The likelihoods of the single data sets
#-------------------------------------------------
# 3.1. Likelihood for population population count data (state-space model)
   # 3.1.1 System process
   for (t in 2:nyears){
      mean1[t] <- f[t-1] / 2 * sjuv[t-1] * Ntot[t-1]
      N1[t] ~ dpois(mean1[t])
      Nad[t] ~ dbin(sad[t-1], Ntot[t-1])
      }
   for (t in 1:nyears){
      Ntot[t] <- Nad[t] + N1[t]
      }
   
   # 3.1.2 Observation process
   for (t in 1:nyears){
      y[t] ~ dnorm(Ntot[t], tauy)
      }

# 3.2 Likelihood for capture-recapture data: CJS model (2 age classes)
# Multinomial likelihood
for (t in 1:2*(nyears-1)){
   m[t,1:nyears] ~ dmulti(pr[t,], r[t])
   }

# Calculate the number of released individuals
for (t in 1:2*(nyears-1)){
   r[t] <- sum(m[t,])
   }

# m-array cell probabilities for juveniles
for (t in 1:(nyears-1)){
   # Main diagonal
   q[t] <- 1-p[t]
   pr[t,t] <- sjuv[t] * p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr[t,j] <- sjuv[t]*prod(sad[(t+1):j])*prod(q[t:(j-1)])*p[j]
      } #j	
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j] <- 0
      } #j
   # Last column: probability of non-recapture
   pr[t,nyears] <- 1-sum(pr[t,1:(nyears-1)])
   } #t

# m-array cell probabilities for adults
for (t in 1:(nyears-1)){
   # Main diagonal
   pr[t+nyears-1,t] <- sad[t] * p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr[t+nyears-1,j] <- prod(sad[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t+nyears-1,j] <- 0
      } #j
   # Last column
   pr[t+nyears-1,nyears] <- 1 - sum(pr[t+nyears-1,1:(nyears-1)])
   } #t

# 3.3. Likelihood for productivity data: Poisson regression
for (t in 1:(nyears-1)){
   J[t] ~ dpois(rho[t])
   rho[t] <- R[t]*f[t]
   }
}
",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(m = m, y = y, J = J, R = R, nyears = dim(m)[2])

# Initial values
inits <- function(){list(mean.sjuv = runif(1, 0, 1), mean.sad = runif(1, 0, 1), mean.p = runif(1, 0, 1), mean.fec = runif(1, 0, 10), N1 = rpois(dim(m)[2], 30), Nad = rpois(dim(m)[2], 30), sigma.y = runif(1 ,0, 10))}  

# Parameters monitored
parameters <- c("mean.sjuv", "mean.sad", "mean.p", "mean.fec", "N1", "Nad", "Ntot", "lambda", "sigma2.y")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 2 min)
ipm <- bugs(bugs.data, inits, parameters, "ipm.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

print(ipm, digits = 3)

# Produce Fig. 11-4
par(cex = 1.2)
lower <- upper <- numeric()
for (i in 1:10){
   lower[i] <- quantile(ipm$sims.list$Ntot[,i], 0.025)
   upper[i] <- quantile(ipm$sims.list$Ntot[,i], 0.975)
   }
plot(ipm$mean$Ntot, type = "b", ylim = c(35, 65), ylab = "Population size", xlab = "Year", las = 1, pch = 16, col = "blue", frame = F, cex = 1.5)
segments(1:10, lower, 1:10, upper, col = "blue")
points(y, type = "b", col = "black", pch = 16, lty = 2, cex = 1.5)
legend(x = 1, y = 65, legend = c("Counts", "Estimates"), pch = c(16, 16), col = c("black", "blue"), lty = c(2, 1), bty = "n")
 

# 11.4. Another example of an IPM: Estimating productivity without explicit productivity data
# Specify model in BUGS language
sink("ipm-prod.bug")
cat("
model {
#-------------------------------------------------
#  Integrated population model
#  - Age structured model with 2 age classes: 
#		1-year old and adult (at least 2 years old)
#  - Age at first breeding = 1 year
#  - Prebreeding census, female-based
#  - All vital rates assumed to be constant
#-------------------------------------------------

#-------------------------------------------------
# 1. Define the priors for the parameters
#-------------------------------------------------
# Observation error
tauy <- pow(sigma.y, -2)
sigma.y ~ dunif(0, 50)
sigma2.y <- pow(sigma.y, 2)

# Initial population sizes
N1[1] ~ dnorm(100, 0.0001)I(0,)     # 1-year
Nad[1] ~ dnorm(100, 0.0001)I(0,)    # Adults

# Survival and recapture probabilities, as well as productivity
for (t in 1:(nyears-1)){
   sjuv[t] <- mean.sjuv
   sad[t] <- mean.sad
   p[t] <- mean.p
   f[t] <- mean.fec
   }

mean.sjuv ~ dunif(0, 1)
mean.sad ~ dunif(0, 1)
mean.p ~ dunif(0, 1)
mean.fec ~ dunif(0, 20)

#-------------------------------------------------
# 2. Derived parameters
#-------------------------------------------------
# Population growth rate
for (t in 1:(nyears-1)){
   lambda[t] <- Ntot[t+1] / Ntot[t]
   }

#-------------------------------------------------
# 3. The likelihoods of the single data sets
#-------------------------------------------------
# 3.1. Likelihood for population population count data (state-space model)
   # 3.1.1 System process
   for (t in 2:nyears){
      mean1[t] <- f[t-1] / 2 * sjuv[t-1] * Ntot[t-1]
      N1[t] ~ dpois(mean1[t])
      Nad[t] ~ dbin(sad[t-1], Ntot[t-1])
      }
   for (t in 1:nyears){
      Ntot[t] <- Nad[t] + N1[t]
      }

   # 3.1.2 Observation process
   for (t in 1:nyears){
      y[t] ~ dnorm(Ntot[t], tauy)
      }

# 3.2 Likelihood for capture-recapture data: CJS model (2 age classes)
# Multinomial likelihood
for (t in 1:2*(nyears-1)){
   m[t,1:nyears] ~ dmulti(pr[t,], r[t])
   }

# Calculate the number of released individuals
for (t in 1:2*(nyears-1)){
   r[t] <- sum(m[t,])
   }

# m-array cell probabilities for juveniles
for (t in 1:(nyears-1)){
   # Main diagonal
   q[t] <- 1-p[t]
   pr[t,t] <- sjuv[t] * p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr[t,j] <- sjuv[t]*prod(sad[(t+1):j])*prod(q[t:(j-1)])*p[j]
      } #j	
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j] <- 0
      } #j
   # Last column: probability of non-recapture
   pr[t,nyears] <- 1-sum(pr[t,1:(nyears-1)])
   } #t

# m-array cell probabilities for adults
for (t in 1:(nyears-1)){
   # Main diagonal
   pr[t+nyears-1,t] <- sad[t] * p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr[t+nyears-1,j] <- prod(sad[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t+nyears-1,j] <- 0
      } #j
   # Last column
   pr[t+nyears-1,nyears] <- 1 - sum(pr[t+nyears-1,1:(nyears-1)])
   } #t
}
",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(m = m, y = y, nyears = dim(m)[2])

# Initial values
inits <- function(){list(mean.sjuv= runif(1, 0, 1), mean.sad = runif(1, 0, 1), mean.p = runif(1, 0, 1), mean.fec = runif(1, 0, 10), N1 = rpois(dim(m)[2], 30), Nad = rpois(dim(m)[2], 30), sigma.y = runif(1, 0, 10))}  

# Parameters monitored
parameters <- c("mean.sjuv", "mean.sad", "mean.p", "mean.fec", "N1", "Nad", "Ntot", "lambda", "sigma2.y")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 1 min)
ipm.prod <- bugs(bugs.data, inits, parameters, "ipm-prod.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(ipm.prod, digits = 3)


# 11.5. IPMs for population viability analysis
# Specify model in BUGS language
sink("ipm-pred.bug")
cat("
model {
#-------------------------------------------------
#  Integrated population model
#  - Age structured model with 2 age classes: 
#		1-year old and adult (at least 2 years old)
#  - Age at first breeding = 1 year
#  - Prebreeding census, female-based
#  - All vital rates assumed to be constant
#-------------------------------------------------

#-------------------------------------------------
# 1. Define the priors for the parameters
#-------------------------------------------------
# Observation error
tauy <- pow(sigma.y, -2)
sigma.y ~ dunif(0, 50)
sigma2.y <- pow(sigma.y, 2)

# Initial population sizes
N1[1] ~ dnorm(100, 0.0001)I(0,)     # 1-year
Nad[1] ~ dnorm(100, 0.0001)I(0,)    # Adults

# Survival and recapture probabilities, as well as productivity
for (t in 1:(nyears-1+t.pred)){
   sjuv[t] <- mean.sjuv
   sad[t] <- mean.sad
   p[t] <- mean.p
   f[t] <- mean.fec
   }

mean.sjuv ~ dunif(0, 1)
mean.sad ~ dunif(0, 1)
mean.p ~ dunif(0, 1)
mean.fec ~ dunif(0, 20)

#-------------------------------------------------
# 2. Derived parameters
#-------------------------------------------------
# Population growth rate
for (t in 1:(nyears-1+t.pred)){
   lambda[t] <- Ntot[t+1] / Ntot[t]
   }

#-------------------------------------------------
# 3. The likelihoods of the single data sets
#-------------------------------------------------
# 3.1. Likelihood for population population count data (state-space model)
   # 3.1.1 System process
   for (t in 2:nyears+t.pred){
      mean1[t] <- f[t-1] / 2 * sjuv[t-1] * Ntot[t-1]
      N1[t] ~ dpois(mean1[t])
      Nad[t] ~ dbin(sad[t-1], Ntot[t-1])
      }
   for (t in 1:nyears+t.pred){
      Ntot[t] <- Nad[t] + N1[t]
      }

   # 3.1.2 Observation process
   for (t in 1:nyears){
      y[t] ~ dnorm(Ntot[t], tauy)
      }

# 3.2 Likelihood for capture-recapture data: CJS model (2 age classes)
# Multinomial likelihood
for (t in 1:2*(nyears-1)){
   m[t,1:nyears] ~ dmulti(pr[t,], r[t])
   }

# Calculate the number of released individuals
for (t in 1:2*(nyears-1)){
   r[t] <- sum(m[t,])
   }

# m-array cell probabilities for juveniles
for (t in 1:(nyears-1)){
   # Main diagonal
   q[t] <- 1-p[t]
   pr[t,t] <- sjuv[t] * p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr[t,j] <- sjuv[t]*prod(sad[(t+1):j])*prod(q[t:(j-1)])*p[j]
      } #j	
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j] <- 0
      } #j
   # Last column: probability of non-recapture
   pr[t,nyears] <- 1-sum(pr[t,1:(nyears-1)])
   } #t

# m-array cell probabilities for adults
for (t in 1:(nyears-1)){
   # Main diagonal
   pr[t+nyears-1,t] <- sad[t] * p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr[t+nyears-1,j] <- prod(sad[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t+nyears-1,j] <- 0
      } #j
   # Last column
   pr[t+nyears-1,nyears] <- 1 - sum(pr[t+nyears-1,1:(nyears-1)])
   } #t

# 3.3. Likelihood for productivity data: Poisson regression
for (t in 1:(nyears-1)){
   J[t] ~ dpois(rho[t])
   rho[t] <- R[t]*f[t]
   }
}
",fill = TRUE)
sink()

# Give the number of future years for which population size shall be estimated
t.pred <- 5

# Bundle data
bugs.data <- list(m = m, y = y, J = J, R = R, nyears = dim(m)[2], t.pred = t.pred)

# Initial values
inits <- function(){list(mean.sjuv = runif(1, 0, 1), mean.sad = runif(1, 0, 1), mean.p = runif(1, 0, 1), mean.fec = runif(1, 0, 10), N1 = rpois(dim(m)[2]+ t.pred, 30), Nad = rpois(dim(m)[2]+ t.pred, 30), sigma.y = runif(1, 0, 10))}  

# Parameters monitored
parameters <- c("mean.sjuv", "mean.sad", "mean.p", "mean.fec", "N1", "Nad", "Ntot", "lambda", "sigma2.y")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 1 min)
ipm.pred <- bugs(bugs.data, inits, parameters, "ipm-pred.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Produce Fig. 11-5
par(cex = 1.2)
lower <- upper <- numeric()
for (i in 1:15){
   lower[i] <- quantile(ipm.pred$sims.list$Ntot[,i], 0.025)
   upper[i] <- quantile(ipm.pred$sims.list$Ntot[,i], 0.975)
   }
plot(ipm.pred$mean$Ntot, type = "b", ylim = c(10, max(upper)), ylab = "Population size", xlab = "Year", las = 1, pch = 16, col = "blue", frame = F)
segments(1:15, lower, 1:15, upper, col = "blue")
points(y, type = "b", col = "black", lty = 2, pch = 16)
legend(x = 1, y = 80, legend = c("Counts", "Estimates"), pch = c(16, 16), col = c("black", "blue"), lty = c(2, 1), bty = "n")
 
mean(ipm.pred$sims.list$Ntot[,15]<30)






# 11.6. Real data example: hoopoe population dynamics 
# Load data
nyears <- 9	  # Number of years
t.pred<-5
# Capture-recapture data: m-array of juveniles and adults (these are males and females together)
marray.j <- matrix (c(15, 3, 0, 0, 0, 0, 0, 0, 198, 0, 34, 9, 1, 0, 0, 0, 0, 287, 0, 0, 56, 8, 1, 0, 0, 0, 455, 0, 0, 0, 48, 3, 1, 0, 0, 518, 0, 0, 0, 0, 45, 13, 2, 0, 463, 0, 0, 0, 0, 0, 27, 7, 0, 493, 0, 0, 0, 0, 0, 0, 37, 3, 434, 0, 0, 0, 0, 0, 0, 0, 39, 405), nrow = 8, ncol = 9, byrow = TRUE)
marray.a <- matrix(c(14, 2, 0, 0, 0, 0, 0, 0, 43, 0, 22, 4, 0, 0, 0, 0, 0, 44, 0, 0, 34, 2, 0, 0, 0, 0, 79, 0, 0, 0, 51, 3, 0, 0, 0, 94, 0, 0, 0, 0, 45, 3, 0, 0, 118, 0, 0, 0, 0, 0, 44, 3, 0, 113, 0, 0, 0, 0, 0, 0, 48, 2, 99, 0, 0, 0, 0, 0, 0, 0, 51, 90), nrow = 8, ncol = 9, byrow = TRUE)

# Population count data
popcount <- c(32, 42, 64, 85, 82, 78, 73, 69, 79,NA,NA,NA,NA,NA)

# Productivity data
J <- c(189, 274, 398, 538, 520, 476, 463, 438) # Number of offspring 
R <- c(28, 36, 57, 77, 81, 83, 77, 72)         # Number of surveyed broods 

# Specify model in BUGS language
sink("ipm.hoopoe.bug")
cat("
model {
#------------------------------------------------------------
#  Integrated population model
#  - Age structured model with 2 age classes: 
#		1-year old and adult (at least 2 years old)
#  - Age at first breeding = 1 year
#  - Prebreeding census, female-based
#  - All vital rates are assumed to be time-dependent (random)
#  - Explicit estimation of immigration
#-------------------------------------------------------------

#----------------------------------------
# 1. Define the priors for the parameters
#----------------------------------------
# Initial population sizes
N1[1] ~ dnorm(100, 0.0001)I(0,)           # 1-year old individuals
NadSurv[1] ~ dnorm(100, 0.0001)I(0,)      # Adults >= 2 years
Nadimm[1] ~ dnorm(100, 0.0001)I(0,)       # Immigrants

# Mean demographic parameters (on appropriate scale)
l.mphij ~ dnorm(0, 0.0001)I(-10,10)       # Bounded to help with convergence
l.mphia ~ dnorm(0, 0.0001)I(-10,10)
l.mfec ~ dnorm(0, 0.0001)I(-10,10)
l.mim ~ dnorm(0, 0.0001)I(-10,10)
l.p ~ dnorm(0, 0.0001)I(-10,10)

# Precision of standard deviations of temporal variability
sig.phij ~ dunif(0, 10)
tau.phij <- pow(sig.phij, -2)
sig.phia ~ dunif(0, 10)
tau.phia <- pow(sig.phia, -2)
sig.fec ~ dunif(0, 10)
tau.fec <- pow(sig.fec, -2)
sig.im ~ dunif(0, 10)
tau.im <- pow(sig.im, -2)

# Distribution of error terms (Bounded to help with convergence)
for (t in 1:((nyears-1)+t.pred)){
   epsilon.phij[t] ~ dnorm(0, tau.phij)I(-15,15)	
   epsilon.phia[t] ~ dnorm(0, tau.phia)I(-15,15)
   epsilon.fec[t] ~ dnorm(0, tau.fec)I(-15,15)
   epsilon.im[t] ~ dnorm(0, tau.im)I(-15,15)
   }

#-------------------------
# 2. Constrain parameters
#-------------------------
for (t in 1:((nyears-1)+t.pred)){
   logit(phij[t]) <- l.mphij + epsilon.phij[t]  # Juv. apparent survival
   logit(phia[t]) <- l.mphia + epsilon.phia[t]  # Adult apparent survival
   log(f[t]) <- l.mfec + epsilon.fec[t]         # Productivity
   log(omega[t]) <- l.mim + epsilon.im[t]       # Immigration
   logit(p[t]) <- l.p                           # Recapture probability
   }

#-----------------------
# 3. Derived parameters
#-----------------------
mphij <- exp(l.mphij)/(1+exp(l.mphij))   # Mean juvenile survival probability
mphia <- exp(l.mphia)/(1+exp(l.mphia))   # Mean adult survival probability
mfec <- exp(l.mfec)                      # Mean productivity
mim <- exp(l.mim)                        # Mean immigration rate

# Population growth rate
for (t in 1:(nyears-1)){
   lambda[t] <- Ntot[t+1] / Ntot[t]
   logla[t] <- log(lambda[t])
   }
mlam <- exp((1/(nyears-1))*sum(logla[1:(nyears-1)]))   # Geometric mean

#--------------------------------------------
# 4. The likelihoods of the single data sets
#--------------------------------------------
# 4.1. Likelihood for population population count data (state-space model)
   # 4.1.1 System process
   for (t in 2:nyears+t.pred){
      mean1[t] <- 0.5 * f[t-1] * phij[t-1] * Ntot[t-1]
      N1[t] ~ dpois(mean1[t])
      NadSurv[t] ~ dbin(phia[t-1], Ntot[t-1])
      mpo[t] <- Ntot[t-1] * omega[t-1]
      Nadimm[t] ~ dpois(mpo[t])
      }

   # 4.1.2 Observation process
   for (t in 1:nyears+t.pred){
      Ntot[t] <- NadSurv[t] + Nadimm[t] + N1[t]
      y[t] ~ dpois(Ntot[t])
      }

# 4.2 Likelihood for capture-recapture data: CJS model (2 age classes)
# Multinomial likelihood
for (t in 1:(nyears-1)){
   marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t])
   marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])
   }

# Calculate number of released individuals
for (t in 1:(nyears-1)){
   r.j[t] <- sum(marray.j[t,])
   r.a[t] <- sum(marray.a[t,])
   }

# m-array cell probabilities for juveniles
for (t in 1:(nyears-1)){
   q[t] <- 1-p[t]
   # Main diagonal
   pr.j[t,t] <- phij[t]*p[t]
   # Above main diagonal
   for (j in (t+1):(nyears-1)){
      pr.j[t,j] <- phij[t]*prod(phia[(t+1):j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      } #j
   # Last column
   pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
   } #t

# m-array cell probabilities for adults
for (t in 1:(nyears-1)){
   # Main diagonal
   pr.a[t,t] <- phia[t]*p[t]
   # above main diagonal
   for (j in (t+1):(nyears-1)){
      pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr.a[t,j] <- 0
      } #j
   # Last column
   pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
   } #t

# 4.3. Likelihood for productivity data: Poisson regression
for (t in 1:(nyears-1)){
   J[t] ~ dpois(rho[t])
   rho[t] <- R[t] * f[t]
   }
}
",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(nyears = nyears, marray.j = marray.j, marray.a = marray.a, y = popcount, J = J, R = R, t.pred=t.pred)

# Initial values
inits <- function(){list(l.mphij = rnorm(1, 0.2, 0.5), l.mphia = rnorm(1, 0.2, 0.5), l.mfec = rnorm(1, 0.2, 0.5), l.mim = rnorm(1, 0.2, 0.5), l.p = rnorm(1, 0.2, 1), 
sig.phij = runif(1, 0.1, 10), sig.phia = runif(1, 0.1, 10), sig.fec = runif(1, 0.1, 10), sig.im = runif(1, 0.1, 10), N1 = round(runif(nyears+t.pred, 1, 30), 0), NadSurv = round(runif(nyears+t.pred, 5, 20), 0), 
Nadimm = round(runif(nyears+t.pred, 1, 20), 0))}

# Parameters monitored
parameters <- c("phij", "phia", "f", "omega", "p", "lambda", "mphij", "mphia", "mfec", "mim", "mlam", "sig.phij", "sig.phia", "sig.fec", "sig.im", "N1", "NadSurv", "Nadimm", "Ntot")

# MCMC settings
ni <- 200
nt <- 4
nb <- 10
nc <- 3

# Call WinBUGS from R (BRT 5 min)
ipm.hoopoe <- bugs(bugs.data, inits, parameters, "ipm.hoopoe.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())



#########################################################################################
#
# 10. Estimation of survival, recruitment and population size using the Jolly-Seber model
# 
##########################################################################################

# 10.3. Fitting the JS model with data augmentation
# 10.3.1. The JS model as a restricted dynamic occupancy model
# Specify model in BUGS language
sink("js-rest.occ.bug")
cat("
model {
# Priors and constraints
for (i in 1:M){
   for (t in 1:(n.occasions-1)){
      phi[i,t] <- mean.phi
      } #t
   for (t in 1:n.occasions){
      p[i,t] <- mean.p
      } #t
   } #i
mean.phi ~ dunif(0, 1)
mean.p ~ dunif(0, 1)

for (t in 1:n.occasions){
   gamma[t] ~ dunif(0, 1)
   } #t

# Likelihood
for (i in 1:M){
   # First occasion
   # State process
   z[i,1] ~ dbern(gamma[1])
   # Observation process
   mu1[i] <- z[i,1] * p[i,1]
   y[i,1] ~ dbern(mu1[i])
   # Subsequent occasions
   for (t in 2:n.occasions){
      # State process
      q[i,t-1] <- 1-z[i,t-1]		# Availability for recruitment
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + gamma[t] * prod(q[i,1:(t-1)]) 
      z[i,t] ~ dbern(mu2[i,t])
      # Observation process
      mu3[i,t] <- z[i,t] * p[i,t]
      y[i,t] ~ dbern(mu3[i,t])
      } #t
   } #i 

# Calculate derived population parameters
for (t in 1:n.occasions){
   qgamma[t] <- 1-gamma[t]
   }
cprob[1] <- gamma[1]
for (t in 2:n.occasions){
   cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
   } #t
psi <- sum(cprob[])            # Inclusion probability
for (t in 1:n.occasions){
   b[t] <- cprob[t] / psi      # Entry probability
   } #t
for (i in 1:M){
   recruit[i,1] <- z[i,1]
   for (t in 2:n.occasions){
      recruit[i,t] <- (1-z[i,t-1]) * z[i,t]
      } #t
   } #i
for (t in 1:n.occasions){
   N[t] <- sum(z[1:M,t])        # Actual population size
   B[t] <- sum(recruit[1:M,t])  # Number of entries
   } #t
for (i in 1:M){
   Nind[i] <- sum(z[i,1:n.occasions])
   Nalive[i] <- 1-equals(Nind[i], 0)
   } #i
Nsuper <- sum(Nalive[])         # Superpopulation size
}
",fill=TRUE)
sink()


# 10.3.2. The JS model as a multistate model
# Specify model in BUGS language
sink("js-ms.bug")
cat("
model {

#--------------------------------------
# Parameters:
# phi: survival probability
# gamma: removal entry probability
# p: capture probability
#--------------------------------------
# States (S):
# 1 not yet entered
# 2 alive
# 3 dead
# Observations (O):
# 1 seen 
# 2 not seen
#--------------------------------------

# Priors and constraints
for (t in 1:(n.occasions-1)){
   phi[t] <- mean.phi
   gamma[t] ~ dunif(0, 1) # Prior for entry probabilities
   p[t] <- mean.p
   }

mean.phi ~ dunif(0, 1)    # Prior for mean survival
mean.p ~ dunif(0, 1)      # Prior for mean capture

# Define state-transition and observation matrices 	
for (i in 1:M){  
   # Define probabilities of state S(t+1) given S(t)
   for (t in 1:(n.occasions-1)){
      ps[1,i,t,1] <- 1-gamma[t]
      ps[1,i,t,2] <- gamma[t]
      ps[1,i,t,3] <- 0
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- phi[t]
      ps[2,i,t,3] <- 1-phi[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      
      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 0
      po[1,i,t,2] <- 1
      po[2,i,t,1] <- p[t]
      po[2,i,t,2] <- 1-p[t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 1
      } #t
   } #i

# Likelihood 
for (i in 1:M){
   # Define latent state at first occasion
   z[i,1] <- 1   # Make sure that all M individuals are in state 1 at t=1
   for (t in 2:n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
   } #i

# Calculate derived population parameters
for (t in 1:(n.occasions-1)){
   qgamma[t] <- 1-gamma[t]
   }
cprob[1] <- gamma[1]
for (t in 2:(n.occasions-1)){
   cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
   } #t
psi <- sum(cprob[])            # Inclusion probability
for (t in 1:(n.occasions-1)){
   b[t] <- cprob[t] / psi      # Entry probability
   } #t

for (i in 1:M){
   for (t in 2:n.occasions){
      al[i,t-1] <- equals(z[i,t], 2)
      } #t
   for (t in 1:(n.occasions-1)){
      d[i,t] <- equals(z[i,t]-al[i,t],0)
      } #t   
   alive[i] <- sum(al[i,])
   } #i

for (t in 1:(n.occasions-1)){
   N[t] <- sum(al[,t])        # Actual population size
   B[t] <- sum(d[,t])         # Number of entries
   } #t
for (i in 1:M){
   w[i] <- 1-equals(alive[i],0)
   } #i
Nsuper <- sum(w[])            # Superpopulation size
}
",fill = TRUE)
sink()


# 10.3.3. The superpopulation parameterization
# Specify model in BUGS language
sink("js-super.bug")
cat("
model {
# Priors and constraints
for (i in 1:M){
   for (t in 1:(n.occasions-1)){
      phi[i,t] <- mean.phi
      } #t
   for (t in 1:n.occasions){
      p[i,t] <- mean.p
      } #t
   } #i

mean.phi ~ dunif(0, 1)         # Prior for mean survival
mean.p ~ dunif(0, 1)           # Prior for mean capture
psi ~ dunif(0, 1)              # Prior for inclusion probability

# Dirichlet prior for entry probabilities
for (t in 1:n.occasions){
   beta[t] ~ dgamma(1, 1)
   b[t] <- beta[t] / sum(beta[1:n.occasions])
   }

# Convert entry probs to conditional entry probs
nu[1] <- b[1]
for (t in 2:n.occasions){
   nu[t] <- b[t] / (1-sum(b[1:(t-1)]))
   } #t

# Likelihood
for (i in 1:M){
   # First occasion
   # State process
   w[i] ~ dbern(psi)                  # Draw latent inclusion
   z[i,1] ~ dbern(nu[1])
   # Observation process
   mu1[i] <- z[i,1] * p[i,1] * w[i]
   y[i,1] ~ dbern(mu1[i])

   # Subsequent occasions
   for (t in 2:n.occasions){
      # State process
      q[i,t-1] <- 1-z[i,t-1]
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + nu[t] * prod(q[i,1:(t-1)]) 
      z[i,t] ~ dbern(mu2[i,t])
      # Observation process
      mu3[i,t] <- z[i,t] * p[i,t] * w[i]
      y[i,t] ~ dbern(mu3[i,t])
      } #t
   } #i 

# Calculate derived population parameters
for (i in 1:M){
   for (t in 1:n.occasions){
      u[i,t] <- z[i,t]*w[i]     # Deflated latent state (u)
      }
   }
for (i in 1:M){
   recruit[i,1] <- u[i,1]
   for (t in 2:n.occasions){
      recruit[i,t] <- (1-u[i,t-1]) * u[i,t]
      } #t
   } #i
for (t in 1:n.occasions){
   N[t] <- sum(u[1:M,t])        # Actual population size
   B[t] <- sum(recruit[1:M,t])  # Number of entries
   } #t
for (i in 1:M){
   Nind[i] <- sum(u[i,1:n.occasions])
   Nalive[i] <- 1-equals(Nind[i], 0)
   } #i
Nsuper <- sum(Nalive[])         # Superpopulation size
}
",fill=TRUE)
sink()



# 10.4. Models with constant survival and time-dependent entry
# Define parameter values
n.occasions <- 7                         # Number of capture occasions
N <- 400                                 # Superpopulation size
phi <- rep(0.7, n.occasions-1)           # Survival probabilities
b <- c(0.34, rep(0.11, n.occasions-1))   # Entry probabilities 
p <- rep(0.5, n.occasions)               # Capture probabilities

PHI <- matrix(rep(phi, (n.occasions-1)*N), ncol = n.occasions-1, nrow = N, byrow = T)
P <- matrix(rep(p, n.occasions*N), ncol = n.occasions, nrow = N, byrow = T)

# Function to simulate capture-recapture data under the JS model
simul.js <- function(PHI, P, b, N){
   B <- rmultinom(1, N, b) # Generate no. of entering ind. per occasion
   n.occasions <- dim(PHI)[2] + 1
   CH.sur <- CH.p <- matrix(0, ncol = n.occasions, nrow = N)
   # Define a vector with the occasion of entering the population
   ent.occ <- numeric()
   for (t in 1:n.occasions){
      ent.occ <- c(ent.occ, rep(t, B[t]))
      }
   # Simulating survival
   for (i in 1:N){
      CH.sur[i, ent.occ[i]] <- 1   # Write 1 when ind. enters the pop.
      if (ent.occ[i] == n.occasions) next
      for (t in (ent.occ[i]+1):n.occasions){
         # Bernoulli trial: has individual survived occasion?
         sur <- rbinom(1, 1, PHI[i,t-1])
         ifelse (sur==1, CH.sur[i,t] <- 1, break)
         } #t
      } #i
   # Simulating capture
   for (i in 1:N){
      CH.p[i,] <- rbinom(n.occasions, 1, P[i,])
      } #i
   # Full capture-recapture matrix
   CH <- CH.sur * CH.p
   
   # Remove individuals never captured
   cap.sum <- rowSums(CH)
   never <- which(cap.sum == 0)
   CH <- CH[-never,]
   Nt <- colSums(CH.sur)    # Actual population size
   return(list(CH=CH, B=B, N=Nt))
   }

# Execute simulation function
sim <- simul.js(PHI, P, b, N)
CH <- sim$CH

# Augment the capture-histories by nz pseudo-individuals
nz <- 500
CH.aug <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz))

# Bundle data
bugs.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1])

# Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = CH.aug)}  

# Parameters monitored
parameters <- c("psi", "mean.p", "mean.phi", "b", "Nsuper", "N", "B", "gamma")

# MCMC settings
ni <- 1000
nt <- 3
nb <- 500
nc <- 3

library(R2WinBUGS)


# Call WinBUGS from R (BRT 11 min)
js.occ <- bugs(bugs.data, inits, parameters, "js-rest.occ.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, working.directory = getwd())

print(js.occ, digits = 3)


# 10.4.2 Analysis of the JS model as a multistate model
# Add dummy occasion
CH.du <- cbind(rep(0, dim(CH)[1]), CH)

# Augment data
nz <- 500
CH.ms <- rbind(CH.du, matrix(0, ncol = dim(CH.du)[2], nrow = nz))

# Recode CH matrix: a 0 is not allowed in WinBUGS!
CH.ms[CH.ms==0] <- 2                     # Not seen = 2, seen = 1

# Bundle data
bugs.data <- list(y = CH.ms, n.occasions = dim(CH.ms)[2], M = dim(CH.ms)[1])

# Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = cbind(rep(NA, dim(CH.ms)[1]), CH.ms[,-1]))}    

# Parameters monitored
parameters <- c("mean.p", "mean.phi", "b", "Nsuper", "N", "B")

# MCMC settings
ni <- 20000
nt <- 3
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 32 min)
js.ms <- bugs(bugs.data, inits, parameters, "js-ms.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

print(js.ms, digits = 3)


# 10.4.3 Analysis of the JS model under the superpopulation parameterization
# Augment capture-histories by nz pseudo-individuals
nz <- 500
CH.aug <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz))

# Bundle data
bugs.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1])

# Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), psi = runif(1, 0, 1), z = CH.aug)}  

# Parameters monitored
parameters <- c("psi", "mean.p", "mean.phi", "b", "Nsuper", "N", "B", "nu")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call WinBUGS from R (BRT 40 min)
js.super <- bugs(bugs.data, inits, parameters, "js-super.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

print(js.super, digits = 3)


# 10.4.4 Comparison of estimates
# Code to produce Fig. 10-4
par(mfrow = c(1,2), mar = c(5, 6, 2, 1), mgp=c(3.4, 1, 0), las = 1)
plot(density(js.occ$sims.list$Nsuper), main = "", xlab = "", ylab = "Density", frame = FALSE, lwd = 2, ylim=c(0, 0.023), col = "blue")
points(density(js.ms$sims.list$Nsuper), type = "l", lty = 2, col = "blue", lwd = 2)
points(density(js.super$sims.list$Nsuper), type = "l", lty = 3, col = "blue", lwd = 2)
abline(v = N, col = "red", lwd = 2)
mtext("Size of superpopulation", 1, line = 3)
text(x = 470, y = 0.02, "(a)")

b1.lower <- b2.lower <- b3.lower <- b1.upper <- b2.upper <- b3.upper <- numeric()
for (t in 1:n.occasions){
   b1.lower[t] <- quantile(js.occ$sims.list$b[,t], 0.025)
   b2.lower[t] <- quantile(js.ms$sims.list$b[,t], 0.025)
   b3.lower[t] <- quantile(js.super$sims.list$b[,t], 0.025)
   b1.upper[t] <- quantile(js.occ$sims.list$b[,t], 0.975)
   b2.upper[t] <- quantile(js.ms$sims.list$b[,t], 0.975)
   b3.upper[t] <- quantile(js.super$sims.list$b[,t], 0.975)
   }

time
js.occ$mean$b

time <- 1:n.occasions
plot(x = time-0.25, y = js.occ$mean$b, xlab = "", ylab = "Entry probability", frame = FALSE, las = 1, xlim = c(0.5, 7.5), pch = 16,  ylim = c(0, max(c(b1.upper, b2.upper))))
segments(time-0.25, b1.lower, time-0.25, b1.upper)
points(x = time, y = js.ms$mean$b, pch = 1)
segments(time, b2.lower, time, b2.upper)
points(x = time+0.25, y = js.super$mean$b, pch = 17)
segments(time+0.25, b3.lower, time+0.25, b3.upper)
points(x = time, y = b, pch = 18, col = "red")
mtext("Year", 1, line = 3)
text(x = 6, y = 0.39, "(b)")


# 10.5. Models with individual capture heterogeneity
# Define parameter values
n.occasions <- 8                         # Number of capture occasions
N <- 300                                 # Size of the superpopulation
phi <- rep(0.75, n.occasions-1)          # Survival probabilities
b <- c(0.37, rep(0.09, n.occasions-1))   # Entry probabilities 
mean.p <- 0.6                            # Mean capture probability
var.p <- 1                               # Indv. Variance of capture prob.
p <- plogis(rnorm(N, qlogis(mean.p), var.p^0.5))

PHI <- matrix(rep(phi, (n.occasions-1)*N), ncol = n.occasions-1, nrow = N, byrow = T)
P <- matrix(rep(p, n.occasions), ncol = n.occasions, nrow = N, byrow = F)

# Execute simulation function
sim <- simul.js(PHI, P, b, N)
CH <- sim$CH

# Specify model in BUGS language
sink("js-super-indran.bug")
cat("
model {
# Priors and constraints
for (i in 1:M){
   for (t in 1:(n.occasions-1)){
      phi[i,t] <- mean.phi
      } #t
   for (t in 1:n.occasions){
      logit(p[i,t]) <- mean.lp + epsilon[i]
      } #t
   } #i

mean.phi ~ dunif(0, 1)              # Prior for mean survival
mean.lp <- log(mean.p / (1-mean.p))
mean.p ~ dunif(0, 1)                # Prior for mean capture
for (i in 1:M){
   epsilon[i] ~ dnorm(0, tau)I(-15,15)
   }
tau <- pow(sigma, -2)
sigma ~ dunif(0, 5)                  # Prior for sd of indv. variation of p
sigma2 <- pow(sigma, 2)
psi ~ dunif(0, 1)                    # Prior for inclusion probability

# Dirichlet prior for entry probabilities
for (t in 1:n.occasions){
   beta[t] ~ dgamma(1, 1)
   b[t] <- beta[t] / sum(beta[1:n.occasions])
   }

# Convert entry probs to conditional entry probs
nu[1] <- b[1]
for (t in 2:n.occasions){
   nu[t] <- b[t] / (1-sum(b[1:(t-1)]))
   } #t

# Likelihood
for (i in 1:M){
   # First occasion
   # State process
   w[i] ~ dbern(psi)                  # Draw latent inclusion
   z[i,1] ~ dbern(nu[1])
   # Observation process
   mu1[i] <- z[i,1] * p[i,1] * w[i]
   y[i,1] ~ dbern(mu1[i])

   # Subsequent occasions
   for (t in 2:n.occasions){
      # State process
      q[i,t-1] <- 1-z[i,t-1]
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + nu[t] * prod(q[i,1:(t-1)]) 
      z[i,t] ~ dbern(mu2[i,t])
      # Observation process
      mu3[i,t] <- z[i,t] * p[i,t] * w[i]
      y[i,t] ~ dbern(mu3[i,t])
      } #t
   } #i 

# Calculate derived population parameters
for (i in 1:M){
   for (t in 1:n.occasions){
      u[i,t] <- z[i,t]*w[i]     # Deflated latent state (u)
      }
   }
for (i in 1:M){
   recruit[i,1] <- u[i,1]
   for (t in 2:n.occasions){
      recruit[i,t] <- (1-u[i,t-1]) * u[i,t]
      } #t
   } #i
for (t in 1:n.occasions){
   N[t] <- sum(u[1:M,t])        # Actual population size
   B[t] <- sum(recruit[1:M,t])  # Number of entries
   } #t
for (i in 1:M){
   Nind[i] <- sum(u[i,1:n.occasions])
   Nalive[i] <- 1-equals(Nind[i], 0)
   } #i
Nsuper <- sum(Nalive[])         # Superpopulation size
}
",fill=TRUE)
sink()

# Augment the capture-histories by nz pseudo-individuals
nz <- 300
CH.aug <- rbind(CH, matrix(0, ncol = dim(CH)[2], nrow = nz))

# Bundle data
bugs.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1])

# Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), sigma = runif(1, 0, 1), z = CH.aug)}  

# Parameters monitored
parameters <- c("sigma2","psi", "mean.p", "mean.phi", "N", "Nsuper", "b", "B")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 179 min)
js.ran <- bugs(bugs.data, inits, parameters, "js-super-indran.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

print(js.ran, digits = 3)


# 10.7. Analysis of a real data set: survival, recruitment and population size of Leisler’s bats
# Specify model in BUGS language
sink("js-tempran.bug")
cat("
model {
# Priors and constraints
for (i in 1:M){
   for (t in 1:(n.occasions-1)){
      logit(phi[i,t]) <- mean.lphi + epsilon[t]
      } #t
   for (t in 1:n.occasions){
      p[i,t] <- mean.p
      } #t
   } #i
mean.p ~ dunif(0, 1)                # Prior for mean capture
mean.phi ~ dunif(0, 1)              # Prior for mean survival
mean.lphi <- log(mean.phi / (1-mean.phi))
for (t in 1:(n.occasions-1)){
   epsilon[t] ~ dnorm(0, tau)
   }
tau <- pow(sigma, -2)
sigma ~ dunif(0, 5)              # Prior for sd of indv. variation of phi
sigma2 <- pow(sigma, 2)

for (t in 1:n.occasions){
   gamma[t] ~ dunif(0, 1)
   } #t

# Likelihood
for (i in 1:M){
   # First occasion
   # State process
   z[i,1] ~ dbern(gamma[1])
   mu1[i] <- z[i,1] * p[i,1]
   # Observation process
   y[i,1] ~ dbern(mu1[i])
   
   # Subsequent occasions
   for (t in 2:n.occasions){
      # State process
      q[i,t-1] <- 1-z[i,t-1]
      mu2[i,t] <- phi[i,t-1] * z[i,t-1] + gamma[t] * prod(q[i,1:(t-1)]) 
      z[i,t] ~ dbern(mu2[i,t])
      # Observation process
      mu3[i,t] <- z[i,t] * p[i,t]
      y[i,t] ~ dbern(mu3[i,t])
      } #t
   } #i

# Calculate derived population parameters
for (t in 1:n.occasions){
   qgamma[t] <- 1-gamma[t]
   }
cprob[1] <- gamma[1]
for (t in 2:n.occasions){
   cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
   } #t
psi <- sum(cprob[])            # Inclusion probability
for (t in 1:n.occasions){
   b[t] <- cprob[t] / psi      # Entry probability
   } #t
for (i in 1:M){
   recruit[i,1] <- z[i,1]
   for (t in 2:n.occasions){
      recruit[i,t] <- (1-z[i,t-1]) * z[i,t]
      } #t
   } #i
for (t in 1:n.occasions){
   N[t] <- sum(z[1:M,t])        # Actual population size
   B[t] <- sum(recruit[1:M,t])  # Number of entries
   } #t
for (i in 1:M){
   Nind[i] <- sum(z[i,1:n.occasions])
   Nalive[i] <- 1-equals(Nind[i], 0)
   } #i
Nsuper <- sum(Nalive[])         # Size of superpopulation
}
",fill=TRUE)
sink()

leis <- as.matrix(read.table("leisleri.txt", sep = " ", header = FALSE))
nz <- 300
CH.aug <- rbind(leis, matrix(0, ncol = dim(leis)[2], nrow = nz))

# Bundle data
bugs.data <- list(y = CH.aug, n.occasions = dim(CH.aug)[2], M = dim(CH.aug)[1])

# Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), sigma = runif(1, 0, 1), z = CH.aug)}  

# Parameters monitored
parameters <- c("psi", "mean.p", "sigma2", "mean.phi", "N", "Nsuper", "b", "B")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 127 min)
nl <- bugs(bugs.data, inits, parameters, "js-tempran.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

print(nl, digits = 3)

# Code to produce Fig. 10-6
# Calculate per-capita recruitment
T <- dim(leis)[2]
f <- matrix(NA, ncol = T, nrow = length(nl$sims.list$B[,1]))
for (t in 1:(T-1)){
   f[,t] <- nl$sims.list$B[,t+1] / nl$sims.list$N[,t+1]
   }
n.lower <- n.upper <- f.lower <- f.upper <- f.mean <- numeric()
for (t in 1:T){
   n.lower[t] <- quantile(nl$sims.list$N[,t], 0.025)
   n.upper[t] <- quantile(nl$sims.list$N[,t], 0.975)
   }
for (t in 1:(T-1)){
   f.lower[t] <- quantile(f[,t], 0.025)
   f.upper[t] <- quantile(f[,t], 0.975)
   f.mean[t] <- mean(f[,t])
   }

par(mfrow = c(1, 2))
plot(nl$mean$N, type = "b", pch = 19, ylab = "Population size", xlab = "", axes = F, cex = 1.5, ylim = c(10, max(n.upper)))
axis(1, at = seq(1, T, 2), labels = seq(1990, 2008, 2))
axis(1, at = 1:T, labels = rep("", T), tcl = -0.25)
axis(2, las = 1)
segments(1:T, n.lower, 1:T, n.upper)

plot(f.mean, type = "b", pch = 19, ylab = "Local per capita recruitment", xlab = "", axes = F, cex = 1.5, ylim = c(0, 0.8))
axis(1, at = seq(1, (T-1), 2), labels = seq(1991, 2008, 2))
axis(1, at = 1:(T-1), labels = rep("", T-1), tcl = -0.25)
axis(2, las = 1)
segments(1:(T-1), f.lower, 1:(T-1), f.upper)



