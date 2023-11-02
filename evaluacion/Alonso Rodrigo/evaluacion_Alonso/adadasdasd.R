rm(list = ls())

# Choose constants
set.seed(1) # 'freeze' RNGs
M <- 250 # Number of sites
T <- 20 # Number of years
lambda <- 100 # Expected abundance at t = 1
gamma <- 1.02 # Population growth rate
p <- 0.6 # Detection probability
# Create array for true abundance and for counts
N <- C <- array(NA, dim = c(M, T))
# Simulate initial conditions of system: true state at t=1
N[,1] <- rpois(M, lambda)
table(N[,1]) # Summarize
# Simulate later true states
for(t in 2:T){
  N[,t] <- rpois(M, N[,t-1] * gamma)
}
# Simulate binomial observation process and generate actual counts
for(t in 1:T){
  C[,t] <- rbinom(M, N[,t], p)
}

par(mfrow = c(1, 3))
hist(N, breaks = 100, main = 'N', col = 'grey')
hist(C, breaks = 100, main = 'C', col = 'grey')
plot(N, C, xlab = 'True N', ylab = 'Observed C', frame = F)
abline(0,1)
lm(c(C) ~ c(N)) # Check slope corresponds to p .... OK !


par(mfrow = c(2, 2))
ylim <- range(N)
matplot(t(N), type = 'l', lty = 1, main = 'Trajectories of true N', frame = F,
        ylim = ylim)
matplot(t(C), type = 'l', lty = 1, main = 'Trajectories of observed C', frame = F,
        ylim = ylim)
plot(table(N[,1]), main = 'Initial N', frame = F)
plot(table(N[,T]), main = 'Final N', frame = F)

# Get Crested Tit data and look at summaries
library(AHMbook)
data(crestedTit)
str(dat <- crestedTit) # Marc prefers short names (comment from Mike)
C <- as.matrix(dat[,6:23]) # grab counts 1999:2016
year <- 1999:2016
