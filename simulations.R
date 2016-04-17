# This file simulates a data set, analyzes it with all the VBCJS algorithms developed so far and compares them with both each other and with a traditional MCMC approach

# Set working directory
setwd("C:/Users/Woody/Desktop/Github Repos/VBCJS/VBCJS")

# Load MCMC results for cliff swallows subset using mid-weights and the data generating functions
require(R2jags)
load('mcmc_csw_subset_midwt.gzip')
source('gen_data.R')

# Define priors
priors.list <- list("mu.b0.prior" = 0,
                    "mu.b1.prior" = 0,
                    "sd.b0.prior" = 1,
                    "sd.b1.prior" = 1,
                    "p.alpha.prior" = 1,
                    "p.beta.prior" = 1,
                    "mu.delta.prior" = 0,
                    "sd.delta.prior" = 1,
                    "tau.alpha.prior" = 1,
                    "tau.beta.prior" = 1)

# Generate Data using parameters estimated from the cliff swallows data
sample.dat <- gen.data( 5, 500, fit.standard.jags$BUGSoutput$mean$beta[1], fit.standard.jags$BUGSoutput$mean$beta[2], fit.standard.jags$BUGSoutput$mean$p, fit.standard.jags$BUGSoutput$mean$delta, fit.standard.jags$BUGSoutput$mean$tau)
X <- sample.dat$X
Y <- sample.dat$Y

t <- ncol(X)
n <- nrow(X)

# Identify first and last captures
f <- array(NA,n)
l <- array(NA,n)

for( i in 1:n){
  f[i] <- min(which(X[i,]==1))
  l[i] <- max(which(X[i,]==1))
}

# Get initial values
source('initial_values.R')
init.list <- initial.values(X,Y,f,l,priors.list)

# Run algorithm using the original VB algorithm where p and phi are assumed to be independent
source('VBCJS_original\\VBCJS_original.R')
system.time( VB.sim.results <- VBCJS(X, Y, f, l, priors.list, init.list, converge.criterion=.0001) )

# Run algorithm using the VB algorithm where p and phi are NOT assumed to be independent using the chi-likelihood
source('VBCJS_chi\\VBCJS_chi.R')
system.time( VB.sim.chi.results <- VBCJS(X, Y, f, l, priors.list, init.list, converge.criterion=.0001) )

# Run algorithm using the VB algorithm where p and phi are NOT assumed to be independent using the chi-likelihood using the marked package to optimize phi and p
source('VBCJS_chi_marked\\VBCJS_chi_marked.R')
system.time( VB.sim.chi.marked.results <- VBCJS(X, Y, f, l, priors.list, init.list, converge.criterion=.0001) )

# Run algorithm using the VB algorithm where p and phi are NOT assumed to be independent using the chi-likelihood, but d is used to compute the missing covariate distributions
source('VBCJS_chi_d\\VBCJS_chi_d.R')
system.time( VB.sim.chi.d.results <- VBCJS(X, Y, f, l, priors.list, init.list, converge.criterion=.0001) )

# Run algorithm using the original approach, then run the chi-likelihood appraoch for one iteration at the end of the algorithm
source('VBCJS_chi_last_iteration\\VBCJS_chi_last_iteration.R')
system.time( VB.sim.chi.last.results <- VBCJS(X, Y, f, l, priors.list, init.list, converge.criterion=.0001) )



# Initialize covariate and death matrix for MCMC algorithm
init.y <- init.list$Y.current
init.y[X==1] <- NA
Y[is.nan(Y)] <- NA
init.y[is.nan(init.y)] <- NA

init.z <- matrix(1, n, t)
for( i in 1:n){
  init.z[i, 1:f[i]] <- NA
  if(l[i]<t){
    init.z[i, (l[i]+1):t] <- 0
  }
}

# Load data and initial values for MCMC algorithm
data <- list(n=n, t=t, f=c(f), X=X, Y=Y) 
inits <- list( list(Y=init.y, z=init.z, beta=c(init.list$mu.b0,init.list$mu.b1), tau=init.list$tau.alpha/init.list$tau.beta, delta=init.list$mu.delta, p=init.list$p.alpha/(init.list$p.alpha+init.list$p.beta)), list(Y=init.y, z=init.z, beta=c(init.list$mu.b0,init.list$mu.b1), tau=init.list$tau.alpha/init.list$tau.beta, delta=init.list$mu.delta, p=init.list$p.alpha/(init.list$p.alpha+init.list$p.beta)), list(Y=init.y, z=init.z, beta=c(init.list$mu.b0,init.list$mu.b1), tau=init.list$tau.alpha/init.list$tau.beta, delta=init.list$mu.delta, p=init.list$p.alpha/(init.list$p.alpha+init.list$p.beta)))
# Use JAGS to fit the MCMC algorithm
require(R2jags)
system.time(fit.sim.jags <- jags(data=data, inits=inits, parameters.to.save=c("p","beta", "delta", "tau"), model.file='CJS_JAGS.R', n.iter=10000) )

# Check traceplots to ensure convergence
traceplot(fit.sim.jags, var="beta")
traceplot(fit.sim.jags, var="delta")

# Compare all methods with posterior density plots
attach.jags(fit.sim.jags)

par(mfrow=c(1,1))

VB.chi.d.final <- VB.sim.chi.d.results[[length(VB.sim.chi.d.results)]]
VB.chi.last.final <- VB.sim.chi.last.results[[length(VB.sim.chi.last.results)]]
VB.chi.final <- VB.sim.chi.results[[length(VB.sim.chi.results)]]
VB.orig <- VB.sim.results[[length(VB.sim.results)]]

curve(dnorm(x, VB.orig$mu.b0, sqrt(VB.orig$cov.beta[1,1])), -1, 1, col="green", ylab='Density', xlab=expression(beta[0]), main=expression("Density Plots of "*beta[0]*" for Simulated Data"))
curve(dnorm(x, VB.chi.final[['mu.b0']], sqrt(VB.chi.final[['cov.beta']][1,1])), add=TRUE, col="blue")
curve(dnorm(x, VB.chi.d.final[['mu.b0']], sqrt(VB.chi.d.final[['cov.beta']][1,1])), add=TRUE, col="purple")
curve(dnorm(x, VB.chi.last.final[['mu.b0']], sqrt(VB.chi.last.final[['cov.beta']][1,1])), add=TRUE, col="gold")
lines(density(beta[,1]), col="red")
legend('topleft', legend=c('MFVB', 'MFVB-Chi','MFVB-Chi-d','MFVB-Chi-Last','MCMC'), col=c('green', 'blue','red', 'purple', 'gold'), lty=1)

curve(dnorm(x, VB.orig$mu.b1, sqrt(VB.orig$cov.beta[2,2])), -0.2, 0.7, col="green", ylab='Density', xlab=expression(beta[1]), main=expression("Density Plots of "*beta[1]*" for Simulated Data"))
curve(dnorm(x, VB.chi.final[['mu.b1']], sqrt(VB.chi.final[['cov.beta']][2,2])), add=TRUE, col="blue")
curve(dnorm(x, VB.chi.d.final[['mu.b1']], sqrt(VB.chi.d.final[['cov.beta']][2,2])), add=TRUE, col="purple")
curve(dnorm(x, VB.chi.last.final[['mu.b1']], sqrt(VB.chi.last.final[['cov.beta']][2,2])), add=TRUE, col="gold")
lines(density(beta[,2]), col="red")
legend('topright', legend=c('MFVB', 'MFVB-Chi','MFVB-Chi-d','MFVB-Chi-Last','MCMC'), col=c('green', 'blue','red', 'purple', 'gold'), lty=1)

par(mfrow=c(2,2))

for(i in 2:5){
  curve(dbeta(x, VB.orig$p.alpha[i], VB.orig$p.beta[i]), col='green', ylab='Density', xlab=bquote(p[.(i)]), main=bquote("Density Plot of "*p[.(i)]), xlim=c(0,1))
  lines(density(plogis(rnorm(10000, VB.chi.final[['mu.eta']][i], sqrt(VB.chi.final[['var.eta']][i])))), col='blue')
  lines(density(plogis(rnorm(10000, VB.chi.d.final[['mu.eta']][i], sqrt(VB.chi.d.final[['var.eta']][i])))), col='purple')
  lines(density(plogis(rnorm(10000, VB.chi.last.final[['mu.eta']][i], sqrt(VB.chi.last.final[['var.eta']][i])))), col='gold')
  lines(density(p[,i]), col='red')
}

legend('topright', legend=c('MFVB', 'MFVB-Chi','MFVB-Chi-d','MFVB-Chi-Last','MCMC'), col=c('green', 'blue','red', 'purple', 'gold'), lty=1)

par(mfrow=c(2,2))

for(i in 2:5){
  curve(dnorm(x, VB.orig$mu.delta[i], VB.orig$sd.delta[i]), n=1000, col='green', ylab='Density', xlab=bquote(Delta[.(i)]), main=bquote("Density Plot of "*Delta[.(i)]), xlim=c(-1,2))
  curve(dnorm(x, VB.chi.final$mu.delta[i], VB.chi.final$sd.delta[i]), n=1000, col='blue', add=TRUE)
  curve(dnorm(x, VB.chi.d.final$mu.delta[i], VB.chi.d.final$sd.delta[i]), n=1000, col='purple', add=TRUE)
  curve(dnorm(x, VB.chi.last.final$mu.delta[i], VB.chi.last.final$sd.delta[i]), n=1000, col='gold', add=TRUE)
  lines(density(delta[,i]), col='red')
}

legend('topright', legend=c('MFVB', 'MFVB-Chi','MFVB-Chi-d','MFVB-Chi-Last','MCMC'), col=c('green', 'blue','red', 'purple', 'gold'), lty=1)

par(mfrow=c(1,1))

curve(dgamma(x, VB.orig$tau.alpha, VB.orig$tau.beta), 0.5, 2, n=1000, col="green", ylab='Density', xlab=expression(tau), main=expression("Density Plot of "*tau*" for Simulated Data"))
curve(dgamma(x, VB.chi.final[['tau.alpha']], VB.chi.final[['tau.beta']]), n=1000, add=TRUE, col="blue")
curve(dgamma(x, VB.chi.d.final[['tau.alpha']], VB.chi.d.final[['tau.beta']]), n=1000, add=TRUE, col="purple")
curve(dgamma(x, VB.chi.last.final[['tau.alpha']], VB.chi.last.final[['tau.beta']]), n=1000, add=TRUE, col="gold")
lines(density(tau), col="red")
legend('topleft', legend=c('MFVB', 'MFVB-Chi','MFVB-Chi-d','MFVB-Chi-Last','MCMC'), col=c('green', 'blue','red', 'purple', 'gold'), lty=1)

detach.jags()
