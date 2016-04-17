# This file analyzes the first 10 years of the cliff swallows data using both VBCJS algorithms and the traditional MCMC approach

# Set working directory
setwd("C:/Users/Woody/Desktop/Github Repos/VBCJS/VBCJS")

# Load data
load('csw_for_woody_1_140715.Rdata')

# Look at only those banded as adults
csw.adult.data <- csw.data[which(csw.data$BANDAGE=='A'),]
# Remove capture histories with NA in them
csw.adult.data <- csw.adult.data[-grep("NA",csw.adult.data$CH),]
# Remove wierd capture history with a "y" in it
csw.adult.data <- csw.adult.data[-131468,]
# Remove wierd covariate value that's extremely high
csw.adult.data <- csw.adult.data[-133252,]

# Create matrix of capture histories
X <- sapply(csw.adult.data$CH, function(x){as.numeric(unlist(strsplit(x,"")))})
dimnames(X) <- NULL
X <- t(X)
# Remove last capture history (because it has no covariate data with it)
X <- X[,-ncol(X)]

# Create matrix of covariates
Y.early <- as.matrix(csw.adult.data[, names(csw.adult.data)[grep("MNEARLYWT", names(csw.adult.data))]])
Y.early <- apply(Y.early, 1, as.numeric)
Y.early <- t(Y.early)
dimnames(Y.early) <- NULL

Y.mid <- as.matrix(csw.adult.data[, names(csw.adult.data)[grep("MNMIDWT", names(csw.adult.data))]])
Y.mid <- apply(Y.mid, 1, as.numeric)
Y.mid <- t(Y.mid)
dimnames(Y.mid) <- NULL

Y.late <- as.matrix(csw.adult.data[, names(csw.adult.data)[grep("MNLATEWT", names(csw.adult.data))]])
Y.late <- apply(Y.late, 1, as.numeric)
Y.late <- t(Y.late)
dimnames(Y.late) <- NULL

# Which has the least missing data?
1 - sum(is.na(Y.early[which(X==1)]))/length(Y.early[which(X==1)])
1 - sum(is.na(Y.mid[which(X==1)]))/length(Y.mid[which(X==1)])
1 - sum(is.na(Y.late[which(X==1)]))/length(Y.late[which(X==1)])
Y.all <- array(c(Y.early, Y.mid, Y.late), dim=c(dim(Y.early),3))
Y.all <- apply(Y.all, c(1,2), function(x){mean(x,na.rm=TRUE)})
1 - sum(is.na(Y.all[which(X==1)]))/length(Y.all[which(X==1)])

# Use only the mid-weights
Y <- Y.mid

# Ensure that each capture has a corresponding covariate
X[is.na(Y)] <- 0
Y[X==0] <- NA

# Remove individuals without any captures
missing.index <- which(apply(X, 1, function(x){sum(x)>0}))
X <- X[missing.index,]
Y <- Y[missing.index,]

t <- ncol(X)
n <- nrow(X)

# Identify first and last captures
f <- array(NA,n)
l <- array(NA,n)
for( i in 1:n){
  f[i] <- min(which(X[i,]==1))
  l[i] <- max(which(X[i,]==1))
}

# Take only the first 10 years
X <- X[which(f<10),1:10]
Y <- Y[which(f<10),1:10]

# Identify first and last captures again
t <- ncol(X)
n <- nrow(X)

f <- array(NA,n)
l <- array(NA,n)

for( i in 1:n){
  f[i] <- min(which(X[i,]==1))
  l[i] <- max(which(X[i,]==1))
}

# Standardize weight covariate
Y <- (Y-mean(Y, na.rm=TRUE))/sd(Y, na.rm=TRUE)

# Generate Priors
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

# Get initial values
source('initial_values.R')
init.list <- initial.values(X,Y,f,l,priors.list)

# Run algorithm using the original VB algorithm where p and phi are assumed to be independent
source('VBCJS_chi\\VBCJS_chi.R')
system.time( VB.csw.results <- VBCJS(X, Y, f, l, priors.list, init.list, converge.criterion=.0001) )

# Initialize covariate and death matrix for MCMC algorithm
init.z <- matrix(1, n, t)
for( i in 1:n){
  init.z[i, 1:f[i]] <- NA
  if(l[i]<t){
    init.z[i, (l[i]+1):t] <- 0
  }
}
init.y <- init.list$Y.current
init.y[X==1] <- NA
Y[is.nan(Y)] <- NA
init.y[is.nan(init.y)] <- NA

data <- list(n=n, t=t, f=c(f), X=X, Y=Y) 
inits <- list( list(Y=init.y, z=init.z, beta=c(init.list$mu.b0,init.list$mu.b1), tau=init.list$tau.alpha/init.list$tau.beta, delta=init.list$mu.delta, p=init.list$p.alpha/(init.list$p.alpha+init.list$p.beta)), list(Y=init.y, z=init.z, beta=c(init.list$mu.b0,init.list$mu.b1), tau=init.list$tau.alpha/init.list$tau.beta, delta=init.list$mu.delta, p=init.list$p.alpha/(init.list$p.alpha+init.list$p.beta)), list(Y=init.y, z=init.z, beta=c(init.list$mu.b0,init.list$mu.b1), tau=init.list$tau.alpha/init.list$tau.beta, delta=init.list$mu.delta, p=init.list$p.alpha/(init.list$p.alpha+init.list$p.beta)))
# Use JAGS to fit model
require(R2jags)
# system.time(fit.standard.jags <- jags(data=data, inits=inits, parameters.to.save=c("p","beta", "delta", "tau"), model.file='CJS_JAGS.R', n.iter=50000) )
load('mcmc_csw_subset_midwt.gzip')

# Check traceplots to ensure convergence
traceplot(fit.standard.jags, var="beta")
traceplot(fit.standard.jags, var="delta")

# Compare MCMC with VBCJS using posterior density plots
attach.jags(fit.standard.jags)

par(mfrow=c(1,2))
VB.final <- VB.csw.results[[length(VB.csw.results)]]
curve(dnorm(x, VB.final$mu.b0, sqrt(VB.final$cov.beta[1,1])), -1, 3, col="blue", ylab='Density', xlab=expression(beta[0]), main=expression("Density Plots of "*beta[0]))
lines(density(beta[,1]), col="red")
legend('topright', legend=c('MFVB','MCMC'), col=c('blue','red'), lty=1)

curve(dnorm(x, VB.final$mu.b1, sqrt(VB.final$cov.beta[2,2])), -0.2, 0.2, col="blue", ylab='Density', xlab=expression(beta[1]), main=expression("Density Plots of "*beta[1]))
lines(density(beta[,2]), col="red")
legend('topright', legend=c('MFVB','MCMC'), col=c('blue','red'), lty=1)

par(mfrow=c(5, 2))

for(i in 2:10){
  plot(density(plogis(rnorm(10000, VB.final$mu.eta[i], sqrt(VB.final$var.eta[i])))), xlim=c(0,1), col='blue', ylab='Density', xlab=bquote(p[.(i)]), main=bquote("Density Plots of "*p[.(i)]))
  lines(density(p[,i]), col='red')
}
plot.new()
legend('center', legend=c('MFVB','MCMC'), col=c('blue','red'), lty=1, bty='n')

par(mfrow=c(5, 2))

for(i in 2:10){
  curve(dnorm(x, VB.final$mu.delta[i], VB.final$sd.delta[i]), -2, 2, col='blue', ylab='Density', xlab=bquote(Delta[.(i)]), main=bquote("Density Plots of "*Delta[.(i)]))
  lines(density(delta[,i]), col='red')
}
plot.new()
legend('center', legend=c('MFVB','MCMC'), col=c('blue','red'), lty=1, bty='n')

par(mfrow=c(1,2))

require(mvtnorm)
beta.samp <- rmvnorm(10000, c(VB.final$mu.b0, VB.final$mu.b1), VB.final$cov.beta)
surv.curve <- function(covar, beta.samp){
  return( plogis(beta.samp[,1] + beta.samp[,2]*covar) )
}
curve(surv.curve(x, t(as.matrix(c(VB.final$mu.b0, VB.final$mu.b1)))), from=-2, to=2, xlab='Weight', ylab='Yearly Survival Probability', main='Yearly Survival Probability by Weight (MFVB)', ylim=c(0,1))
conf.band <- apply(as.array(seq(-2,2,by=.1)), 1, function(x){quantile(surv.curve(x, beta.samp), c(0.025,0.975))})
polygon(x=c(seq(-2,2,by=.1), rev(seq(-2,2,by=.1)) ), y=c(conf.band[1,],rev(conf.band[2,])), col='gray', border=NA)
curve(surv.curve(x, t(as.matrix(c(VB.final$mu.b0, VB.final$mu.b1)))), from=-2, to=2, add=TRUE)

beta.samp <- beta
surv.curve <- function(covar, beta.samp){
  return( plogis(beta.samp[,1] + beta.samp[,2]*covar) )
}
curve(surv.curve(x, t(as.matrix(colMeans(beta)))), from=-2, to=2, xlab='Weight', ylab='Yearly Survival Probability', main='Yearly Survival Probability by Weight (MCMC)', ylim=c(0,1))
conf.band <- apply(as.array(seq(-2,2,by=.1)), 1, function(x){quantile(surv.curve(x, beta.samp), c(0.025,0.975))})
polygon(x=c(seq(-2,2,by=.1), rev(seq(-2,2,by=.1)) ), y=c(conf.band[1,],rev(conf.band[2,])), col='gray', border=NA)
curve(surv.curve(x, t(as.matrix(colMeans(beta)))), from=-2, to=2, add=TRUE)


detach.jags()

