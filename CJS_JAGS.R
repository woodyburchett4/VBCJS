model{
  
  ## Priors
  for( i in 1:(t)){
    p[i] ~ dbeta(1,1)
    delta[i] ~ dnorm(0,1)
  }
  beta[1] ~ dnorm(0,1)
  beta[2] ~ dnorm(0,1)
  
  tau ~ dgamma(1,1)
  
  ## Likelihood
  for(i in 1:n){
    z[i,f[i]] ~ dbern(1)
    X[i,f[i]] ~ dbern(1)
    
    for(j in 1:(f[i]-1)){
      Y[i,j] <- 0
    }
    
    for(j in (f[i]+1):t){
      X[i,j] ~ dbern(mu1[i,j])
      mu1[i,j] <- p[j]*z[i,j]
      z[i,j] ~ dbern(mu2[i,j])
      logit(phi[i,(j-1)]) <- beta[1] + beta[2]*Y[i,(j-1)]
      mu2[i,j] <- phi[i,(j-1)]*z[i,(j-1)]
      Y[i,j] ~ dnorm(mu3[i,j], tau)
      mu3[i,j] <- Y[i,(j-1)] + delta[j]
    }
  }
  
}