# Run VB Model using the algorithm where p and phi are NOT assumed to be independent using the chi-likelihood, but d is used to compute the missing covariate distributions

VBCJS <- function(X, Y, f, l, priors, inits, converge.criterion=1e-5){
  
  # Load priors
  for (v in 1:length(priors)) assign(names(priors)[v], priors[[v]])
  
  # Initialize parameters
  for (v in 1:length(inits)) assign(names(inits)[v], inits[[v]])
  
  t <- ncol(X)
  n <- nrow(X)
  D.probs <- matrix(NA, n, t)
  Y.var <- matrix(as.numeric(is.na(Y)), n, t)
  
  # Load objective function and gradient for the Laplace Approximation of beta and eta
  source('VBCJS_chi_d\\beta_functions.R')
  
  # Load objective function and gradient for the Laplace Approximation of the missing covariates
  source('VBCJS_chi_d\\covariate_functions.R')
  
  # Track number of iteration
  num.iter <- 0
  result.list <- list()
  euclidean.norm.current <- 10000000
  euclidean.norm.last <- 0
  VB.time <- system.time({
    # Begin iterative algorithm
    while( abs(euclidean.norm.current - euclidean.norm.last) > converge.criterion){
      
      # Update Laplace approximation for beta and eta
      res <- optim(par=c(mu.eta[2:t], mu.b0, mu.b1), fn=beta.eta.objective, gr=beta.eta.grad, hessian=TRUE, method='BFGS', Y.current=Y.current, mu.b0.prior=mu.b0.prior, mu.b1.prior=mu.b1.prior, sd.b0.prior=sd.b0.prior, sd.b1.prior=sd.b1.prior)
      mu.b0 <- res$par[t]
      mu.b1 <- res$par[t+1]
      total.covar <- solve(res$hessian)
      cov.beta <- total.covar[t:(t+1),t:(t+1)]
      mu.eta <- c(NA,res$par[1:(t-1)])
      var.eta <- c(NA,diag(total.covar)[1:(t-1)])
      
      # Compute categorical probabilities for D
      for( i in 1:n){
        for( j in 1:t){
          if( j < l[i]){
            D.probs[i, j] <- 0
          }else if(l[i]==t & j==t){
            D.probs[i,j] <- 1
          }else if( j == l[i]){
            D.probs[i, j] <- exp( log(1-plogis(mu.b0 + mu.b1*Y.current[i,j])) + 0.5*(plogis(mu.b0 + mu.b1*Y.current[i,j])*(1-plogis(mu.b0 + mu.b1*Y.current[i,j]))*(cov.beta[1,1] + 2*cov.beta[1,2]*Y.current[i,j] + cov.beta[2,2]*Y.current[i,j]^2  + Y.var[i,j]*mu.b1^2)) )
          }else if( j == t){
            D.probs[i, j] <- exp( sum( log(1-plogis(mu.eta[(l[i]+1):t]))) + sum(log(plogis(mu.b0 + mu.b1*Y.current[i,(l[i]:(t-1))])) + 0.5*plogis(mu.b0 + mu.b1*Y.current[i,(l[i]:(t-1))])*(1-plogis(mu.b0 + mu.b1*Y.current[i,(l[i]:(t-1))]))*(cov.beta[1,1] + 2*cov.beta[1,2]*Y.current[i,(l[i]:(t-1))] + cov.beta[2,2]*Y.current[i,(l[i]:(t-1))]^2 + Y.var[i,(l[i]:(t-1))]*mu.b1^2)) )
          }else{
            D.probs[i, j] <- exp( sum( log(1-plogis(mu.eta[(l[i]+1):t]))) + sum(log(plogis(mu.b0 + mu.b1*Y.current[i,(l[i]:(j-1))])) + 0.5*plogis(mu.b0 + mu.b1*Y.current[i,(l[i]:(j-1))])*(1-plogis(mu.b0 + mu.b1*Y.current[i,(l[i]:(j-1))]))*(cov.beta[1,1] + 2*cov.beta[1,2]*Y.current[i,(l[i]:(j-1))] + cov.beta[2,2]*Y.current[i,(l[i]:(j-1))]^2 + Y.var[i,(l[i]:(j-1))]*mu.b1^2)) + log(1-plogis(mu.b0 + mu.b1*Y.current[i,j])) + 0.5*(plogis(mu.b0 + mu.b1*Y.current[i,j])*(1-plogis(mu.b0 + mu.b1*Y.current[i,j]))*(cov.beta[1,1] + 2*cov.beta[1,2]*Y.current[i,j] + cov.beta[2,2]*Y.current[i,j]^2  + Y.var[i,j]*mu.b1^2)) )
          }
        }
        D.probs[i,] <- D.probs[i,]/sum(D.probs[i,])
      }
      
      # Update Laplace approximation for missing covariates
      cov.curr <- 1
      cov.last <- 0
      cov.diff <- cov.curr - cov.last
      while(abs(cov.diff) > .1){
        cov.last <- cov.curr
        
        for( i in 1:n){
          for(j in f[i]:(t-1)){
            if(is.na(Y[i,j]) & f[i]!=t){
              Y.current[i,j] <- nlm(f=Y.nlm.objective, p=Y.current[i,j], i=i, j=j, gradtol=1e-6, Y.current=Y.current, D.probs=D.probs, mu.b0=mu.b0, mu.b1=mu.b1, mu.delta=mu.delta, tau.alpha=tau.alpha, tau.beta=tau.beta)$estimate
            }
          }
          if(is.na(Y[i,t]) & f[i]!=t){
            Y.current[i,t] <- Y.current[i,(t-1)] + mu.delta[t]
          }
        }
        
        Y.var <- return.Y.var(Y.current, D.probs=D.probs, mu.b0=mu.b0, mu.b1=mu.b1, mu.delta=mu.delta, tau.alpha=tau.alpha, tau.beta=tau.beta)
        
        cov.curr <- mean(Y.current, na.rm=TRUE)
        cov.diff <- cov.curr - cov.last
      }
      
      # Update tau
      tau.alpha <- tau.alpha.prior + sum(abs(f-t))/2
      tau.beta <- tau.beta.prior
      for(i in 1:n){
        for(j in (f[i]+1):t){
          if(f[i] != t){
            tau.beta <- tau.beta + ((Y.current[i,j] - Y.current[i,(j-1)] - mu.delta[j])^2 + sd.delta[j]^2 + Y.var[i,j] + Y.var[i,(j-1)] )/2
          }
        }
      }
      
      # Update delta
      for( j in 2:t){
        mu.delta[j] <- ((tau.alpha/tau.beta)*sum((Y.current[,j]-Y.current[,(j-1)])*as.numeric(j>f), na.rm=TRUE) + ((1/sd.delta.prior)^2)*mu.delta.prior)/((tau.alpha/tau.beta)*sum(j>f) + (1/sd.delta.prior)^2)
        sd.delta[j] <- 1/sqrt(((tau.alpha/tau.beta)*sum(j>f) + (1/sd.delta.prior)^2))
      }
      
      # Update convergence criterion and store results
      euclidean.norm.last <- euclidean.norm.current
      param.vector <- c(mu.b0, mu.b1, tau.alpha, tau.beta, mu.eta[2:t], mu.delta[2:t])
      euclidean.norm.current <- sqrt(sum(param.vector^2))
      print(euclidean.norm.current)
      num.iter <- num.iter + 1
      result.list[[num.iter]] <- list( mu.delta = mu.delta, sd.delta = sd.delta, mu.b0 = mu.b0, mu.b1 = mu.b1, cov.beta = cov.beta, tau.alpha = tau.alpha, tau.beta = tau.beta, mu.eta = mu.eta, var.eta = var.eta, num.iter=num.iter, euclidean.diff=abs(euclidean.norm.current - euclidean.norm.last))
    }
  })
  
  result.list[[num.iter]] <- list( mu.delta = mu.delta, sd.delta = sd.delta, mu.b0 = mu.b0, mu.b1 = mu.b1, cov.beta = cov.beta, tau.alpha = tau.alpha, tau.beta = tau.beta, mu.eta = mu.eta, var.eta = var.eta, time=VB.time, num.iter=num.iter, covariates=Y.current, covariates.variance=Y.var, euclidean.diff=abs(euclidean.norm.current - euclidean.norm.last))
  
  # Return results
  return( result.list)
}
