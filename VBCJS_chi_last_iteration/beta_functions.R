# Define objective function to maximize for Laplace approximation of the betas
beta.objective <- function(par, Y.current, D.probs, mu.b0.prior, mu.b1.prior, sd.b0.prior, sd.b1.prior){
  b0 <- par[1]
  b1 <- par[2]
  
  result <- 0
  for(i in 1:n){
    for(j in f[i]:(t-1)){
      if(f[i] != t){
        result <- result + sum(D.probs[i, (j+1):t])*log(plogis(b0+b1*Y.current[i,j])) + D.probs[i, j]*log(1-plogis(b0+b1*Y.current[i,j]))
      }
    }
  }
  result <- result + log(dnorm(b0, mu.b0.prior, sd.b0.prior)) + log(dnorm(b1, mu.b1.prior, sd.b1.prior))
  
  return(result)
}

# Define gradient of the objective function to maximize for Laplace approximation of the betas
beta.grad <- function(par, Y.current, D.probs, mu.b0.prior, mu.b1.prior, sd.b0.prior, sd.b1.prior){
  b0 <- par[1]
  b1 <- par[2]
  
  result.b0 <- 0
  result.b1 <- 0
  for(i in 1:n){
    for(j in f[i]:(t-1)){
      if(f[i] != t){
        result.b0 <- result.b0 + sum(D.probs[i, (j+1):t])*(1-plogis(b0+b1*Y.current[i,j])) - D.probs[i, j]*(plogis(b0+b1*Y.current[i,j]))
        result.b1 <- result.b1 + sum(D.probs[i, (j+1):t])*Y.current[i,j]*(1-plogis(b0+b1*Y.current[i,j])) - D.probs[i, j]*Y.current[i,j]*(plogis(b0+b1*Y.current[i,j]))
      }
    }
  }
  result.b0 <- result.b0 - (1/sd.b0.prior^2)*(b0-mu.b0.prior)
  result.b1 <- result.b1 - (1/sd.b1.prior^2)*(b1-mu.b1.prior)
  
  return(result)
}

# Define function to compute the proper variance-covariance matrix for the joint distribution of the betas
return.cov.beta <- function(b0, b1, Y.current, D.probs, mu.b0.prior, mu.b1.prior, sd.b0.prior, sd.b1.prior){
  b00 <- 0
  b01 <- 0
  b11 <- 0
  
  for(i in 1:n){
    for(j in f[i]:(t-1)){
      if(f[i] != t){
        b00 <- b00 - sum(D.probs[i, j:t])*plogis(b0+b1*Y.current[i,j])*(1-plogis(b0+b1*Y.current[i,j]))
        b01 <- b01 - Y.current[i,j]*sum(D.probs[i, j:t])*plogis(b0+b1*Y.current[i,j])*(1-plogis(b0+b1*Y.current[i,j]))
        b11 <- b11 - Y.current[i,j]^2*sum(D.probs[i, j:t])*plogis(b0+b1*Y.current[i,j])*(1-plogis(b0+b1*Y.current[i,j]))
      }
    }
  }
  b00 <- b00 - 1/sd.b0.prior^2
  b11 <- b11 - 1/sd.b1.prior^2
  
  return( -solve(matrix(c(b00,b01,b01,b11), 2, 2)))
}

# Define objective function to maximize for Laplace approximation of the betas and etas
beta.eta.objective <- function(param, Y.current, mu.b0.prior, mu.b1.prior, sd.b0.prior, sd.b1.prior){
  eta <- c(NA, param[1:(t-1)])
  b0 <- param[t]
  b1 <- param[t+1]
  
  result <- 0
  
  for(i in 1:n){
    phi.vec <- plogis(b0+b1*Y.current[i,])
    phi.vec[t] <- 0
    p.vec <- plogis(eta)
    if(f[i]==l[i]){
      result <- result + 0
    }else{
      result <- result + sum(log(phi.vec[f[i]:(l[i]-1)]) + X[i,(f[i]+1):l[i]]*log(p.vec[(f[i]+1):l[i]]) + (1-X[i,(f[i]+1):l[i]])*log(1-p.vec[(f[i]+1):l[i]]))
    }
    
    if(l[i]==t){
      chi <- 1
    }else{
      chi <- (1-phi.vec[l[i]])
      for(j in (l[i]+1):t){
        chi <- chi + prod(phi.vec[l[i]:(j-1)]*(1-p.vec[(l[i]+1):j]))*(1-phi.vec[j])
      }
    }
    
    result <- result + log(chi)
    
  }
  
  result <- result + log(dnorm(b0, mu.b0.prior, sd.b0.prior))
  
  result <- result + log(dnorm(b1, mu.b1.prior, sd.b1.prior))
  
  result <- result + sum(log(dnorm(eta[2:t], 0, 1.6)))
  
  return(-result)
}

# Define gradient of the objective function to maximize for Laplace approximation of the betas and etas
beta.eta.grad <- function(param, Y.current, mu.b0.prior, mu.b1.prior, sd.b0.prior, sd.b1.prior){
  eta <- c(NA, param[1:(t-1)])
  b0 <- param[t]
  b1 <- param[t+1]
  
  result <- array(0, length(param))
  
  for(i in 1:n){
    phi.vec <- plogis(b0+b1*Y.current[i,])
    phi.vec[t] <- 0
    p.vec <- plogis(eta)
    if(f[i]==l[i]){
      result <- result + 0
    }else{
      result[1:(t-1)] <- result[1:(t-1)] + as.numeric(2:t >= f[i]+1)*as.numeric(2:t <= l[i])*(X[i,2:t]*(1-p.vec[2:t]) + (1-X[i,2:t])*(-p.vec[2:t]))
      result[t] <- result[t] + sum(1-phi.vec[f[i]:(l[i]-1)])
      result[t+1] <- result[t+1] + sum(Y.current[i,f[i]:(l[i]-1)]*(1-phi.vec[f[i]:(l[i]-1)]))
    }
    
    if(l[i]<t){
      chi.2 <- array(0, length(param))
      chi <- (1-phi.vec[l[i]])
      chi.2[t] <- -phi.vec[l[i]]*(1-phi.vec[l[i]])
      chi.2[t+1] <- -Y.current[i,l[i]]*phi.vec[l[i]]*(1-phi.vec[l[i]])
      for(j in (l[i]+1):t){
        chi <- chi + prod(phi.vec[l[i]:(j-1)]*(1-p.vec[(l[i]+1):j]))*(1-phi.vec[j])
        chi.2[1:(t-1)] <- chi.2[1:(t-1)] + as.numeric(2:t <= j)*as.numeric(2:t >= l[i]+1)*prod(phi.vec[l[i]:(j-1)]*(1-p.vec[(l[i]+1):j]))*(1-phi.vec[j])*(-p.vec[2:t])
        chi.2[t] <- chi.2[t] + prod(phi.vec[l[i]:(j-1)]*(1-p.vec[(l[i]+1):j]))*(1-phi.vec[j])*(sum(1-phi.vec[l[i]:(j-1)]) - phi.vec[j])
        chi.2[t+1] <- chi.2[t+1] + prod(phi.vec[l[i]:(j-1)]*(1-p.vec[(l[i]+1):j]))*(1-phi.vec[j])*(sum(Y.current[i,l[i]:(j-1)]*(1-phi.vec[l[i]:(j-1)])) - Y.current[i,j]*phi.vec[j])
      }
      result[1:(t-1)] <- result[1:(t-1)] + as.numeric((2:t)>l[i])*(1/chi)*(chi.2[1:(t-1)])
      result[t] <- result[t] + (1/chi)*(chi.2[t])
      result[t+1] <- result[t+1] + (1/chi)*(chi.2[t+1])
    }
    
  }
  
  result[1:(t-1)] <- result[1:(t-1)] - param[1:(t-1)]/(1.6^2)
  result[t] <- result[t] - b0
  result[t+1] <- result[t+1] - b1
  
  return(-result)
  
}