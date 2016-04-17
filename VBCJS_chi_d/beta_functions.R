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