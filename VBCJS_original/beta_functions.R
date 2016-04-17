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
