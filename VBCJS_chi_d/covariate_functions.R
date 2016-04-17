# Define objective function to maximize for Laplace approximation of the missing covariates
Y.objective <- function(Y.mis, i, j, Y.current, D.probs, mu.b0, mu.b1, mu.delta, tau.alpha, tau.beta){
  result <- 0
  result <- result - (tau.alpha/tau.beta)*(Y.mis - (-mu.delta[j+1] + Y.current[i,(j+1)] + Y.current[i,(j-1)] + mu.delta[j])/2 )^2
  result <- result + sum(D.probs[i, (j+1):t])*log(plogis(mu.b0+mu.b1*Y.mis)) + D.probs[i, j]*log(1-plogis(mu.b0+mu.b1*Y.mis))
  
  return(result)
}

# Define gradient of the objective function to maximize for Laplace approximation of the missing covariates
Y.grad <- function(Y.mis, i, j, Y.current, D.probs, mu.b0, mu.b1, mu.delta, tau.alpha, tau.beta){
  result <- 0
  result <- result - 2*(tau.alpha/tau.beta)*Y.mis - (tau.alpha/tau.beta)*(mu.delta[j+1] - Y.current[i,(j+1)] - Y.current[i,(j-1)] - mu.delta[j])
  result <- result + sum(D.probs[i, (j+1):t])*(mu.b1/(exp(mu.b0+mu.b1*Y.mis)+1)) + D.probs[i, j]*(mu.b1/(exp(mu.b0+mu.b1*Y.mis)+1) - mu.b1)
}

# Wrap the objective and gradient functions for the Laplace approximation of the missing covariates into a single function for nlm
Y.nlm.objective <- function(Y.mis, i, j, Y.current, D.probs, mu.b0, mu.b1, mu.delta, tau.alpha, tau.beta){
  res <- -Y.objective(Y.mis, i, j, Y.current, D.probs, mu.b0, mu.b1, mu.delta, tau.alpha, tau.beta)
  attr(res, "gradient") <- -Y.grad(Y.mis, i, j, Y.current, D.probs, mu.b0, mu.b1, mu.delta, tau.alpha, tau.beta)
  return(res)
}

# Define function to compute variances of the missing covariates
return.Y.var <- function(Y.current, D.probs, mu.b0, mu.b1, mu.delta, tau.alpha, tau.beta){
  Y.var <- matrix(NA, nrow(Y.current), ncol(Y.current))
  
  for( i in 1:n){
    for( j in (f[i]):(t-1)){
      if(!is.na(Y[i,j])){
        Y.var[i,j] <- 0
      }else{
        Y.var[i,j] <- -(tau.alpha/tau.beta)*2 - sum(D.probs[i, j:t])*mean(plogis(mu.b0+mu.b1*Y.current[i,j])*(mu.b1^2)*(1-plogis(mu.b0+mu.b1*Y.current[i,j])))
        Y.var[i,j] <- -1/Y.var[i,j]
      }
    }
    if(!is.na(Y[i,t])){
      Y.var[i,t] <- 0
    }else{
      Y.var[i,t] <- 1/(tau.alpha/tau.beta)
    }
  }
  
  return( Y.var)
}
