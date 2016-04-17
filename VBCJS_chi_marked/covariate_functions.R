# Define objective function to maximize for Laplace approximation of the missing covariates
Y.objective <- function(Y.mis, i, j, Y.current, mu.b0, mu.b1, mu.eta, mu.delta, tau.alpha, tau.beta){
  result <- 0
  result <- result - (tau.alpha/tau.beta)*(Y.mis - (-mu.delta[j+1] + Y.current[i,(j+1)] + Y.current[i,(j-1)] + mu.delta[j])/2 )^2
  
  if(j < l[i]){
    result <- result + log(plogis(mu.b0+mu.b1*Y.mis))
  }else{
    phi.vec <- plogis(mu.b0+mu.b1*Y.current[i,])
    phi.vec[j] <- plogis(mu.b0+mu.b1*Y.mis)
    p.vec <- plogis(mu.eta)
    if(l[i]==t){
      chi <- 1
    }else{
      chi <- (1-phi.vec[l[i]])
      for(k in (l[i]+1):t){
        chi <- chi + prod(phi.vec[l[i]:(k-1)]*(1-p.vec[(l[i]+1):k]))*(1-phi.vec[k])
      }
    }
    result <- result + log(chi)
  }
  
  return(result)
}

# Define gradient of the objective function to maximize for Laplace approximation of the missing covariates
Y.grad <- function(Y.mis, i, j, Y.current, mu.b0, mu.b1, mu.eta, mu.delta, tau.alpha, tau.beta){
  result <- 0
  result <- result - 2*(tau.alpha/tau.beta)*Y.mis - (tau.alpha/tau.beta)*(mu.delta[j+1] - Y.current[i,(j+1)] - Y.current[i,(j-1)] - mu.delta[j])
  
  if(j < l[i]){
    result <- result + mu.b1*(1-plogis(mu.b0+mu.b1*Y.mis))
  }else{
    phi.vec <- plogis(mu.b0+mu.b1*Y.current[i,])
    phi.vec[j] <- plogis(mu.b0+mu.b1*Y.mis)
    p.vec <- plogis(mu.eta)
    if(l[i]<t){
      chi <- (1-phi.vec[l[i]])
      tmp <- -mu.b1*phi.vec[l[i]]*(1-phi.vec[l[i]])*as.numeric(j==l[i])
      for(k in (l[i]+1):t){
        chi <- chi + prod(phi.vec[l[i]:(k-1)]*(1-p.vec[(l[i]+1):k]))*(1-phi.vec[k])
        tmp <- tmp + prod(phi.vec[l[i]:(k-1)]*(1-p.vec[(l[i]+1):k]))*(1-phi.vec[k])*(as.numeric(j==k)*(-phi.vec[k])*mu.b1 + as.numeric(j<k)*(1-phi.vec[j])*mu.b1)
      }
    }
    result <- result + (1/chi)*tmp
  }
  
  return(result)
}

# Wrap the objective and gradient functions for the Laplace approximation of the missing covariates into a single function for nlm
Y.nlm.objective <- function(Y.mis, i, j, Y.current, mu.b0, mu.b1, mu.eta, mu.delta, tau.alpha, tau.beta){
  res <- -Y.objective(Y.mis, i, j, Y.current, mu.b0, mu.b1, mu.eta, mu.delta, tau.alpha, tau.beta)
  attr(res, "gradient") <- -Y.grad(Y.mis, i, j, Y.current, mu.b0, mu.b1, mu.eta, mu.delta, tau.alpha, tau.beta)
  return(res)
}

# Define function to compute variances of the missing covariates
return.Y.var <- function(Y.current, mu.b0, mu.b1, mu.eta, mu.delta, tau.alpha, tau.beta){
  Y.var <- matrix(NA, nrow(Y.current), ncol(Y.current))
  
  for( i in 1:n){
    for( j in (f[i]):(t-1)){
      if(!is.na(Y[i,j])){
        Y.var[i,j] <- 0
      }else{
        if(j < l[i]){
          Y.var[i,j] <- -(tau.alpha/tau.beta)*2 - plogis(mu.b0+mu.b1*Y.current[i,j])*(mu.b1^2)*(1-plogis(mu.b0+mu.b1*Y.current[i,j]))
        }else{
          
          phi.vec <- plogis(mu.b0+mu.b1*Y.current[i,])
          p.vec <- plogis(mu.eta)
          
          chi <- (1-phi.vec[l[i]])
          tmp <- -mu.b1*phi.vec[l[i]]*(1-phi.vec[l[i]])*as.numeric(j==l[i])
          tmp.2 <- (mu.b1^2)*(phi.vec[l[i]]*(1-phi.vec[l[i]])^2 - (phi.vec[l[i]]^2)*(1-phi.vec[l[i]]))*as.numeric(j==l[i])
          for(k in (l[i]+1):t){
            chi <- chi + prod(phi.vec[l[i]:(k-1)]*(1-p.vec[(l[i]+1):k]))*(1-phi.vec[k])
            tmp <- tmp + prod(phi.vec[l[i]:(k-1)]*(1-p.vec[(l[i]+1):k]))*(1-phi.vec[k])*(as.numeric(j==k)*(-phi.vec[k])*mu.b1 + as.numeric(j<k)*(1-phi.vec[j])*mu.b1)
            tmp.2 <- tmp.2 + prod(phi.vec[l[i]:(k-1)]*(1-p.vec[(l[i]+1):k]))*(1-phi.vec[k])*(as.numeric(j==k)*(phi.vec[k]^2)*(mu.b1^2) + as.numeric(j<k)*(1-phi.vec[j])^2*mu.b1^2) + prod(phi.vec[l[i]:(k-1)]*(1-p.vec[(l[i]+1):k]))*(1-phi.vec[k])*(as.numeric(j==k)*(-phi.vec[k])*(1-phi.vec[k])*mu.b1^2 + as.numeric(j<k)*(-phi.vec[j])*(1-phi.vec[j])*mu.b1^2)
          }
          
          
          Y.var[i,j] <- -(tau.alpha/tau.beta)*2 - (1/chi^2)*(tmp^2) + (1/chi)*tmp.2
          
        }
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