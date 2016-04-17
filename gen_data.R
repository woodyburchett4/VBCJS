# Function that generates a capture and a covariate matrix

gen.data <- function(t, n, true.b0, true.b1, true.p, true.delta, true.tau){
  # Declare death matrix, observation matrix, and covariate matrix
  Z <- matrix(0, n, t)
  X <- matrix(0, n, t)
  Y <- matrix(NA, n, t)
  
  # Find first sampling occasion
  first.vec <- sample(1:(t-1),n,replace=TRUE)
  for( i in 1:n){
    
    # Record first capture
    X[i, first.vec[i]] <- 1
    Z[i, first.vec[i]] <- 1
    Y[i, first.vec[i]] <- rnorm(1)
    
    if(first.vec[i]!=t){
      # Generate covariates
      for(j in (first.vec[i]+1):t){
        Y[i, j] <- rnorm(1, Y[i, (j-1)] + true.delta[j],  sqrt(1/true.tau))
      }
      
      # Compute death matrix
      for(j in (first.vec[i]+1):t){
        true.phi <- plogis(true.b0 + true.b1*Y[i,(j-1)])
        Z[i, j] <- Z[i, (j-1)]*rbinom(1, 1, true.phi)
      }
      
      # Compute observation matrix
      for(j in (first.vec[i]+1):t){
        X[i, j] <- Z[i, j]*rbinom(1, 1, true.p[j])
      }
    }
  }
  
  # Remove covariates on non-observed occasions
  Y[which(X==0)] <- NA
  
  return(list(X=X,Y=Y))
}