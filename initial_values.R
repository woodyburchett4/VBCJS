# This file generates initial values for all of the VBCJS algorithms

initial.values <- function(X, Y, f, l, priors){
  
  # Load priors
  for (v in 1:length(priors)) assign(names(priors)[v], priors[[v]])
  
  # Initialize parameters
  t <- ncol(X)
  n <- nrow(X)
  D.probs <- matrix(NA, n, t)
  Y.var <- matrix(as.numeric(is.na(Y)), n, t)
  
  # Initialize delta
  diff.store <- list()
  length(diff.store) <- t
  for( i in 1:n){
    for(j in 1:(t-1)){
      if(!is.na(Y[i,j]) & !is.na(Y[i,(j+1)]) ){
        diff.store[[j+1]] <- as.vector( c(diff.store[[j+1]] , Y[i,j+1] - Y[i,j]))
      }
    }
  }
  mu.delta <- unlist(lapply(diff.store, mean))
  mu.delta[1] <- 0
  mu.delta[is.na(mu.delta)] <- 0
  
  sd.delta <- unlist(lapply(diff.store, sd))/sqrt(unlist(lapply(diff.store, length)))
  sd.delta[is.na(sd.delta)] <- 1
  sd.delta[1] <- 1
  
  # Initialize tau
  tau.alpha <- tau.alpha.prior + sum(abs(f-t))/2
  tau.beta <- tau.beta.prior + (sum(abs(f-t))*mean(unlist(lapply(diff.store, sd))^2,na.rm=TRUE))/2
  
  # Initialize covariate matrix
  Y.current <- Y
  for(i in 1:n){
    for(j in (f[i]+1):t){
      if(is.na(Y.current[i,j])){
        if(mean(is.na(Y.current[i,j:t]))==1 ){
          Y.current[i,j] <- Y.current[i,j-1] + mu.delta[j]
        }else{
          next.value <- min(which(!is.na(Y.current[i,]) & (1:t)>j))
          Y.current[i,j] <- Y.current[i,j-1] + (Y.current[i,next.value]-Y.current[i,j-1])*(1/(next.value-j+1))
        }
      }
    }
  }
  
  # Use marked to get estimates for beta and p using the initial covariate values as data
  require(marked)  
  # Initialize coefficients and capture probabilities
  tmp <- data.frame(ch=apply(X, 1, function(x){paste0(x,collapse='')}), stringsAsFactors=FALSE)
  tmp <- cbind(tmp,Y.current)
  colnames(tmp)[2:(ncol(Y)+1)] <- paste0('td',1:ncol(Y))
  tmp[is.na(tmp)] <- 0
  
  tmp.proc <- process.data(tmp, model='cjs')
  design.Phi=list(time.varying=c('td'))
  design.p=list()
  design.parameters=list(Phi=design.Phi,p=design.p)
  ddl=make.design.data(tmp.proc,parameters=design.parameters)
  names(ddl$Phi)
  names(ddl$p)
  
  Phi.mod=list(formula=~td)
  p.mod=list(formula=~time)
  model=crm(tmp.proc, ddl, hessian=TRUE, model.parameters=list(Phi=Phi.mod,p=p.mod))
  
  mu.b0 <- unname(model$results$beta$Phi)[1]
  mu.b1 <- unname(model$results$beta$Phi)[2]
  cov.beta <- unname(model$results$beta.vcv[1:2,1:2])
  
  p.alpha <- array(NA, t)
  p.beta <- array(NA, t)
  mu.eta <- array(NA, t)
  var.eta <- array(NA, t)
  for(j in 2:t){
    mu.eta[j] <- unname(model$results$beta$p[1] + model$results$beta$p[j-1]*as.numeric(j>2))
    var.eta[j] <- diag(model$results$beta.vcv)[3] + (diag(model$results$beta.vcv)[j+1] + model$results$beta.vcv[3,j+1])*as.numeric(j>2)
    p.alpha[j] <- p.alpha.prior + sum(X[,j]*as.numeric(f<j)*as.numeric(j<=l))
    p.beta[j] <- p.beta.prior + sum((1-X[,j])*as.numeric(f<j)*as.numeric(j<l))
  }
  
  return(list(mu.b0=mu.b0, mu.b1=mu.b1, cov.beta=cov.beta, Y.current=Y.current, p.alpha=p.alpha, p.beta=p.beta, mu.delta=mu.delta, sd.delta=sd.delta, tau.alpha=tau.alpha, tau.beta=tau.beta, mu.eta=mu.eta, var.eta=var.eta))
  
}
  