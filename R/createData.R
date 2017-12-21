createData = function(n, Phi, nu, T, sigmabeta, missingprop){
  X = sample(c(0,1,2), n, replace = TRUE)
  #pi = rbeta(1,0.5,0.5)
  pi = rbeta(1,1,1)
  gamma = rbinom(n=T, size=1, prob=pi)
  Sigma = round(riwish(nu, Phi*nu), 5)
  beta = rep(0,T)
  ind = which(gamma==1)
  for (i in 1:length(ind)){
    beta[ind[i]] = rnorm(1, 0, sqrt(sigmabeta*Sigma[ind[i], ind[i]]))
  }
  Xbeta = X %*% t(beta)
  epsilon = mvrnorm(n, rep(0,T), Sigma)
  Y = Xbeta + epsilon
  Y = as.numeric(Y)
  missing.ind = sample(1:length(Y), round(missingprop * length(Y)))
  Y[missing.ind] = NA
  Y = matrix(Y, ncol=T)
  return(list(X=X, pi=pi, gamma=gamma, beta=beta, Sigma=Sigma, Y=Y))
}
