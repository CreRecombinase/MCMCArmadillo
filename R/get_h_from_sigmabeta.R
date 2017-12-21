get_h_from_sigmabeta = function(X,sigmabeta,Sigma,gamma,n,T){
  ind = which(gamma==1)
  num = sum(X^2)/n * sum(diag(Sigma)[ind]) * sigmabeta
  denom = num + sum(diag(Sigma))
  return(num/denom)
}
