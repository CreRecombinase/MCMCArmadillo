update_Sigma = function(n,nu,X,beta,Phi,Y){
  r = Y - X %*% t(beta)
  emp = em_with_zero_mean_c(r, maxit=100)
  return(riwish(n+nu, emp*n + Phi*nu))
}
