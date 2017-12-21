update_h=function(initialh, hiter,gam, beta, Sig,X,T){
  outh = rep(0,hiter)
  outh[1] = initialh
  outsigbeta = rep(0,hiter)
  outsigbeta[1] = get_sigmabeta_from_h(initialh, gam,Sig,X,T)
  lik = rep(0,hiter)
  for (i in 2:hiter){
    h1 = outh[i-1]
    h2 = h1
    r = runif(1, -0.1, 0.1)
    h2 = h2 + r
    if (h2 < 0){h2 = (abs(h2))}
    if (h2 > 1){h2 = 2-h2}
    ind = which(gam==1)
    sigmabeta1 = get_sigmabeta_from_h_c(h1,gam,Sig,X,T)
    sigmabeta2 = get_sigmabeta_from_h_c(h2,gam,Sig,X,T)
    lik1 = sum(dnorm(beta[ind], rep(0,length(ind)),
                     sqrt(sigmabeta1*diag(Sig)[ind]), log=TRUE))
    lik2 = sum(dnorm(beta[ind], rep(0,length(ind)),
                     sqrt(sigmabeta2*diag(Sig)[ind]), log=TRUE))
    acceptanceprob = exp(lik2-lik1)
    e = runif(1,0,1)
    if(e < acceptanceprob){
      outh[i] = h2
      outsigbeta[i] = sigmabeta2
      lik[i] = lik2
    }else{
      outh[i] = h1
      outsigbeta[i] = sigmabeta1
      lik[i] = lik1
    }
  }
  return(list(h=outh, sigbeta = outsigbeta, lik = lik))
}
