doMCMC_R = function(X,Y,n,T,Phi,nu,initialbeta,
                    initialgamma, initialSigma,
                    initialsigmabeta,
                    marcor,Vbeta,
                    niter,bgiter,hiter,switer){
  outbeta = matrix(0, niter, T)
  outgam = matrix(0, niter, T)
  outSigma = array(0,dim=c(T,T,niter))
  outsb = rep(0,niter)
  outh = rep(0,niter)
  tar = matrix(0,niter,3)
  outbeta[1,] = initialbeta
  outgam[1,] = initialgamma
  outSigma[,,1] = initialSigma
  outsb[1] = initialsigmabeta
  tar[1,] = rep(0,3)
  outh[1] = get_h_from_sigmabeta(X,outsb[1], outSigma[,,1],
                                 outgam[1,],n,T)
  for (i in 2:niter){
    gam1 = outgam[i-1,]
    beta1 = outbeta[i-1,]
    Sigma1 = outSigma[,,i-1]
    sigmabeta1 = outsb[i-1]
    h1 = outh[i-1]
    bg = update_betagam_sw(X,Y,gam1,beta1,Sigma1,sigmabeta1,
                           Vbeta,bgiter,switer)
    gam2 = bg$gam[bgiter,]
    beta2 = bg$beta[bgiter,]
    Sigma2 = update_Sigma(n,nu,X,beta2,Phi,Y)
    hsig = update_h(h1,hiter,gam2,beta2,Sigma2,X,T)
    outh[i] = hsig$h[hiter]
    outsb[i] = hsig$sigbeta[hiter]
    outgam[i,] = gam2
    outbeta[i,] = beta2
    outSigma[,,i] = Sigma2
    if(is.nan(outsb[i])){outsb[i]=1000}
    tar[i,] = get_target(X,Y,outsb[i],Sigma2,gam2,beta2)
    print(i)
  }
}
