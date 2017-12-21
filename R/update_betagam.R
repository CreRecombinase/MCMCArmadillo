update_betagam = function(X,
                          Y,
                          gam1,
                          beta1,
                          Sigma,
                          sigmabeta,
                          Vbeta,
                          bgiter)
  {
  T = length(gam1)
  outgamma = matrix(0, bgiter, T)
  outbeta = matrix(0, bgiter, T)
  outgamma[1,] = gam1
  outbeta[1,] = beta1
  tar = rep(0,bgiter)
  for (i in 2:bgiter){
    gam1 = outgamma[i-1, ]; beta1 = outbeta[i-1,]
    temp = update_gamma(X, Y, outgamma[i-1,])
    gam2 = as.numeric(temp$newgamma);     beta2 = beta1*gam2
    ind = which(gam2==1)
    beta2[ind] = beta1[ind] + rnorm(length(ind), 0, sqrt(Vbeta))
    changeind = temp$changeind
    change = gam2[changeind]
    A = betagam_accept(X, Y, sigmabeta, Sigma, Vbeta, gam1, beta1, gam2, beta2, changeind, change)
    check = runif(1,0,1)
    if(exp(A[1])>check){
      tar[i] = A[2]
      outgamma[i,] = gam2; outbeta[i,] = beta2;
    }else{
      tar[i] = A[3]
      outgamma[i,] = gam1; outbeta[i,] = beta1;
    }
  }
  return(list(gam = outgamma, beta = outbeta, tar = tar))
}

update_gamma = function(X,
                        Y,
                        gamma)
  {
  newgamma = gamma
  T = length(gamma)
  p1 = p2 = 0.5;
  ind0 = which(gamma==0)
  ind1 = which(gamma==1)
  s = length(ind1)
  if(s==0){
    ##nothing to remove
    p1=1; p2=0
  }else if(s==T){
    #nothing to add
    p1=0; p2=1
  }
  case = sample(c(1,2), size=1, prob = c(p1,p2))
  if (case==1){
    marcor = abs(colMeans(Y*X, na.rm=TRUE)[ind0])
    if(s<(T-1)){
      add = sample(ind0, size=1, prob=marcor)
    }else{
      add = ind0
    }
    newgamma[add] = 1
    changeind = add
  }
  if (case==2){
    remove = sample(ind1, size=1)
    newgamma[remove] = 0
    changeind = remove
  }
  return(list(newgamma=newgamma, changeind=changeind))
}

betagam_accept = function(X,
                           Y,
                           sigmabeta1,
                           inputSigma,
                           Vbeta,
                           gam1,
                           beta1,
                           gam2,
                           beta2,
                           changeind,
                           change)
{
  newtarget = sum(get_target(X, Y, sigmabeta1, inputSigma, gam2, beta2))
  oldtarget = sum(get_target(X, Y, sigmabeta1, inputSigma, gam1, beta1))
  proposal_ratio = dnorm(beta1[changeind]-beta2[changeind],
                         0, sqrt(Vbeta), log = TRUE)
  s1 = sum(gam1==1)
  s2 = sum(gam2==1)
  marcor = abs(colMeans(X*Y, na.rm=TRUE))
  if(change==1){
    tempadd = marcor[changeind] / sum(marcor[gam1==0])
    proposal_ratio = -log(tempadd)-log(s2)-proposal_ratio
  }
  if(change==0){
    tempadd = marcor[changeind] / sum(marcor[gam2==0])
    proposal_ratio = log(tempadd) + log(s1) + proposal_ratio
  }
  final_ratio = newtarget - oldtarget + proposal_ratio
  return(c(final_ratio, newtarget, oldtarget, proposal_ratio))
}
