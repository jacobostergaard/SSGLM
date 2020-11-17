
log.logistic <- function(x,a,b,c=0, density=FALSE){
  # log-logistic distribution function
  # https://en.wikipedia.org/wiki/Log-logistic_distribution
  x = x/a+c
  rtn = (1+x^-b)^-1
  if(density){
    rtn = (b/a)*x^(b-1)/(1+x^b)^2
  }
  return(rtn)
}


sim_glm <- function(nsim, b0, r0, rk, bk, xbst, xrst, dt = 1){

  lag_r  = length(rk) # resting kernel
  lag_b  = length(bk) # bursting kernel
  x      = numeric(nsim)-1                          # latent process
  y      = numeric(nsim)                            # observed process
  lam    = numeric(nsim)                            # intensity process
  maxlag = max(lag_r,lag_b)

  # Simulate some data
  for(i in (maxlag+1):nsim){
    #i=lag0+1

    probx = px_x(x[i-1],c(xbst,xrst))

    x[i] = ifelse(runif(1) > probx,sign(x[i-1])*(abs(x[i-1])+1),-sign(x[i-1]))

    if(x[i] <= 0){
      kern = sum(rk*y[(i-1):(i-lag_r)])
      lam[i] = exp(r0+kern)
    }else{
      kern = sum(bk*y[(i-1):(i-lag_b)])
      lam[i] = exp(b0+kern)
    }
    y[i] = rpois(1,lam[i]*dt)#min(rpois(1,lam[i]*dt),1)
  }

  return(list(y=y,x=x,lam=lam))
}


sim_izhikevich <- function(nsim, pars){

  if(is.null(pars$V0)){
    V0 = -70
  } else{
    V0 = pars$V0
  }
  if(is.null(pars$U0)){
    U0 = pars$b*V0
  } else{
    U0 = pars$U0
  }

  V = matrix(V0, nr = nsim)
  U = matrix(U0, nr = nsim)
  I = matrix(pars$I, nr=nsim) + rnorm(nsim)*pars$sig

  tmp = izhikevich(V,U,I,c(pars$dt,pars$a,pars$b,pars$c,pars$d))

  u = tmp[(nsim+1):(2*nsim),]
  v = tmp[1:nsim,]
  s = 1*(v>=30)

  out = data.frame(spk = s, v = v, u = u, i = I)

  return(out)
}
