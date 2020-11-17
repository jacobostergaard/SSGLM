# Functions for analyzing data


get_intensity <- function(y, x, B0.bst, k.bst, B0.rst, k.rst){
  # Calculate the estimated intensity using:
  #   y:      observed spike train
  #   x:      decoded latent state
  #   B0.bst: base line parameter for bursting
  #   k.bst:  bursting kernel (non-exponential domain)
  #   B0.rst: base line parameter for resting
  #   k.rst:  resting kernel (non-exponential domain)

  if(length(y)!=length(x)){
    stop("Spike train y and latent state x must be the same length")
  }

  if(is.null(k.bst)){ # if no bursting kernel, use baseline parameter only
    tmp.b = B0.bst
  } else{
    tmp.b = as.numeric( filter(y,  k.bst, sides=1) ) + B0.bst
  }

  if(is.null(k.rst)){ # if no resting kernel, use baseline parameter only
    tmp.r = B0.rst
  } else{
    tmp.r = as.numeric( filter(y,  k.rst, sides=1) ) + B0.rst
  }

  tmp.mode = 1*(x>0) # spike train mode according to the latent state

  out = tmp.b*tmp.mode+tmp.r*(1-tmp.mode) # estimated log-intensity
  out = exp(out)                          # output is intensity

  return(out)

}


rescale_waiting_times <- function(spk, lambda, dt=1){
  # Use time rescaling to find rescaled waiting times. These should correspond to a unit exponential
  # Inputs:
  #   spk:    spike times
  #   lambda: intensity
  #   dt:     time resolution (optional)

  Nt   = length(spk)
  z    = numeric(Nt)

  z[1] = sum(lambda[1:(spk[1]-1)] + lambda[2:spk[1]], na.rm = TRUE)/2
  for(i in 2:Nt){
    # Simple integration
      # z[i] = sum(  lambda[  (spk[i-1]+1) : spk[i]  ])

    # Use over- and under-sums and take average
      # hi = sum(  lambda[  (spk[i-1]+1) : spk[i]  ])
      # lo = sum(  lambda[  (spk[i-1]+1) : (spk[i]-1)  ])
      # z[i] = (hi+lo)/2

    # Trapezoidal rule
      idx1 = (spk[i-1]+1)
      idx2 = spk[i]
      z[i] = sum(lambda[idx1:(idx2-1)]  + lambda[(idx1+1):idx2])/2
  }

  z = z*dt # rescaled spike times and scaled in accordance with time resolution

  return(z)
}



get_cov <- function(W, p, frac=NULL){
  # Extract covariances from estimated covariance matrices
  # Input:
  #   W:  estimated covariance matrices at each timestep (array which is p x p x N)
  #   p:  vector of length 2 with the indices for the desired parameter combination

  covs = W[p[1],p[2],]
  out  = covs

  if(!is.null(frac)){
    N    = length(out)
    idx  = round(N*(1-frac)):N
    if(frac ==1){
      idx = 1:N
    }
    out = mean(out[idx])
  }

  return(out)
}

get_cor <- function(W,p, frac=NULL){
  # Extract correlations from estiamted covariance matrices
  # Input:
  #   W:  estimated covariance matrices at each timestep (array which is p x p x N)
  #   p:  vector of length 2 with the indices or names for the desired parameter combination

  if(class(p)=="character"){
    p = which(colnames(W[,,1]) %in% p)
  }

  covs = W[p[1],p[2],]
  std1 = sqrt(W[p[1],p[1],])
  std2 = sqrt(W[p[2],p[2],])

  out  = covs/(std1*std2)

  if(!is.null(frac)){
    N    = length(out)
    idx  = round(N*(1-frac)):N
    if(frac ==1){
      idx = 1:N
    }
    out = mean(out[idx])
  }

  return(out)
}


P_est <- function(P, frac){
  # Calculate parameter estimates for the last 'frac' ratio of the data
  # Input:
  #   P:    parameter estimates for all time steps
  #   frac: fraction of estimates to use, if =0 then only last estimate is used, if =1 all estimates are used

  N    = nrow(P)
  idx  = round(N*(1-frac)):N
  if(frac ==1){
    idx = 1:N
  }

  out = apply(P[idx,],2,mean)

  out = data.frame(t(out))
  return(out)

}




getF <- function(tmphist){
  # Return density estimate from a histogram in (x,y) coordinates format
  return(data.frame(x=tmphist$mids,y=cumsum(tmphist$density/sum(tmphist$density))))
}

Z.cdf <- function(z, a, b){
  # Return the cumulative distribution function values at z for the latent state transitions

  tmp1 = log( log_logistic(2:(max(z)+1),a,b) )              # log-prob to switch mode at X[k+1] given X[k]  = z
  tmp2 = cumsum( log(1-log_logistic(2:(max(z)+1),a,b)) )    # log-prob to reach X[k] = z
  z.pdf = exp( tmp1 + tmp2 )
  z.cdf = cumsum(z.pdf)

  out = data.frame(z=z, cdf=z.cdf[z], pdf = z.pdf[z])

  return(out)
}

optim_xpars <- function(init,tmp.hist){
  # Find x parameters by optimization to emipirical distributions

  tmp.x = getF(tmp.hist)$x
  tmp.y = getF(tmp.hist)$y

  minf <- function(pars){
    tmp1 = Z.cdf(tmp.x,pars[1],pars[2])$cdf
    tmp2 = tmp.y
    return(sum((tmp1-tmp2)^2))
  }

  tmp.par = optim(init,minf,method = "Nelder-Mead")$par
  tmp.par
}




burst_info <- function(spk, cutoff=NULL, refrac = 0){
  # Extract information on observed bursts: length in timesteps and number of spikes in each burst
  # Inputs:
  #   spk:      observed spike times
  #   cutoff:   cutoff interval size
  #   refrac:   include a refractory period

  isi = diff(spk)

  if(is.null(cutoff)){
    tmp     = sort(isi) # sort isi's according to size
    tmp2    = sort(diff(tmp))
    # cutoff = tmp[which(diff(tmp) == max(diff(tmp)) )+1]-1 # distribution is bimodal, find where second mode starts, subtract 1
    cutoff = tmp2[which(diff(tmp2)==max(diff(tmp2)))+1]-1
  }

  idx = which(isi>cutoff)
  if(length(idx)>0){
    bst = list()
    bst[[1]]   = spk[1:idx[1]]
    if(length(idx)>1){
      for(i in 2:length(idx)){
        bst[[i]]   = spk[(idx[i-1]+1):idx[i]]
      }
    }
    bst_length = unlist(lapply(bst, function(x) diff(range(x))))+1 + refrac # take the interval length, add 1 to include the end point spike, add refractory period
    no_spikes  = unlist(lapply(bst, length)) # number of spikes in each burst
  } else{
    bst_length = NULL
    no_spikes  = NULL
  }

  out = list(burstlength = bst_length, no_spikes = no_spikes)

  return(out)
}




find_xpars       <- function(spk, cutoff=NULL, refrac=0, init=c(1000,1,1000,1)){
  # Estimate X-state transition parameters
  # Inputs:
  #   spk:      spike train
  #   cutoff:   cutoff interval size
  #   refrac:   include a refractory period
  #   init:     initial parameter estimates


  if(is.null(cutoff)){
    tmp     = sort(isi) # sort isi's according to size
    tmp2    = sort(diff(tmp))
    # cutoff = tmp[which(diff(tmp) == max(diff(tmp)) )+1]-1 # distribution is bimodal, find where second mode starts, subtract 1
    cutoff = tmp2[which(diff(tmp2)==max(diff(tmp2)))+1]-1
  }

  bst   = burst_info(spk,cutoff,refrac)$burstlength
  isi   = diff(spk)
  rst   = isi[which(isi>cutoff)]

  bst.h = hist(bst, breaks=seq(.5,max(bst)+.5,1), plot=FALSE)
  rst.h = hist(rst, breaks=seq(.5,max(rst)+.5,1), plot=FALSE)

  bst.pars = optim_xpars(init[1:2],bst.h)
  rst.pars = optim_xpars(init[3:4],rst.h)
  out = list()
  out$pars = data.frame(bst = bst.pars, rst = rst.pars, row.names = c("alpha","beta"))
  out$bst = getF(bst.h)
  out$rst = getF(rst.h)
  out$pars

  return(out)
}






make_mpfpp_par_input <- function(Bb0, Bb=NULL, Sb, Br0, Br=NULL, Sr, xpars, M=10, q0=1e-4, q1 =1e-3, s_lag = 1e4, dt=1){
  # Set the input par list for the mpfpp function

  pars    = list()
  pars$w0 = 1
  pars$Neffpct = .1
  pars$Bb = c(Bb0, Bb) #c(log(length(tmp)/diff(range(tmp))),Bb)
  pars$Br = c(Br0, Br) #log(length(tmp)/diff(range(tmp)))#Br
  pars$Sb = Sb
  pars$Sr = Sr
  # if(length(q1)==1){
  #   pars$q = c(q0,rep(q1, length(pars$Bb)-1),q0)
  # } else{
  #   pars$q = c(q0,q1,q0)
  # }
  pars$q = q0

  pars$M = M
  pars$dt = dt
  pars$s_lag = s_lag
  pars$a1 = xpars$bst[1]
  pars$b1 = xpars$bst[2]
  pars$a2 = xpars$rst[1]
  pars$b2 = xpars$rst[2]
  return(pars)
}




get_true_glm <- function(glm, N=NULL){
  # Extract results from glm simulation object

  if(is.null(N)){
    N = glm$nsim
  }

  P    = unlist(glm$simpars)[-6]
  B0b  = P[1]
  Kb   = glm$S%*%as.matrix(P[2:4])
  B0r  = P[5]
  Kr   = NULL
  L    = glm$sim$lam[1:N]
  X    = glm$sim$x[1:N]
  z    = rescale_waiting_times(glm$spk, L, glm$dt)
  ks   = suppressWarnings(ks.test(z, 'pexp'))

  out  = list(P=P, B0b=B0b, Kb = Kb, B0r = B0r, Kr=Kr, L=L, X=X, z = z, ks = ks)

  return(out)
}


get_est <- function(res, est_ratio=0.5){
  # Extract results from mpfpp simulation object

  P   = P_est(res$P, est_ratio)
  p   = length(P)
  B0b = P$b0
  B0r = P$r0

  pb  = ncol(res$Sb)
  pr  = ncol(res$Sr)

  Kb  = res$Sb%*%as.matrix(t(P[2:(2+pb-1)]))
  Kr  = res$Sr%*%as.matrix(t(P[(2+pb):(2+pb+pr-1)]))

  L   = get_intensity(as.numeric(res$y), res$X, B0b, Kb, B0r, NULL)
  X   = res$X
  z   = rescale_waiting_times(res$spk, L, res$dt)
  ks  = suppressWarnings(ks.test(z, 'pexp'))

  N    = dim(res$W)[3]
  idx  = round(N*(1-est_ratio)):N
  bCov = apply(res$W[,,idx], 1:2, mean)
  bCor = cov2cor(bCov)

  out  = list(P=P, B0b=B0b, Kb = Kb, B0r = B0r, Kr=Kr, L=L, X=X, z = z, ks = ks, cov = bCov, cor = bCor)

  return(out)
}





make_y <- function(spk,compress=100,N=NULL,verbose=TRUE){

  spk_tmp   = downsample(spk,compress)

  if(is.null(N) | N < max(spk_tmp)){
    N = max(spk_tmp)
  }
  sp_spk    = sparseVector(rep(1,length(spk_tmp)),spk_tmp, N)
  # sp_spk    = sp_spk[-(1:(spk_tmp[1]-1))] # Remove leading 0's (start at first spike)

  # if(is.null(N)){
  #   N = length(sp_spk)
  # }
  # sp_spk[1:N]
  spk_out = list(y = sp_spk[1:N], x = 1)

  if(verbose){
    cat("Data summary:",
        "\n Numer of spikes in full dataset:   ", length(spk),
        "\n Numer of spikes in reduced dataset:", length(spk_tmp),
        "\n Time compression:                  ", compress,
        "\n Number of datapoints:              ", sprintf("%.0f",N),
        "\n Total Spikes in data to analyze:   ", sum(spk_out$y))

  }
  return(spk_out)
}
