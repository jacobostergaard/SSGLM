mpfpp <- function(y,x, pars, usetrue = FALSE, verbose = TRUE){
  # Wrapper function for estimation


  if(class(y)!="dgTMatrix"){
    y = Matrix(y, sparse = TRUE)
  }

  # Input data
  if(is.null(x)){     # No latent states are provided
    x = -1
  }

  if(is.null(pars$Neffpct)){     # No latent states are provided
    Neffpct = .1
  } else{
    Neffpct = pars$Neffpct
  }

  # Initial parameters
  Bb = pars$Bb        # Baseline init
  Br = pars$Br
  w0 = pars$w0        # Initial variance estimate

  # Parameter settings
  Sb    = pars$Sb        # Spline parameters, bursting
  Sr    = pars$Sr        # Spline parameters, resting
  xpars = c(pars$a1, pars$b1, pars$a2, pars$b2)
  q     = pars$q         # Variance parameters for linear state filter
  M     = pars$M         # Number of particle to use in non-linear state filter
  dt    = pars$dt        # Time resolution in input data
  s_lag = pars$s_lag  # Fixed lag smoother
  if(is.null(s_lag)){ # If no fixed lag smoother parameter is provided, set to 0
    s_lag = 0
  }

  nobs  = length(y)   # Number of observations
  p = length(Bb)+length(Br)

  # Run marginalized stochastic point process filter algorithm
  PFcpp = mpfpp_cpp(y = y, x = as.numeric(x), Bb = Bb, Br = Br, Sb = Sb, Sr = Sr, xpars = xpars, s_lag_in=s_lag, Neffpct = Neffpct, M = M, w0 = w0, q = q, dt = pars$dt, usetrue = usetrue, verbose = verbose)


  # Organize output
  bstnames = paste0("b",1:length(pars$Bb)-1)
  rstnames = paste0("r",1:length(pars$Br)-1)
  parnames = c(bstnames,rstnames)

  out     = list()
  out$spk = which(y>0)
  out$Sb  = Sb
  out$Sr  = Sr
  out$N   = nobs
  out$p   = p
  out$Pbin = pars$Bb
  out$Prin = pars$Br

  out$lag = max(pars$s_lag, nrow(pars$S))

  out$P   = t(PFcpp[3+1:p,])
  out$W   = array(PFcpp[-(1:(p+3)),], dim = c(p,p,ncol(PFcpp)))

  out$L   = as.numeric(PFcpp[1,])
  out$X   = as.numeric(PFcpp[2,])
  out$prob = as.numeric(PFcpp[3,])
  out$dt  = dt

  names(out$Pbin) = bstnames
  names(out$Prin) = rstnames
  colnames(out$Sb) = bstnames[-1]
  colnames(out$Sr) = rstnames[-1]
  colnames(out$P) = parnames
  colnames(out$W) = rownames(out$W) = parnames

  if(usetrue){
    out$P = c(out$Pbin,out$Prin)
    out$W = NULL
  }

  return(out)
}


