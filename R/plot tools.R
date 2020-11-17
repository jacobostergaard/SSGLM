# Functions to plotting results

ks_plot <- function(z, a = .05, col1='red', col2='dodgerblue', axlabs = TRUE){
  # Transform rescaled spike times into U(0,1) variables and plot against identity
  # Input:
  #   z:      rescaled spike times
  #   a:      confidence level for bounds
  #   col1:   color for the U(0,1) variables
  #   col2:   color for the confidence region
  #   axlabs: if TRUE, then labels both axes

  z = sort(z)
  N = length(z)
  u = 1-exp(-z)
  r = (1:N-.5)/N

  crit_val = ks_crit(a)

  x = seq(-1,2,.1)

  plot(0,0, type='n', xlim=c(0,1), ylim=c(0,1), bty='n', ann=FALSE)
  polygon(c(x,rev(x)), c(x+crit_val/sqrt(N),rev(x-crit_val/sqrt(N))), col=col2, border=NA)

  abline(0,1,col=add.alpha('black', .60), lwd=1.5)
  lines(r,u, col=add.alpha(col1, .75), lwd=1.5)

  if(axlabs){
    mtext(text = "Theoretical", side=1, line= 2.5, cex=1)
    mtext(text = "Empirical", side=2, line=3, cex=1, las=0)
  }

}


plot_par <- function(P,W=NULL,idx=NULL, a=.05, col1 = 'red', col2 = 'dodgerblue', ...){
  # Plot parameter convergence
  # Input:
  #   P:      vector of the estimated parameter values at each time step
  #   W:      variance estimates at each time step (optional), if NULL, no confidence region is given
  #   a:      confidence level for the confidence region
  #   col1:   color for the parameter estimates
  #   col2:   color for the confidence region
  #   ...:    additional input to plot() function, such as ylim.

  N = length(P)
  if(is.null(idx)){
    idx = seq(1,N,length=5000)
  }


  if(is.null(W)){
    plot(idx,P[idx], type='l', xlim=range(idx), col=col1, ...)
  } else{
    plot(idx,P[idx], type='n', xlim=range(idx), ...)
    SE = sqrt(W)
    q  = qnorm(1-a/2)
    hi = P+q*SE
    lo = P-q*SE
    polygon(c(idx,rev(idx)), c(hi[idx],rev(lo[idx])), col=col2, border=NA)
    lines(idx,P[idx], col=col1)
  }


}








plot_xpars <- function(spk, cutoff=NULL, refrac=0, xpars=c(1000,1,1000,1)){

  # Plot X-state transition parameters
  # Inputs:
  #   spk:      spike train
  #   cutoff:   cutoff interval size (to make the right burst lengths)
  #   refrac:   include a refractory period
  #   xpars:    estimated transition parameters


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

    bst.pars = xpars$bst
    rst.pars = xpars$rst

    layout(1)
    par(mar=c(2,2,2,2), oma=c(1,1,1,1), bty='n', mfrow=c(3,2))
    layout(matrix(1:4,nr=2, byrow=FALSE))

    idx = min(bst):max(bst)
    plot(idx,getF(bst.h)$y[idx], col=add.alpha(par()$fg,.5), cex=.5, pch=16)
    lines(idx,Z.cdf(idx, bst.pars[1],bst.pars[2])$cdf, col=add.alpha('red',.75), lwd=3)


    hist(bst,breaks=seq(0.75*min(bst),1.25*max(bst),length=length(bst)/5), prob=TRUE, border=NA, col=add.alpha(par()$fg,.5))
    lines(idx[-1]+.5, diff(Z.cdf(idx, bst.pars[1],bst.pars[2])$cdf), col=add.alpha('red',.75), lwd=3)

    idx = min(rst):max(rst)
    plot(idx,getF(rst.h)$y[idx], col=add.alpha(par()$fg,.5), cex=.5, pch=16)
    lines(idx,Z.cdf(idx, rst.pars[1],rst.pars[2])$cdf, col=add.alpha('red',.75), lwd=3)


    hist(rst,breaks=seq(0.75*min(rst),1.25*max(rst),length=length(rst)/5), prob=TRUE, border=NA, col=add.alpha(par()$fg,.5))
    lines(idx[-1]+.5, diff(Z.cdf(idx, rst.pars[1],rst.pars[2])$cdf), col='red', lwd=3)

    layout(1)
}

