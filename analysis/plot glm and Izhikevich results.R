misc::clean_up()
library(SSGLM)
datalib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/plots/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/TeX/SSGLM plos-latex-template/gfx/"

# Load results from simulations and mpfpp runs

fn = paste0(datalib, "izh_sim.Rda")
load(fn) # Load simulated Izhikevich neuron

fn = paste0(datalib,"izh_mpfpp_res.Rda")
load(fn) # Load mpfpp results for Izhikevich neuron

fn = paste0(datalib,"ssglm_sim.Rda")
load(fn) # Load simulated SSGLM neuron

fn  = paste0(datalib,"glm_mpfpp.Rda")
load(fn) # Load SSGLM mpfpp outcome

fn  = paste0(datalib,"glm_mpfpp_res.Rda")
load(fn) # Load mpfpp results for SSGLM neuron

tru = get_true_glm(glm)





# Plot of qualitative comparison of Izhikevich and SSGLM neurons

    pdf(paste0(plotlib,"izh_glm_isi.pdf"), width=8, height=4)
      col1 = add.alpha('black',.5)
      col2 = add.alpha('red',.5)
      layout(1)
      hist(izh$isi, breaks=seq(0,10000,10), border=NA, col=col1, ylim=c(0,0.008), prob=TRUE, ann=FALSE, xlim=c(0,4000), xaxt='n')
      par(new=TRUE)
      hist(glm$isi, breaks=seq(0,10000,10), border=NA, col=col2, ylim=c(0,0.008), prob=TRUE, ann=FALSE, xlim=c(0,4000), xaxt='n')
      axis(1, at=pretty(c(0,4000)), labels = pretty(c(0,40)))
      legend("topright", c("Izhikevich neuron","SSGLM neuron"), fill=c(col1,col2), bty='n')
      mtext("Density", side=2, line=2.5)
      mtext("Interspike interval (ms)", side=1, line=2)
    dev.off()
      # quantile(izh$isi, probs=seq(0,1,.1))
      # quantile(glm$isi, probs=seq(0,1,.1))
    cat("\n# spikes from Izhikevich neuron:",length(izh$spk))
    cat("\n# spikes from SSGLM neuron:     ",length(glm$spk))


# Plot fit of transition parameters for Izhikevich neuron

  spk    = izh$spk
  cutoff = 15/izh$dt
  refrac = min(izh$isi)

  bst   = burst_info(spk,cutoff,refrac)$burstlength
  isi   = diff(spk)
  rst   = isi[which(isi>cutoff)]

  bst.h = hist(bst, breaks=seq(.5,4000+.5,1), plot=FALSE)
  rst.h = hist(rst, breaks=seq(.5,4000+.5,1), plot=FALSE)

  bst.pars = izh$xpars$bst
  rst.pars = izh$xpars$rst

  idx = 1:4000

  col1 = add.alpha('black',.25)
  col2 = add.alpha('red',.75)
  col3 = add.alpha('black',.75)

  fn = paste0(plotlib, "izh_xfit_bst.pdf")
  pdf(fn, height=3.5, width=4)
    par(mar=c(2,2,2,2), oma=c(1,1,1,1), bty='n', mfrow=c(1,1))
    plot(0,0, type='n', xlim=c(0,3000), ylim=c(0,1), xaxt='n')
    axis(1, at=pretty(c(0,3000)), pretty(c(0,30)))
    lines(idx,getF(bst.h)$y[idx], col=col3, lty=1, lwd=1)
    points(idx,getF(bst.h)$y[idx], col=col1, cex=.5, pch=16)
    lines(idx,Z.cdf(idx, bst.pars[1],bst.pars[2])$cdf, col=col2, lwd=3)
    legend("topleft", bty='n', legend=c("Observed","Fitted"), pch=c(16,NA), lty=c(NA,1), lwd=c(NA,3), col=c(col3,col2), seg.len = .85)
    mtext("CDF", side = 2, line=2)
    mtext("Burst length (ms)", side=1, line=2)
    mtext("A", side = 3, line=0, adj=0, family = 'sans', cex=1.25)
  dev.off()

  fn = paste0(plotlib, "izh_xfit_rst.pdf")
  pdf(fn, height=3.5, width=4)
    par(mar=c(2,2,2,2), oma=c(1,1,1,1), bty='n', mfrow=c(1,1))
    plot(0,0, type='n', xlim=c(1000,4000), ylim=c(0,1), xaxt='n')
    axis(1, at=pretty(c(1000,4000)), pretty(c(10,40)))
    lines(idx,getF(rst.h)$y[idx], col=col3, lty=1, lwd=1)
    points(idx,getF(rst.h)$y[idx], col=col1, cex=.5, pch=16, xlim=c(1000,4000))
    lines(idx,Z.cdf(idx, rst.pars[1],rst.pars[2])$cdf, col=col2, lwd=3)
    legend("topleft", bty='n', legend=c("Observed","Fitted"), pch=c(16,NA), lty=c(NA,1), lwd=c(NA,3), col=c(col3,col2), seg.len = .95)
    mtext("CDF", side = 2, line=2)
    mtext("Rest length (ms)", side=1, line=2)
    mtext("B", side = 3, line=0, adj=0, family = 'sans', cex=1.25)
  dev.off()





# Plot parameter convergence for SSGLM and Izhikevich mpfpp runs
  N  = izh$nsim
  fn = paste0(plotlib,"izh_par_convergence.pdf")
  pdf(fn, width=3.5, height=6)
      layout(1:4)
      par(mar=c(3,4,1,1), las=1, oma=c(0,0,1,0), bty='n')
      est.col = add.alpha('red',.75)
      ci.col = add.alpha('dodgerblue',.35)
      idx = seq(4001,N,length=5000)
      P = izh$mpfpp$P
      W = izh$mpfpp$W
      P.est = as.numeric(izh.res$P)
      for(i in 1:4){
        plot_par(P[,i], W[i,i,], idx = idx, ylim=as.numeric(P.est)[i]+c(-2,2), col1 = est.col, col2 = ci.col, xaxt='n', ann=FALSE)
        axis(1, pretty(c(0,N)), pretty(c(0,N/100)))
        mtext(bquote(beta[.(i-1)]),2,line=2.5, las=1)
        mtext("ms",1, line=2)
      }
      mtext("B", side = 3, line=-1, adj=0, family = 'sans', cex=1.05, outer=TRUE)
  dev.off()


  fn = paste0(plotlib,"glm_par_convergence.pdf")
  pdf(fn, width=3.5, height=6)
      layout(1:4)
      par(mar=c(3,4,1,1), las=1, oma=c(0,0,1,0), bty='n')
      est.col = add.alpha('red',.75)
      ci.col = add.alpha('dodgerblue',.35)
      idx = seq(4001,N,length=5000)
      P = res$P
      W = res$W
      P.est = as.numeric(glm.res$P)
      for(i in 1:4){
        plot_par(P[,i], W[i,i,], idx = idx, ylim=as.numeric(P.est)[i]+c(-2,2), col1 = est.col, col2 = ci.col, xaxt='n', ann=FALSE)
        abline(h=unlist(glm$simpars)[i], lty=3)
        axis(1, pretty(c(0,N)), pretty(c(0,N/100)))
        mtext(bquote(beta[.(i-1)]),2,line=2.5, las=1)
        mtext("ms",1, line=2)
      }
      mtext("A", side = 3, line=-1, adj=0, family = 'sans', cex=1.05, outer=TRUE)
  dev.off()




# Plot the kernel estimates

  col1 = add.alpha('black',.75)
  col2 = add.alpha('red',.75)
  layout(1)
  pdf(paste0(plotlib,"sim_kernel_nobase.pdf"), width=5, height=3)

      par(mar=c(4,4,1,1), oma=c(0,0,0,0), las=1)
      plot(exp(glm.res$Kb), type='l', bty='n', lwd=2, ylim=c(0,2.5), col=col2, xlab="",ylab="", xaxt='n')
      lines(exp(izh.res$Kb),col=col1,lwd=2, lty=1)

      # lines(exp(tru$Kb+ tru$B0b),col=col1,lwd=2, lty=1)
      # lines(exp(tru$Kb),col=col1,lwd=2, lty=3)
      axis(1, at = pretty(c(0,2000)), labels = pretty(c(0,20)))

      abline(h=1, lty=3, col=add.alpha('black',.5))
      mtext("ms",1, outer=FALSE, line=2, cex=1.25)
      # legend("topright", c("Izhikevich","SSGLM","With baseline", "Without baseline"), lty=c(NA,NA,1,3), col=c(col1,col2, col1,col1), pch=c(15,15,NA,NA), lwd=3, bty='n')
      legend("topright", c("Estimated", "True (Izhikevich)"), lty=c(1,1), col=c(col2,col1), lwd=3, bty='n')
      mtext("A", side = 3, line=0, adj=0, family = 'sans', cex=1.25)
  dev.off()

  layout(1)
  pdf(paste0(plotlib,"sim_kernel_base.pdf"), width=5, height=3)

      par(mar=c(4,4,1,1), oma=c(0,0,0,0), las=1)
      plot(exp(izh.res$Kb+izh.res$B0b), type='l', bty='n', lwd=2, ylim=c(0,7), col=col1, xlab="",ylab="", xaxt='n')
      lines(exp(glm.res$Kb+glm.res$B0b),col=col2,lwd=2, lty=1)
      axis(1, at = pretty(c(0,2000)), labels = pretty(c(0,20)))

      abline(h=1, lty=3, col=add.alpha('black',.5))
      mtext("ms",1, outer=FALSE, line=2, cex=1.25)
      # legend("topright", c("With baseline", "Without baseline"), lty=c(1,3), lwd=3, bty='n')
      legend("topright", c("Estimated", "True (Izhikevich)"), lty=c(1,1), col=c(col2,col1), lwd=3, bty='n')
      mtext("B", side = 3, line=0, adj=0, family = 'sans', cex=1.25)
  dev.off()







  col1 = add.alpha('black',.75)
  col2 = add.alpha('red',.75)
  col3 = add.alpha('dodgerblue',.75)
  pdf(paste0(plotlib,"izh_decode.pdf"), width=8, height=3.5)
      par(mar=c(0,0,0,0), oma=c(4,0,1,0))
      layout(1:3, heights = c(4,1,4))

      idx = seq(475e3,N,1)
      plot(idx,izh.res$L[idx], type='l', bty='n', axes=FALSE, col=col2, lwd=2, ann=FALSE, xlim=c(475e3,500e3))
      plot(idx,izh.res$L[idx], type='n', bty='n', axes=FALSE, col=col2, lwd=2, ylim=c(-.25,.25), ann=FALSE)

      tmpspk = izh$spk[izh$spk>min(idx)]
      points(tmpspk, rep(0,length(tmpspk)), pch=124, cex=2, col=col1)

      idx2 = idx
      idx2[which(idx2 %in% which(abs(diff(izh.res$X))>300))] = NA

      plot(idx2,izh.res$X[idx2], type='l', col=col3, bty='n', lwd=2, yaxt='n', ann=FALSE, xaxt='n')
      abline(h=0, lty=3)
      axis(1, at=pretty(c(475e3,500e3)), labels = pretty(c(4750,5000)))
      mtext("ms",1, outer=FALSE, line=3, cex=1)

      mtext("B", side = 3, line=-1, adj=0, family = 'sans', cex=1.25, outer=TRUE)
  dev.off()




  pdf(paste0(plotlib,"glm_decode.pdf"), width=8, height=3.5)
      par(mar=c(0,0,0,0), oma=c(4,0,1,0))
      layout(1:3, heights = c(4,1,4))
      idx = seq(475e3,N,1)

      plot(idx,glm.res$L[idx], type='l', bty='n', axes=FALSE, col=col2, lwd=2, ann=FALSE, xlim=c(475e3,500e3))
      plot(idx,glm.res$L[idx], type='n', bty='n', axes=FALSE, col=col2, lwd=2, ylim=c(-.25,.25), ann=FALSE)
      tmpspk = glm$spk[glm$spk>min(idx)]
      points(tmpspk, rep(0,length(tmpspk)), pch=124, cex=2, col=col1)

      idx2 = idx
      idx2[which(idx2 %in% which(abs(diff(glm.res$X))>300))] = NA
      idx3 = idx
      idx3[which(idx3 %in% which(abs(diff(tru$X))>300))] = NA

      plot(idx2, glm.res$X[idx2], type='l', col=col3, bty='n', lwd=2, yaxt='n', ann=FALSE, xaxt='n')
      lines(idx3, tru$X[idx3], type='l', col=col1, lwd=2, lty=1)
      abline(h=0, lty=3)
      axis(side = 1, at=pretty(c(475e3,500e3)), labels = pretty(c(4750,5000)))
      mtext("ms",1, outer=FALSE, line=3, cex=1.25)
      mtext("A", side = 3, line=-1, adj=0, family = 'sans', cex=1.25, outer=TRUE)
  dev.off()



  col1 = add.alpha('red',.75)
  col2 = add.alpha('dodgerblue',.25)
  pdf(paste0(plotlib,"glm_ks.pdf"), width=4, height=4)
      ks_plot(glm.res$z, col1 = col1, col2 = col2)
      mtext("A", side = 3, line=0, adj=0, family = 'sans', cex=1.25)
  dev.off()

  pdf(paste0(plotlib,"izh_ks.pdf"), width=4, height=4)
      ks_plot(izh.res$z, col1 = col1, col2 = col2)
      mtext("B", side = 3, line=0, adj=0, family = 'sans', cex=1.25)
  dev.off()


# UV burst prob plot

  pdf(paste0(plotlib,"uv_prob.pdf"), width=10, height=8)

      par(mfrow=c(2,1), bty='n', mar=c(0,0,0,0), oma=c(5,5,1,1))
      # par(mfrow=c(2,3), mar=c(0,4,1,4), oma=c(4,4,0,1))
      layout(matrix(1:2, nr=1, nc=2), widths = c(5,1))
      Cols = rev(colorRampPalette(RColorBrewer::brewer.pal(n = 11,name = "RdBu"))(1300))
      Cols = Cols[151:1150]
      Cols = add.alpha(Cols,.75)

      idx = seq(475e3,500e3)
      pct = izh$mpfpp$prob[idx]
      pct[pct==0] = 0.001
      colidx = round(pct*1e3)

      # Plot (U,V) with P(burst) coloring
      plot(0,0, type='n', ann=FALSE, las=1, xlim = c(0,15), ylim=c(-80,40), cex.axis=1.25)

      points(izh$sim$u[idx],izh$sim$v[idx], type='p', pch=16, cex=.65, col=Cols[colidx])
      mtext(expression(u[t]), side=1, line=3, cex=1.5)
      mtext(expression(v[t]), side=2, line=3, cex=1.5)

      izh_par = list(V0 = -70, U0 = 10, sig=20, a = 0.02, b = 0.2, c = -50, d = 2, I = 20, dt = izh$dt, type = "tonic bursting")


      # Add fixed points
      f <- function(v) 0.04*v^2+5*v+140+round(mean(izh$sim$i))
      g <- function(v) 0.08*v+5
      v = seq(-80,20,.1)
      idx.stable = which(g(v) < 0)
      idx.unstable = which(g(v) >= 0)

      lines(f(v)[idx.stable],v[idx.stable], ylim=c(-100,100), lwd=3, col=add.alpha('orange',.95))
      lines(f(v)[idx.unstable],v[idx.unstable], lty=2, lwd=3, col=add.alpha('orange',.95))

      izh_par2 = list(V0 = -70, U0 = 10, sig=0, a = 0.02, b = 0.2, c = -50, d = 2, I = 20, dt = .01*.1, type = "tonic bursting")
      izh2     = sim_izhikevich(500e3, izh_par2)
      u = izh2$u[-(1:10000)]
      v = izh2$v[-(1:10000)]
      id = which(diff(v) < -20)
      lines(u[1:id[1]],v[1:id[1]], lwd=3, col=add.alpha('black',.75))
      for(i in 2:length(id)){
        lines(u[(id[i-1]+1):id[i]],v[(id[i-1]+1):id[i]], lwd=3, col=add.alpha('black',.75))
      }

      # idx = seq(1,nsim,5)
      par(mar=c(0,0,2,1))
      # Colorbar
      plot(rep(1,1000),seq(0,1,len=1000), col=Cols, pch=15,cex=3, axes=FALSE, xlab=NA, ylab=NA)
      axis(side = 4, at = seq(0,1,.25), labels = paste0(seq(0,1,.25)*100,"%"), cex.axis=1.25, line = -2.5, lwd = 0, las=1)
      mtext(expression("P("~X[t]>0~")"), side = 3, line = -.2, cex=1.5)
  dev.off()

  glm.res$ks
  izh.res$ks
