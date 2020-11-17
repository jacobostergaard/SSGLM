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

col1 = add.alpha('black',.75)
col2 = add.alpha('red',.75)
col3 = add.alpha('dodgerblue',.25)
fn = paste0(plotlib, "izh_glm_res.pdf")
pdf(fn, height=12, width=9)

    par(mar=c(3,3,2,2), oma=c(1,1.5,1,1), bty='n', mfrow=c(3,2))

    # Plot fit of transition parameters
      plot(0,0, type='n', xlim=c(0,3000), ylim=c(0,1), xaxt='n')
          axis(1, at=pretty(c(0,3000)), pretty(c(0,30)))
          points(idx,getF(bst.h)$y[idx], col=col1, cex=.5, pch=16)
          lines(idx,Z.cdf(idx, bst.pars[1],bst.pars[2])$cdf, col=col2, lwd=3)
          legend("topleft", bty='n', legend=c("Observed","Fitted CDF"), pch=c(16,NA), lty=c(NA,1), pt.cex = c(1.25,NA), lwd=c(NA,3), col=c(col1,col2))
          mtext("CDF", side = 2, line=2)
          mtext("Burst length [ms]", side=1, line=2.25)
          mtext("A", side = 3, line=0, adj=0, family = 'sans', cex=1.25)

      plot(0,0, type='n', xlim=c(1000,4000), ylim=c(0,1), xaxt='n')
          axis(1, at=pretty(c(1000,4000)), pretty(c(10,40)))
          points(idx,getF(rst.h)$y[idx], col=col1, cex=.5, pch=16, xlim=c(1000,4000))
          lines(idx,Z.cdf(idx, rst.pars[1],rst.pars[2])$cdf, col=col2, lwd=3)
          legend("topleft", bty='n', legend=c("Observed","Fitted CDF"), pch=c(16,NA), lty=c(NA,1), pt.cex = c(1.25,NA), lwd=c(NA,3), col=c(col1,col2))
          mtext("CDF", side = 2, line=2)
          mtext("Rest length [ms]", side=1, line=2.25)
          mtext("B", side = 3, line=0, adj=0, family = 'sans', cex=1.25)

    # Plot the kernel estimates

      plot(exp(glm.res$Kb), type='l', bty='n', lwd=2, ylim=c(0,2.5), col=col2, xlab="",ylab="", xaxt='n')
          lines(exp(izh.res$Kb),col=col1,lwd=2, lty=1)
          axis(1, at = pretty(c(0,2000)), labels = pretty(c(0,20)))
          abline(h=1, lty=3, col=add.alpha('black',.5))
          mtext(expression(e^{S~beta}), side = 2, line=2)
          mtext("History lag [ms]",1, outer=FALSE, line=2.25, cex=1)
          legend("topright", c("Estimated", "True"), lty=c(1,1), col=c(col2,col1), lwd=3, bty='n')
          mtext("C", side = 3, line=0, adj=0, family = 'sans', cex=1.25)

      plot(exp(izh.res$Kb+izh.res$B0b), type='l', bty='n', lwd=2, ylim=c(0,7), col=col1, xlab="",ylab="", xaxt='n')
          lines(exp(glm.res$Kb+glm.res$B0b),col=col2,lwd=2, lty=1)
          axis(1, at = pretty(c(0,2000)), labels = pretty(c(0,20)))
          abline(h=1, lty=3, col=add.alpha('black',.5))
          mtext(expression(e^{S~beta+beta[0]}), side = 2, line=2)
          mtext("History lag [ms]",1, outer=FALSE, line=2.25, cex=1)
          legend("topright", c("Estimated", "True"), lty=c(1,1), col=c(col2,col1), lwd=3, bty='n')
          mtext("D", side = 3, line=0, adj=0, family = 'sans', cex=1.25)

    # Plot the KS tests

      ks_plot(glm.res$z, col1 = col2, col2 = col3)
          legend("topleft", c("Identity line", "Quantiles","95% CI limit"), lty=c(1,1,1), col=c(col1,col2,col3), lwd=c(2,2,8), bty='n')
          mtext("E", side = 3, line=0, adj=0, family = 'sans', cex=1.25)

      ks_plot(izh.res$z, col1 = col2, col2 = col3)
          legend("topleft", c("Identity line", "Quantiles","95% CI limit"), lty=c(1,1,1), col=c(col1,col2,col3), lwd=c(2,2,8), bty='n')
          mtext("F", side = 3, line=0, adj=0, family = 'sans', cex=1.25)

dev.off()
