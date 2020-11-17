misc::clean_up()
library(SSGLM)
datalib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/plots/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/TeX/SSGLM plos-latex-template/gfx/"

savePDF = FALSE
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


# Plot the kernel estimates

col1 = add.alpha('red',.75)
col2 = add.alpha('red',.5)

col3 = add.alpha('black',.5)
col4 = add.alpha('black',.5)

layout(1)
 # pdf(paste0(plotlib,"sim_kernel_nobase.pdf"), width=5, height=3)


tmp = t(glm.res$P)
glm.P = cbind(tmp, tmp-1.96*sqrt(diag(glm.res$cov)), tmp+1.96*sqrt(diag(glm.res$cov)))
colnames(glm.P) = c("est","lo","hi")
glm.P = as.data.frame(glm.P)
glm.kb.est = glm$S%*%as.numeric(glm.P$est[2:4])
glm.kb.lo = glm$S%*%as.numeric(glm.P$lo[2:4])
glm.kb.hi = glm$S%*%as.numeric(glm.P$hi[2:4])

tmp = t(izh.res$P)
izh.P = cbind(tmp, tmp-1.96*sqrt(diag(izh.res$cov)), tmp+1.96*sqrt(diag(izh.res$cov)))
colnames(izh.P) = c("est","lo","hi")
izh.P = as.data.frame(izh.P)
izh.kb.est = izh$mpfpp$S%*%as.numeric(izh.P$est[2:4])
izh.kb.lo = izh$mpfpp$S%*%as.numeric(izh.P$lo[2:4])
izh.kb.hi = izh$mpfpp$S%*%as.numeric(izh.P$hi[2:4])


par(mfrow=c(2,1),mar=c(4,4,1,1), oma=c(0,0,0,0), las=1)
layout(1)
if(savePDF) pdf(paste0(plotlib,"sim_kernel_nobase.pdf"), width=5, height=3)
    par(mar=c(4,4,1,1), oma=c(0,0,0,0), las=1)
    plot(exp(glm.kb.est), type='n', ylim=c(0,3), bty='n', xlab="",ylab="", xaxt='n')
    # polygon( c(1:length(izh.kb.est), length(izh.kb.est):1), exp(c(izh.kb.lo, rev(izh.kb.hi))), col=col4, border=NA)
    polygon( c(1:nrow(glm$S), nrow(glm$S):1), exp(c(glm.kb.lo, rev(glm.kb.hi))), col=col2, border=NA)
    # lines(exp(glm.kb.est), lwd=2, col=col1)
    lines(exp(izh.kb.est), lwd=2, col=col3)
    axis(1, at = pretty(c(0,2000)), labels = pretty(c(0,20)))
    abline(h=1, lty=2, col=add.alpha('black',.75))
    mtext("ms",1, outer=FALSE, line=2, cex=1.25)
    # legend("topright", c("Izhikevich","SSGLM","With baseline", "Without baseline"), lty=c(NA,NA,1,3), col=c(col1,col2, col1,col1), pch=c(15,15,NA,NA), lwd=3, bty='n')
    legend("topright", c("Estimated", "True (Izhikevich)"), fill=c(col2,col4), bty='n')
    mtext("A", side = 3, line=0, adj=0, family = 'sans', cex=1.25)
if(savePDF) dev.off()


layout(1)
if(savePDF) pdf(paste0(plotlib,"sim_kernel_base.pdf"), width=5, height=3)
    par(mar=c(4,4,1,1), oma=c(0,0,0,0), las=1)
    plot(exp(glm.kb.est+glm.P$est[1]), type='n', ylim=c(0,10), bty='n', xlab="",ylab="", xaxt='n')
    # polygon( c(1:length(izh.kb.est), length(izh.kb.est):1), exp(c(izh.kb.lo+izh.P$lo[1], rev(izh.kb.hi+izh.P$hi[1]))), col=col4, border=NA)
    polygon( c(1:nrow(glm$S), nrow(glm$S):1), exp(c(glm.kb.lo+glm.P$lo[1], rev(glm.kb.hi+glm.P$hi[1]))), col=col2, border=NA)
    # lines(exp(glm.kb.est+glm.P$est[1]), lwd=2, col=col1)
    lines(exp(izh.kb.est+izh.P$est[1]), lwd=2, col=col3)
    axis(1, at = pretty(c(0,2000)), labels = pretty(c(0,20)))
    abline(h=1, lty=2, col=add.alpha('black',.75))
    mtext("ms",1, outer=FALSE, line=2, cex=1.25)
    # legend("topright", c("With baseline", "Without baseline"), lty=c(1,3), lwd=3, bty='n')
    legend("topright", c("Estimated", "True (Izhikevich)"), fill=c(col2,col4), bty='n')
    mtext("B", side = 3, line=0, adj=0, family = 'sans', cex=1.25)
if(savePDF) dev.off()
