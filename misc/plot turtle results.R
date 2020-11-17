misc::clean_up()
library(SSGLM)
datalib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/plots/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/TeX/SSGLM plos-latex-template/gfx/"

savePDF = FALSE



all.res = list(n1=list(),n2=list())
for(i in 1:2){
  for(j in 1:7){
    fn = paste0(datalib,"turtle_neuron",i,"_trial",j,".Rda")
    load(fn)
    all.res[[i]][[j]] = res
    cat("\nLoaded neuron",i,"trial",j)
  }
}
res = all.res


if(savePDF) pdf(file = paste0(plotlib, "turtle_decode.pdf"), height=5, width=8)
    par(mar=c(3,3,1,1), mfrow=c(1,1))
    plot(0,0, type='n', xlim=c(0,10*1e3), ylim=c(1,7), yaxt='n', bty='n', ann=FALSE, bty='n')
    axis(side = 2, at = 1:7, labels=7:1, tick = FALSE, las=1)

    offs = -.1
    axis(side = 2, at = 1:7+offs, labels = rep("(2)",7), tick = FALSE, las=1, cex.axis=.5, line = -1)
    offs = .1
    axis(side = 2, at = 1:7+offs, labels = rep("(1)",7), tick = FALSE, las=1, cex.axis=.5, line = -1)


    mtext("Trial (neuron)", 2,line=2)
    mtext("Time from stimulus (ms)", 1,line=2)
    for(i in 1:7){
      offs = .1
      x = res$n1[[i]]$spk
      tmp  = (8-i+offs)*(res$n1[[i]]$prob>0.5)
      tmp[tmp==0] = NA

      segments(0,8-i+offs,1e4,8-i+offs, col=add.alpha('black',.15), lwd=5, lend=1)
      lines(tmp, col=add.alpha('black',.6), lwd=5, lend=1)
      points(x,rep(8-i,length(x))+offs, pch=124, col=add.alpha('black',.75))

      offs = -.1
      x = res$n2[[i]]$spk

      tmp = (8-i+offs)*(res$n2[[i]]$prob>0.5)
      tmp[tmp==0] = NA
      segments(0,8-i+offs,1e4,8-i+offs, col=add.alpha('black',.15), lwd=5, lend=1)
      lines(tmp, col=add.alpha('black',.6), lwd=5, lend=1)
      points(x,rep(8-i,length(x))+offs, pch=124, col=add.alpha('black',.75))
    }
    abline(h=1:7, lty=3)
if(savePDF) dev.off()


if(savePDF) pdf(file = paste0(plotlib, "turtle_kernel.pdf"), height=3, width=5)
par(mar=c(3,3,1,1), mfrow=c(1,1))
    col1 = 'red'
    col2 = 'dodgerblue'

    maxk  = nrow(res$n1[[1]]$Sb)
    P1 = P2 = 0
    plot(0,0, type='n', xlim=c(0,maxk), ylim=c(0,.031), bty='n')

    for(i in 1:7){
      tmpres = res$n1[[i]]
      tmpP1  = as.numeric(tmpres$est$P)
      P1     = P1+tmpP1/7
      lines(exp(tmpres$est$Kb+tmpres$est$B0b), col=add.alpha(col1,.75), lwd=1, lty=3)

      tmpres = res$n2[[i]]
      tmpP2  = as.numeric(tmpres$est$P)
      P2     = P2 + tmpP2/7
      lines(exp(tmpres$est$Kb+tmpres$est$B0b), col=add.alpha(col2,.75), lwd=1, lty=3)
    }

    S1 = res$n1[[1]]$S
    S2 = res$n2[[1]]$S
    p1 = length(P1)
    p2 = length(P2)

    k1 = S1%*%P1[2:(p1-1)]
    k2 = S2%*%P2[2:(p2-1)]

    lines(exp(k1+P1[1]), col=add.alpha(col1,.85), lwd=2)
    lines(exp(k2+P2[1]), col=add.alpha(col2,.85), lwd=2)

    mtext("ms", 1,line=2)

    legend("topright", c("Neuron 1", "Neuron 2", "average", "individual"), col=add.alpha(c(col1,col2,'black','black'),.65), lty=c(NA,NA,1,3), lwd=c(NA,NA,2,1), pch=c(15,15,NA,NA),bty='n')
if(savePDF) dev.off()



if(savePDF) pdf(file = paste0(plotlib, "turtle_kernel_nobase.pdf"), height=3, width=5)
par(mar=c(3,3,1,1), mfrow=c(1,1))
    plot(0,0, type='n', xlim=c(0,maxk), ylim=c(0,2), bty='n')

    for(i in 1:7){
      tmpres = res$n1[[i]]
      lines(exp(tmpres$est$Kb), col=add.alpha(col1,.75), lwd=1, lty=3)
      tmpres = res$n2[[i]]
      lines(exp(tmpres$est$Kb), col=add.alpha(col2,.75), lwd=1, lty=3)
    }

    lines(exp(k1), col=add.alpha(col1,.85), lwd=2)
    lines(exp(k2), col=add.alpha(col2,.85), lwd=2)

    mtext("ms", 1,line=2)
if(savePDF) dev.off()



ks.results <- function(neuron,trial, P, k){
  res  = res[[neuron]][[trial]]
  p    = length(P)
  L    = get_intensity(as.numeric(res$y), res$X, B0.bst = P[1], k.bst = k, B0.rst = P[p],k.rst = NULL)
  z    = rescale_waiting_times(res$spk, L, 1)
  ks   = suppressWarnings(ks.test(z, 'pexp'))

  return(list(z=z, pval = ks$p.value))
}

z     = list(n1=numeric(),n2=numeric())
P     = list(P1, P2)
k     = list(k1, k2)
pvals = matrix(nr=7,nc=2)
for(i in 1:2){
  for(j in 1:7){
    tmp  = ks.results(i,j,P[[i]],k[[i]])
    z[[i]] = c(z[[i]],tmp$z)
    pvals[j,i] = tmp$pval
  }
}

round(pvals,3)

ks.test(z$n1, "pexp")
ks.test(z$n2, "pexp")

col1 = add.alpha('red',.75)
col2 = add.alpha('dodgerblue',.25)

if(savePDF) pdf(file = paste0(plotlib, "turtle_ks1.pdf"), height=4, width=4)
    ks_plot(z$n1, a = .05, col1 = col1, col2 = col2)
if(savePDF) dev.off()

if(savePDF) pdf(file = paste0(plotlib, "turtle_ks2.pdf"), height=4, width=4)
    ks_plot(z$n2, a = .05, col1 = col1, col2 = col2)
if(savePDF) dev.off()

length(z$n1)
length(z$n2)
round(exp(P1),3)
round(exp(P2),3)
296/163
1/round(exp(P2[1]-tmpP1[1]),3)




z2     = numeric()
i=2
  for(j in 1:7){
    if(j != 2){
      tmp  = ks.results(i,j,P[[i]],k[[i]])
      z2 = c(z2,tmp$z)
    }
  }
ks.test(z2, 'pexp')
