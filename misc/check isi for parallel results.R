idxn = 1:249
spks = trials[idxn]
Hst = numeric(0)
tmp.isi = lapply(spks, function(x) 1e3*unlist(lapply(x,diff)))
for(i in 1:249){
  tmp.hst = hist(tmp.isi[[i]], breaks=c(seq(0,100,1),1e6), plot=FALSE)
  x = tmp.hst$mids
  n = length(x)
  x = x[-n]
  y = tmp.hst$density[-n]
  Hst = rbind(Hst,y)
}

tick.at = pretty(1:ncol(Hst))
boxplot(Hst, xaxt='n', xlim=c(0,ncol(Hst)), ylim=c(0,.02))
axis(1, at=tick.at, seq(0,100,length=length(tick.at)))

par(new=TRUE)
hist(1e3*isi.all, xlim=c(0,100), breaks=c(seq(0,200,1)-.5,1e6), border=NA, col=add.alpha('red',.2), ylim=c(0,.02))


kts = c(0,0,0,seq(0,100,length=5))
# kts = c(0,0,0,0,10,15,20,20,20,20)

Sb        = splines::splineDesign(kts, x=1:max(kts), ord=4, outer.ok=TRUE)
par(new=TRUE)
matplot(Sb, lty=1, col=1, type='l', xlim=c(0,100), yaxt='n')
ncol(Sb)

B = c(-10,.5,2,-2)
K = exp(Sb%*%B)
abline(h=.8/max(K), lty=3)
K = K*.8/max(K)
lines(K, col='red', lwd=2)


M = (Ps[,1,])
ord = order(apply(M,1,na.median),decreasing = TRUE)

par.idx = 8
dimnames(Ps)[[2]][par.idx]
M = (Ps[,par.idx,])
# ord = order(apply(M,1,na.median),decreasing = TRUE)
boxplot(t(M[ord,]), xaxt='n')#, ylim=c(-10,.1))
dim(Ps)
dimnames(Ps)[[2]]


