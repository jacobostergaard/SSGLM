misc::clean_up()
library(SSGLM)
datalib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/plots/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/TeX/SSGLM plos-latex-template/gfx/"
set.seed(1234)

# Load turtle spike trians
turtles <- read.delim(paste0(datalib,"turtles.txt"))
names(turtles) = paste0("s",1:50)

# Convert spike times to seconds
tosec   = 2.5e-5
turtles = turtles*tosec

# Pick neurons to analyze
ids = 12:13
get_spiketimes <- function(i){
  out = list()
  tmp = turtles[,i]
  tmp = tmp[tmp>0]    # remove trailing zeros
  for(j in 1:10){
    x = tmp[tmp>40*(j-1) & tmp<40*j]-(j-1)*40
    out[[j]] = x
  }
  return(out)
}


trials = list()
for(i in 1:ncol(turtles)){
  trials[[i]] = get_spiketimes(i)
}
names(trials) = paste0("n",1:ncol(turtles))

spks_in_trials = matrix(unlist(lapply(trials,function(x) unlist(lapply(x,length)))), nr=50,nc=10, byrow=TRUE)

# quantile(spks_in_trials)

# apply(spks_in_trials, 2, function(x) sum(x<10))

# quantile(apply(spks_in_trials[,-9], 1, mean))
tmp = apply(spks_in_trials[,-9], 1, mean)
# hist(tmp, breaks=c(seq(0,500,10)))

# quantile(tmp,c(.25,.75))

idx = which(tmp<90 & tmp>25)
# sum(spks_in_trials[idx,-9] < 20)

# hi = 84
# lo = 19
# idx = 1*(spks_in_trials >= lo & spks_in_trials <= hi)

# tmp = apply(spks_in_trials,1,mean)
# tmplo = apply(spks_in_trials,1,min)
# tmphi = apply(spks_in_trials,1,max)
# boxplot(spks_in_trials)
# boxplot(spks_in_trials[tmp<100,])
# boxplot(spks_in_trials[tmp<100 & tmp>10,-9])
#
# spks_in_trials
# idx = which(tmphi<200 & tmplo>0)

# boxplot(spks_in_trials[idx,-9])

# idx
# length(idx)
plot(0,0, type='n', xlim=c(0,40), ylim=c(0,50))
for(i in 1:length(idx)){
  for(j in 1:10){
    pts = trials[[idx[i]]][[j]]
    n = length(pts)

    if(j!=9){
      points(pts,rep(idx[i]+j/10,n), pch='.', cex=1.5, col = ifelse(j==9, "red","black"))
    }

  }
}

trials = trials[idx]
# length(trials)

# trials = list(n1=get_spiketimes(12),n2=get_spiketimes(13))
# rmv = c(3,8,9) # remove these trials: either little to no spikes, or very irregular behavior
# trials$n1 = trials$n1[-rmv]
# trials$n2 = trials$n2[-rmv]



# lapply(trials,length)
# trials$n1

# subtract initial 10s quiescent period (stimulus at 10s) and convert to milliseconds
# trials$n1 = lapply(trials$n1, function(x) (x-10)*1e3)
# trials$n2 = lapply(trials$n2, function(x) (x-10)*1e3)


# trials can be rounded to nearest millisecond without losing any spikes, i.e. same number of unique spikes after rounding

# Downsample measurements
downsample <- function(x, freq=1){
  x = round(x/freq)

  if(any(diff(x)==0)){
    # remove any double counts
    x = x[-which(diff(x) == 0)]
  }
  return(x)
}

# trials$n1 = lapply(trials$n1,downsample)
# trials$n2 = lapply(trials$n2,downsample)


for(i in 1:length(trials)){
  trials[[i]] = lapply(trials[[i]], function(x) x*1e3)
  # trials[[i]] = lapply(trials[[i]], function(x) (x-10)*1e3)

  trials[[i]] = lapply(trials[[i]],downsample)
}

isi = unlist( lapply(trials, function(x) lapply(x,function(y) diff(sort(y)))) )

# quantile(isi)
# hist(isi, breaks=1000, xlim=c(0,1000))
# plot(isi, pch=16, col=add.alpha('red',.5))
# abline(h=250, lty=3)

maxb = 250 # choose max isi for spikes in a burst

# Find all burst lenghts for all trials for both neurons
bst = numeric(0)
nneurons = length(trials)
ntrials = length(trials[[1]])
for(i in 1:nneurons){
  for(j in 1:ntrials){
    # i=j=1
    tmp_spk = trials[[i]][[j]]

    tmp = burst_info(tmp_spk,cutoff = maxb, refrac = min(isi))
    bst = c(bst,tmp$burstlength)
  }
}

bst = bst[bst!=(1+min(isi))]
rst = as.numeric(isi[isi>maxb])

# Find X transition parameterrs using all isi's
bst.h = hist(bst, breaks=seq(.5,max(bst)+.5,1), plot=FALSE)
rst.h = hist(rst, breaks=seq(.5,max(rst)+.5,1), plot=FALSE)

init = c(60e3,1,100e3, 1)
bst.pars = optim_xpars(init[1:2],bst.h)
rst.pars = optim_xpars(init[3:4],rst.h)

tmp = list()
tmp$pars = data.frame(bst = bst.pars, rst = rst.pars, row.names = c("alpha","beta"))
tmp$bst = getF(bst.h)
tmp$rst = getF(rst.h)
xpars = tmp$pars
# xpars

col1 = add.alpha('black',.25)
col2 = add.alpha('red',.75)
col3 = add.alpha('black',.75)

idx = 1:10000

fn = paste0(plotlib, "turtle_xfit_bst.pdf")
# pdf(fn, height=3.5, width=4)
par(mar=c(2,2,2,2), oma=c(1,1,1,1), bty='n', mfrow=c(2,1))
plot(0,0, type='n', xlim=c(0,3000), ylim=c(0,1))
lines(idx,getF(bst.h)$y[idx], col=col3, lty=1, lwd=1)
points(idx,getF(bst.h)$y[idx], col=col1, cex=.5, pch=16)
lines(idx,Z.cdf(idx, bst.pars[1],bst.pars[2])$cdf, col=col2, lwd=3)
legend("topleft", bty='n', legend=c("Observed","Fitted"), pch=c(16,NA), lty=c(NA,1), lwd=c(NA,3), col=c(col3,col2), seg.len = .85)
mtext("CDF", side = 2, line=2)
mtext("Burst length (ms)", side=1, line=2)
# dev.off()

fn = paste0(plotlib, "turtle_xfit_rst.pdf")
# pdf(fn, height=3.5, width=4)
# par(mar=c(2,2,2,2), oma=c(1,1,1,1), bty='n', mfrow=c(1,1))
plot(0,0, type='n', xlim=c(0,4000), ylim=c(0,1))
lines(idx,getF(rst.h)$y[idx], col=col3, lty=1, lwd=1)
points(idx,getF(rst.h)$y[idx], col=col1, cex=.5, pch=16, xlim=c(1000,4000))
lines(idx,Z.cdf(idx, rst.pars[1],rst.pars[2])$cdf, col=col2, lwd=3)
legend("topleft", bty='n', legend=c("Observed","Fitted"), pch=c(16,NA), lty=c(NA,1), lwd=c(NA,3), col=c(col3,col2), seg.len = .95)
mtext("CDF", side = 2, line=2)
mtext("Rest length (ms)", side=1, line=2)
# dev.off()


# Find quantiles of bursting isi's
g <- function(spk, maxb=250, prbs = NULL){
  if(is.null(prbs)){
    prbs = c(.1,.25,.5,.75,.9)
  }
  isi       = diff(spk)
  qts       = quantile(isi[isi<maxb], probs = prbs)

  return(qts)
}

qts = list()
for(i in 1:nneurons){
  tmp = matrix(nr=ntrials, nc=5)
  for(j in 1:ntrials){
    tmp[j,] = g(trials[[i]][[j]],maxb)
  }
  colnames(tmp) = c(paste0(100*c(.1,.25,.5,.75,.9),"%"))
  qts[[i]] = tmp
}

names(qts) = paste0("n",1:nneurons)

tmp = lapply(qts, function(x) apply(x, 2, function(y) mean(y,na.rm = TRUE)))
names(tmp) = paste0(names(qts),".mean")

qts = c(qts,tmp)


f <- function(neuron, trial){
  # neuron = 1; trial = 1;

  kts = c(25,50,75,100);
  kts = c(20,40,60,80,100,120);
  Bb0 = 0; Br0 = 0; verbose = TRUE
  spk = trials[[neuron]][[trial]]

  # kts = quantile(diff(spk), prob=c(0.1,.25))

  # kts = quantile(diff(spk), prob=0.1)

  if(all(kts<=1)){
    isi       = diff(spk)
    kts       = quantile(isi[isi<maxb], probs = kts)
  }


  kts       = c(kts,maxb)    # add auxiliary knot at max burst 250 ms
  kts       = c(0,0,0,0,kts) # add auxiliary knot at 0 ms
  Sb        = splines::splineDesign(kts, x=1:max(kts), ord=4, outer.ok=TRUE)
  Bb        = rep(0,ncol(Sb))
  Bb[1]     = -10 # First weight is (large) negative to ensure refractoriness
  Br        = Bb[1]
  Sr        = splines::splineDesign(c(rep(0,4),max(kts)), x=1:max(kts), ord=4, outer.ok=TRUE)
  maxk      = max(kts)
  Sr        = 1-log_logistic(1:(5*maxk), maxk, 10)
  # Sr        = matrix(Sr[Sr>1e-6], ncol=1)

  M         = 200
  pars      = SSGLM::make_mpfpp_par_input(Bb0, Bb, Sb, Br0, Br=Br, Sr, xpars, M=M, q0 = 1e-8, q1 = 0, s_lag = 2000, dt=1)



  pars$q    = c(1e-5,rep(1e-4,ncol(Sb)), 1e-5, rep(1e-4,ncol(Sr)))
  y         = as.numeric(SSGLM::make_y(spk,1,N = 40e3, verbose = FALSE)$y)
  # y = y[1:13800]
  X         = sample(c(-1,1),M,replace=TRUE) # Sample random initial states +/- 1 for X

  # Run algorithm

  res = SSGLM::mpfpp(y,X,pars, usetrue = FALSE, verbose = verbose)
  res$y = y
  res$est = get_est(res)

  return(res)
}

# length(trials)
totres = list()


totres = list()

i=j=1
# res = f(3,1)
# for(i in 1:length(trials)){
# for(i in 9:length(trials)){
  # for(i in 1:2){
  tmpres = list()
  cat("Analyzing neuron", i,"\n")
  # for(j in 1:10){
    # for(j in 1:2){
    cat("  - trial",j,"\n")
    if(j != 9 & length(trials[[i]][[j]]) > 10){ # exclude trials with less than 10 spikes
      tmpres[[j]] = f(i,j)
    } else{
      tmpres[[j]] = NA
    }
  # }
  totres[[i]] = tmpres
  # fn = paste0(datalib,"turtle_",names(trials)[i],".Rda")
  # save(tmpres, file = fn)
# }




par(mar=c(3,3,1,1), mfrow=c(1,1))

layout(matrix(c(1,1,2,3),nc=2))
plot(0,0, type='n', xlim=c(10*1e3,20*1e3), ylim=c(1,10), yaxt='n', bty='n', ann=FALSE, bty='n')
axis(side = 2, at = 1:10, labels=10:1, tick = FALSE, las=1)

# offs = -.1
# axis(side = 2, at = 1:7+offs, labels = rep("(2)",7), tick = FALSE, las=1, cex.axis=.5, line = -1)
# offs = .1
# axis(side = 2, at = 1:7+offs, labels = rep("(1)",7), tick = FALSE, las=1, cex.axis=.5, line = -1)

mtext("Trial (neuron)", 2,line=2)
mtext("Time from stimulus (ms)", 1,line=2)
for(i in 1:10){
  res = totres[[i]]
  if(!is.null(res)){
    x = res$spk
    tmp = (11-i)*(res$prob>0.5)
    tmp[tmp==0] = NA

    segments(0,11-i,2e4,11-i, col=add.alpha('black',.15), lwd=5, lend=1)
    lines(tmp, col=add.alpha('black',.6), lwd=5, lend=1)
    points(x,rep(11-i,length(x)), pch=124, col=add.alpha('black',.75))

    cat("Trial",i,"has pval", round(res$est$ks$p.value,3))
  } else{
    cat("Trial",i,"has pval NULL")
  }
  cat("\n")
}
abline(h=1:10, lty=3)


Pest = numeric(0)
kb = matrix(NA, nr=10, nc=nrow(res$Sb))
kr = matrix(NA, nr=10, nc=nrow(res$Sr))
for(i in 1:10){
  if(i!=9){
    res = totres[[i]]

    tmp = res$est$P
    Pest = rbind(Pest,tmp)
    np = length(tmp)
    kb[i,] = exp(res$Sb%*%as.numeric(res$est$P[2:(np-2)])+res$est$P$b0)
    kr[i,] = exp(res$Sr%*%as.numeric(res$est$P[np])+res$est$P$r0)
  }
}



# par(mfrow=c(2,1))
plot(0,0,  type='n', xlim=c(0,300), ylim=c(0,.2))
for(i in 1:10){
  lines(kb[i,])
}
tmp = apply(Pest,2,mean)
lines(exp(res$Sb%*%as.numeric(tmp[2:(np-2)])+tmp[1]),lwd=2, col='red')

plot(0,0,  type='n', xlim=c(0,300), ylim=c(0,0.003))
for(i in 1:10){
  lines(kr[i,])
}
lines(exp(res$Sr%*%as.numeric(tmp[np])+tmp[np-1]),lwd=2, col='red')

Pest
