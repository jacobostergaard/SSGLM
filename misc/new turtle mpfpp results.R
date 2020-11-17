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

quantile(spks_in_trials)
# hi = 84
# lo = 19
# idx = 1*(spks_in_trials >= lo & spks_in_trials <= hi)

tmp = apply(spks_in_trials,1,mean)
tmplo = apply(spks_in_trials,1,min)
tmphi = apply(spks_in_trials,1,max)
boxplot(spks_in_trials)
boxplot(spks_in_trials[tmp<100,])
boxplot(spks_in_trials[tmp<100 & tmp>10,-9])

spks_in_trials
idx = which(tmphi<300 & tmplo>0)
boxplot(spks_in_trials[idx,-9])

idx
length(idx)
plot(0,0, type='n', xlim=c(0,40), ylim=c(0,50))
for(i in 1:length(idx)){
  for(j in 1:10){
    pts = trials[[idx[i]]][[j]]
    n = length(pts)

    points(pts,rep(idx[i]+j/10,n), pch='.', cex=1.5, col = ifelse(j==9, "red","black"))
  }
}


trials[[20]]
trials[[34]]

trials = trials[idx]
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

quantile(isi)
hist(isi, breaks=1000, xlim=c(0,1000))
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
xpars

col1 = add.alpha('black',.25)
col2 = add.alpha('red',.75)
col3 = add.alpha('black',.75)

idx = 1:10000

fn = paste0(plotlib, "turtle_xfit_bst.pdf")
# pdf(fn, height=3.5, width=4)
par(mar=c(2,2,2,2), oma=c(1,1,1,1), bty='n', mfrow=c(1,1))
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
par(mar=c(2,2,2,2), oma=c(1,1,1,1), bty='n', mfrow=c(1,1))
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


# On to particle filter analysis

# Try a spiketrain spk
f <- function(neuron, trial, kts=.5, Bb0=0, Br0=0, maxb=250, verbose=FALSE, save_data = TRUE){
  # neuron = 1; trial = 1; kts = 52; Bb0 = -5; Br0 = -5; verbose = TRUE
  spk = trials[[neuron]][[trial]]
  if(all(kts<=1)){
    isi       = diff(spk)
    kts       = quantile(isi[isi<maxb], probs = kts)
  }


  kts       = c(kts,maxb)    # add auxiliary knot at max burst 250 ms
  kts       = c(0,0,0,0,kts) # add auxiliary knot at 0 ms
  Sb        = splines::splineDesign(kts, x=1:max(kts), ord=4, outer.ok=TRUE)
  Bb        = rep(0,ncol(Sb))
  Bb[1]     = -10 # First weight is (large) negative to ensure refractoriness
  Sr        = splines::splineDesign(c(rep(0,4),max(kts)), x=1:max(kts), ord=4, outer.ok=TRUE)
  M         = 200
  pars      = SSGLM::make_mpfpp_par_input(Bb0, Bb, Sb, Br0, Br=Bb[1], Sr, xpars, M=M, q0 = 1e-8, q1 = 0, s_lag = 2000, dt=1)

  pars$q    = c(1e-5,rep(1e-4,ncol(Sb)), 1e-5, rep(1e-4,ncol(Sr)))
  y         = as.numeric(SSGLM::make_y(spk,1,N = 10e3, verbose = FALSE)$y)
  X         = sample(c(-1,1),M,replace=TRUE) # Sample random initial states +/- 1 for X

  # Run algorithm

  res = SSGLM::mpfpp(y,X,pars, usetrue = FALSE, verbose = verbose)
  res$y = y
  res$est = get_est(res)

  fn = paste0(datalib,"turtle_",names(trials)[i],"_t",j,".Rda")
  if(save_data) save(res,file = fn)

  invisible(res)
}


# i=1;j=1
# res = f(neuron = i, trial = j, kts=50, Br0=-15, Bb0=-5, verbose=TRUE)
# print(res$est$ks)
#
# layout(matrix(c(1,1,2,3),nr=2,nc=2, byrow=TRUE))
# idx = seq(1,1e4,length=5000)
# plot(idx,res$X[idx], type='l', col=add.alpha('dodgerblue',.75),lwd=2)
# points(trials[[i]][[j]],rep(0,length(trials[[i]][[j]])), pch=124)
# par(new=TRUE)
# plot(idx,res$prob[idx], type='l', col=add.alpha('red',.75), ylim=c(0,1), axes=FALSE)
#
# np = length(res$est$P)
# idxb = grepl("b",names(res$est$P)) & !grepl("0",names(res$est$P))
# idxr = grepl("r",names(res$est$P)) & !grepl("0",names(res$est$P))
# k.tmp = res$Sb%*%as.numeric(res$est$P[idxb])
# plot(exp(k.tmp), type='l')
# k.tmp = res$Sr%*%as.numeric(res$est$P[idxr])
# lines(exp(k.tmp), col='red')
# ks_plot(res$est$z, col1 = add.alpha('red',.75), col2=add.alpha('dodgerblue',.25))
#
# Sys.sleep(1)

# Use 50% quantiles for bursting isi's as kts
# qts.50 = c(qts$n1.mean[3],qts$n2.mean[3])
qts.50 = lapply(qts[grepl("mean",names(qts))], function(x) x[3])
# quantile(unlist(qts.50))

i=j=1
# Run all trials
for(i in 1:nneurons){

  # kts = as.numeric(qts.50[i])
  # kts = c(25,50,75)
  kts = c(25,75)

  for(j in 1:ntrials){

    cat("Neuron",names(trials)[i],"- trial",j,"\n")
    res.pval = 0
    niter = 0
    # if(j != c(7,8,9)){
    if(j > 0){
      # while(res.pval < 0.05 & niter < 10){
        res = f(neuron = i, trial = j, kts = kts, Bb0 = -5, Br0 = -5, verbose = TRUE)
        res.pval = res$est$ks$p.value
        # niter = niter+1
        cat("\n")
      # }
    }

    print(res$est$ks)
    cat("-----------------------------------------------\n")
  }
}


