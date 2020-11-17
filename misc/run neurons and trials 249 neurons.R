library(misc)
library(SSGLM)
library(parallel)
clean_up()
datalib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"
set.seed(1234)


#### Scrip below ####


# Load turtle spike trians
dat = R.matlab::readMat(icloud_lib("GitHub/Source/R/Extracellular triggered EPSP/ExtracellularUnits.mat")) #okay
turtles = list()
for(i in 1:249)
  eval(parse( text=paste0("turtles$n",i," <- dat$units[[1]][[i]][[1]]") ))


ntrains = length(turtles)

# Convert spike times to seconds
tosec   = 2.5e-5
for(i in 1:ntrains)
  turtles[[i]] = turtles[[i]]*tosec

# Pick neurons to analyze
get_spiketimes <- function(i){
  out = list()
  tmp = turtles[[i]]
  tmp = tmp[tmp>0]    # remove trailing zeros
  for(j in 1:10){
    x = tmp[tmp>40*(j-1) & tmp<40*j]-(j-1)*40
    out[[j]] = x
  }
  return(out)
}


trials = list()
for(i in 1:ntrains){
  trials[[i]] = get_spiketimes(i)
}
names(trials) = paste0("n",1:ntrains)


spks_in_trials = matrix(unlist(lapply(trials,function(x) unlist(lapply(x,length)))), nr=ntrains,nc=10, byrow=TRUE)

isi = unlist(lapply(trials, function(x) lapply(x, diff)))

tmp = apply(spks_in_trials[,-9], 1, mean)
# quantile(tmp)

idx = which(tmp<100 & tmp>25) # roughly 25% and 75% quantiles
trials = trials[idx]
# length(trials)



# trials can be rounded to nearest millisecond without losing any spikes, i.e. same number of unique spikes after rounding

# Downsample measurements
downsample <- function(x, freq=1e-3){
  x = round(x/freq)

  if(any(diff(x)==0)){
    # remove any double counts
    x = x[-which(diff(x) == 0)]
  }
  return(x)
}

for(i in 1:length(trials)){
  trials[[i]] = lapply(trials[[i]],downsample)
}

isi = unlist( lapply(trials, function(x) lapply(x,function(y) diff(sort(y)))) )


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

# # fn = paste0(plotlib, "turtle_xfit_bst.pdf")
# # pdf(fn, height=3.5, width=4)
# par(mar=c(2,2,2,2), oma=c(1,1,1,1), bty='n', mfrow=c(2,1))
# plot(0,0, type='n', xlim=c(0,3000), ylim=c(0,1))
# lines(idx,getF(bst.h)$y[idx], col=col3, lty=1, lwd=1)
# points(idx,getF(bst.h)$y[idx], col=col1, cex=.5, pch=16)
# lines(idx,Z.cdf(idx, bst.pars[1],bst.pars[2])$cdf, col=col2, lwd=3)
# legend("topleft", bty='n', legend=c("Observed","Fitted"), pch=c(16,NA), lty=c(NA,1), lwd=c(NA,3), col=c(col3,col2), seg.len = .85)
# mtext("CDF", side = 2, line=2)
# mtext("Burst length (ms)", side=1, line=2)
# # dev.off()

# # fn = paste0(plotlib, "turtle_xfit_rst.pdf")
# # pdf(fn, height=3.5, width=4)
# # par(mar=c(2,2,2,2), oma=c(1,1,1,1), bty='n', mfrow=c(1,1))
# plot(0,0, type='n', xlim=c(0,4000), ylim=c(0,1))
# lines(idx,getF(rst.h)$y[idx], col=col3, lty=1, lwd=1)
# points(idx,getF(rst.h)$y[idx], col=col1, cex=.5, pch=16, xlim=c(1000,4000))
# lines(idx,Z.cdf(idx, rst.pars[1],rst.pars[2])$cdf, col=col2, lwd=3)
# legend("topleft", bty='n', legend=c("Observed","Fitted"), pch=c(16,NA), lty=c(NA,1), lwd=c(NA,3), col=c(col3,col2), seg.len = .95)
# mtext("CDF", side = 2, line=2)
# mtext("Rest length (ms)", side=1, line=2)
# # dev.off()


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


# unlist(lapply(trials, function(x) unlist(lapply(x,diff))))
isi250 = isi[isi < 250]
# median(isi250)
isi50 = isi[isi < 50] # 50 is about the median of the ISIs < 250

# hist(isi250, breaks=100)
# abline(v=quantile(isi250, probs = c(.1,.25,.5)), col='red')

# quantile(isi250, probs = c(.1,.25,.5))
# round(ecdf(isi250)(c(20,40,60,80,100,120)),2)

# layout(1:2)
# kts = c(quantile(isi50,probs=c(.1,.25,.5,.75)),50)     # overall knots by probs
# hist(isi50, breaks=50)
# seq(0,50,length=5)
# kts       = c(0,0,0,0,kts) # add auxiliary knot at 0 ms
# Sb        = splines::splineDesign(kts, x=1:max(kts), ord=4, outer.ok=TRUE)
# matplot(Sb, lty=1, col=1, type='l', xlim=c(0,50))
#
#
# hist(isi50, breaks=50)
# Ps = c(-10,0,2,2,0)
# ktmp = Sb%*%Ps
# plot(exp(ktmp))
# abline(v=kts)


f <- function(neuron, trial){
  # neuron = 1; trial = 1;

  spk = trials[[neuron]][[trial]]
  if(length(spk)>10){
        Bb0 = 0; Br0 = 0; verbose = FALSE

        # kts = c(25,50,75,100);
        # kts = c(20,40,60,80,100,120);                       # specific knots
        # prbs = ecdf(isi[isi < 250])(kts)

        # kts = quantile(diff(spk), prob=c(0.1,.25,.5))       # individual knots by probs
        # kts = quantile(isi[isi < 250],probs=c(.1,.25,.5))     # overall knots by probs

        # isi50 = isi[isi < 50]
        kts = c(quantile(isi50,probs=c(.1,.25,.5)),50)     # overall knots by probs
        # kts = quantile(isi[isi < 250],probs=c(.1,.25,.5,.75)) # overall knots by probs
        # kts = quantile(diff(spk), prob=0.1)
        kts = c(0,0,0,seq(0,50,length=6))


        # if(all(kts<=1)){
        #   isi       = diff(spk)
        #   kts       = quantile(isi[isi<maxb], probs = kts)
        # }


        # kts       = c(kts,maxb)    # add auxiliary knot at max burst 250 ms



        # kts       = c(0,0,0,0,kts) # add auxiliary knot at 0 ms


        Sb        = splines::splineDesign(kts, x=1:max(kts), ord=4, outer.ok=TRUE)
        Bb        = rep(0,ncol(Sb))
        Bb[1]     = -10 # First weight is (large) negative to ensure refractoriness
        Br        = Bb[1]
        Sr        = splines::splineDesign(c(rep(0,4),max(kts)), x=1:max(kts), ord=4, outer.ok=TRUE)
        maxk      = max(kts)
        Sr        = 1-log_logistic(1:(3*maxb), maxb, 10)
        # Sr        = matrix(Sr[Sr>1e-6], ncol=1)
        # matplot(Sb, lty=1, col=1, type='l')
        # matplot(Sr, lty=1, col=1, type='l')

        M         = 200
        pars      = SSGLM::make_mpfpp_par_input(Bb0, Bb, Sb, Br0, Br=Br, Sr, xpars, M=M, q0 = 1e-8, q1 = 0, s_lag = 2000, dt=1)



        pars$q    = c(1e-5,rep(1e-4,ncol(Sb)), 1e-5, rep(1e-4,ncol(Sr)))
        y         = as.numeric(SSGLM::make_y(spk,1,N = 40e3, verbose = FALSE)$y)
        # y = y[1:13800]
        X         = sample(c(-1,1),M,replace=TRUE) # Sample random initial states +/- 1 for X

        # Run algorithm

        # y = y[1:5000]
        res = SSGLM::mpfpp(y,X,pars, usetrue = FALSE, verbose = verbose)
        res$y = y
        res$est = get_est(res)
  } else{
    res = NA
  }
  return(res)
}

# length(trials)
totres = list()


i=j=1
# res = f(3,1)
for(i in 1:length(trials)){
  # tmpres = list()
  cat("Analyzing neuron", i,"\n")
  # for(j in 1:10){
  # # for(j in 1:2){
  #   cat("  - trial",j,"\n")
  #   if(j != 9 & length(trials[[i]][[j]]) > 10){ # exclude trials with less than 10 spikes
  #     tmpres[[j]] = f(i,j)
  #   } else{
  #     tmpres[[j]] = NA
  #   }
  # }
  # totres[[i]] = tmpres
  totres[[i]] = lapply(1:10, function(j) f(neuron=i,trial=j))
  # totres[[i]] = mclapply(1:2, function(j) f(neuron=i,trial=j), mc.cores = 7)

  fn = paste0(datalib,"turtle2_",names(trials)[i],"_qts3.Rda")
  save(tmpres, file = fn)
}



rest = mclapply(1:3, g, mc.cores = 3)


i
g(1)
i=1
rest = list()
rest[[1]] = g(1)
rest[[2]] = g(2)


cat("Analyzing neuron", i,"\n")
for(j in 1:10){
  # for(j in 1:2){
  cat("  - trial",j,"\n")
  if(j != 9 & length(trials[[i]][[j]]) > 10){ # exclude trials with less than 10 spikes
    tmpres[[j]] = f(i,j)
  } else{
    tmpres[[j]] = NA
  }
}



# spk = tmpres[[1]]$spk
#
# plot(spk, rep(1,length(spk)), pch=124)
# tmpprob = tmpres[[1]]$prob>.5
# tmpprob[tmpprob == 0] = NA
# lines(tmpprob, col=add.alpha('orange',.5), lwd=25, lend=3)


