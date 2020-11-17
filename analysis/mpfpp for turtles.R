library(misc)
library(SSGLM)
library(parallel)
clean_up()
datalib = "/Users/jacob/Data/" #"Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/TeX/SSGLM plos-latex-template/gfx/"

set.seed(1234)

# list.files(datalib)

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
# boxplot(spks_in_trials)

# tmp = which(spks_in_trials<10,arr.ind = TRUE)
# idx = sort(unique( tmp[tmp[,2]!=9,1] )) # neurons with a trial with >500 spikes
# spks_in_trials[idx,]
# idx = which(tmp<100 & tmp>25) # roughly 25% and 75% quantiles
idx = 1:249
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

# fn = paste0(plotlib, "turtle_xfit_bst.pdf")
#     pdf(fn, height=3.5, width=4)
#     par(mar=c(2,2,2,2), oma=c(1,1,1,1), bty='n', mfrow=c(1,1))
#     plot(0,0, type='n', xlim=c(0,3000), ylim=c(0,1))
#     lines(idx,getF(bst.h)$y[idx], col=col3, lty=1, lwd=1)
#     points(idx,getF(bst.h)$y[idx], col=col1, cex=.5, pch=16)
#     lines(idx,Z.cdf(idx, bst.pars[1],bst.pars[2])$cdf, col=col2, lwd=3)
#     legend("topleft", bty='n', legend=c("Observed","Fitted"), pch=c(16,NA), lty=c(NA,1), lwd=c(NA,3), col=c(col3,col2), seg.len = .85)
#     mtext("CDF", side = 2, line=2)
#     mtext("Burst length (ms)", side=1, line=2)
# dev.off()
#
# fn = paste0(plotlib, "turtle_xfit_rst.pdf")
# pdf(fn, height=3.5, width=4)
#     par(mar=c(2,2,2,2), oma=c(1,1,1,1), bty='n', mfrow=c(1,1))
#     plot(0,0, type='n', xlim=c(0,4000), ylim=c(0,1))
#     lines(idx,getF(rst.h)$y[idx], col=col3, lty=1, lwd=1)
#     points(idx,getF(rst.h)$y[idx], col=col1, cex=.5, pch=16, xlim=c(1000,4000))
#     lines(idx,Z.cdf(idx, rst.pars[1],rst.pars[2])$cdf, col=col2, lwd=3)
#     legend("topleft", bty='n', legend=c("Observed","Fitted"), pch=c(16,NA), lty=c(NA,1), lwd=c(NA,3), col=c(col3,col2), seg.len = .95)
#     mtext("CDF", side = 2, line=2)
#     mtext("Rest length (ms)", side=1, line=2)
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

isi250 = isi[isi < 250]
isi50  = isi[isi < 50] # 50 is about the median of the ISIs < 250: median(isi250)

hst = hist(isi, breaks = seq(0,max(isi),1),plot = FALSE)

# hst$mids[which(hst$counts== max(hst$counts))]



layout(1:2)
kts = c(quantile(isi50,probs=c(.1,.25,.5,.75)),50)     # overall knots by probs
hist(isi50, breaks=50)
kts = c(0,0,0,seq(0,50,length=4))
# kts = c(0,0,0,0,10,15,20,20,20,20)
Sb        = splines::splineDesign(kts, x=1:max(kts), ord=4, outer.ok=TRUE)
# matplot(Sb, lty=1, col=1, type='l', xlim=c(0,50))

# kts = c(rep(0,4),2^(1:10)[-(1:3)])
# Sb        = splines::splineDesign(kts, x=1:1000, ord=4, outer.ok=TRUE)

Bb0 = 0; Br0 = 0; verbose = FALSE
# kts       = c(0,0,0,seq(0,50,length=6))
# kts       = c(0,0,0,0,10,20,20,20) # very short kernel captures only refractoryness: no stochastic adaptation

# kts = c(0,0,0,seq(0,100,length=5)) # longer kernel with relatively few parameters

# Sb        = splines::splineDesign(kts, x=1:max(kts), ord=4, outer.ok=TRUE)
Bb        = rep(0,ncol(Sb))
Bb[1]     = -10 # First weight is (large) negative to ensure refractoriness
Br        = Bb[1]
# Sr        = splines::splineDesign(c(rep(0,4),max(kts)), x=1:max(kts), ord=4, outer.ok=TRUE)
maxk      = max(kts)
# Sr        = 1-log_logistic(1:(3*maxb), maxb, 10)
# Sr        = 1-log_logistic(1:(3*maxb), .5*maxb, 10)
Sr        = 1-log_logistic(1:(3*maxb), 2*maxb, 10)
M         = 200
pars      = SSGLM::make_mpfpp_par_input(Bb0, Bb, Sb, Br0, Br=Br, Sr, xpars, M=M, q0 = 1e-8, q1 = 0, s_lag = 2000, dt=1)
pars$q    = c(1e-5,rep(1e-4,ncol(Sb)), 1e-5, rep(1e-4,ncol(Sr)))

f <- function(neuron, trial){
  # neuron = 1; trial = 1;

  spk = trials[[neuron]][[trial]]
  if(length(spk)>10){
        y         = as.numeric(SSGLM::make_y(spk,1,N = 40e3, verbose = FALSE)$y)
        X         = sample(c(-1,1),M,replace=TRUE) # Sample random initial states +/- 1 for X
        # Run algorithm
        # y = y[1:5000]s
        res = SSGLM::mpfpp(y,X,pars, usetrue = FALSE, verbose = verbose)
        res$y = y
        res$est = get_est(res)
  } else{
    res = NA
  }
  return(res)
}

totres = list()

# res = f(3,1)

i=j=1

h = function(i){
  tmpres = lapply(1:10, function(j) f(neuron=i,trial=j))
  tmpres$neuron = names(trials)[i]

  # fn = paste0(datalib,"turtle_X",i,".Rda")
  fn = paste0(datalib,"turtle_",names(trials)[i],"_x2.Rda")
  save(tmpres, file = fn)
  # return(tmpres)
}

# tmpres = h(i)
# tmpres
# length(trials)
mclapply(1:length(trials), h, mc.cores = 6)

# Manual parallelization...:
# for(i in 1:30){
# for(i in 31:60){
# for(i in 61:90){
# for(i in 91:length(trials)){
#   # tmpres = list()
#   cat("Analyzing neuron ", i," (",names(trials)[i],")\n", sep="")
#   tmpres = lapply(1:10, function(j) f(neuron=i,trial=j))
#   # tmpres = mclapply(1:10, function(j) f(neuron=i,trial=j), mc.cores = 7)
#   totres[[i]] = tmpres
#   fn = paste0(datalib,"turtle_",names(trials)[i],".Rda")
#   save(tmpres, file = fn)
# }

log.logistic
