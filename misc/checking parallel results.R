library(misc)
library(SSGLM)
library(parallel)
clean_up()
datalib = "/Users/jacob/Data/" #"Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"
set.seed(1234)


na.mean = function(x) {
  mean(x, na.rm=TRUE)
}
na.median = function(x) {
  median(x, na.rm=TRUE)
}
na.sum = function(x) {
  sum(x, na.rm=TRUE)
}
na.min = function(x){
  min(x, na.rm=TRUE)
}

dat = R.matlab::readMat(misc::icloud_lib("GitHub/Source/R/Extracellular triggered EPSP/ExtracellularUnits.mat")) #okay
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


# format(Sys.time(), "%a %b %d %H:%M:%S %Y")
  tmp = Sys.time()
# List current files in datalib
  fns = list.files(datalib)
  fns = fns[!grepl("turtle2",fns)]
  # fns = fns[!grepl("x",fns)]
  fns = fns[grepl("x",fns)]
  fns
  pct = round(100*length(fns)/249,2)
  msg = paste0(pct,"% neurons analyzed at ",tmp)
  cat("\n",msg)

  # 5.22%   neurons analyzed at 2020-10-14 11:00:34
  # 13.25%  neurons analyzed at 2020-10-14 12:00:36
  # 24.5%   neurons analyzed at 2020-10-14 13:00:37
  # 33.73%  neurons analyzed at 2020-10-14 14:00:39
  # 42.17%  neurons analyzed at 2020-10-14 15:00:40
  # 51%     neurons analyzed at 2020-10-14 16:00:41
  # 60.64%  neurons analyzed at 2020-10-14 17:00:43
  # 69.48%  neurons analyzed at 2020-10-14 18:00:44
  # 78.71%  neurons analyzed at 2020-10-14 19:00:45
  # 87.55%  neurons analyzed at 2020-10-14 20:00:46
  # 96.79%  neurons analyzed at 2020-10-14 21:00:48
  # 100%    neurons analyzed at 2020-10-14 22:00:49%

  tmp = substr(fns,9,nchar(fns)-4)
  if(any(grepl("_",tmp))){
    tmp = substr(tmp, 1, nchar(tmp)-2)
  }
  idx = as.numeric(tmp)
  fns = fns[order(idx)]

  g <- function(x){
    if(!all(is.na(x))){
      p = x$est$ks$p.value
      P = rbind(x$est$P,sqrt(diag(x$est$cov)))
    } else{
      p = NA
      P = NA
    }
    return(list(pval=p, pars=P))
  }

  f <- function(i){
    fn = paste0(datalib,fns[i])
    load(fn)
    res = list()
    for(j in 1:10){
      res[[j]] = g(tmpres[[j]])
    }
    return(res)
  }


  # res = lapply(1:2, f)
  res = mclapply(1:length(fns), f, mc.cores = 1)
  # res = mclapply(1:length(fns), f, mc.cores = 6)

  load(paste0(datalib,fns[1]))

  tmp = tmpres[[1]]
  rm("tmpres")
  Sb = tmp$Sb
  Sr = tmp$Sr
  rm("tmp")

  nPb = ncol(Sb)+1
  nPr = ncol(Sr)+1
  nP  = nPb+nPr
  Pvals = matrix(NA,nr=length(idx),nc=10)
  Ps = SEs = array(NA,dim=c(length(idx),nP,10), dimnames = list(paste("n",1:length(idx)),c(paste0("b",0:(nPb-1)),paste0("r",0:(nPr-1))),paste("trial",1:10)))

  i=j=1
  for(i in 1:length(res)){
    tmp = res[[i]]
    for(j in 1:10){
      x = tmp[[j]]
      Pvals[i,j] = x$pval
      if(!all(is.na(x$pars))){
        Ps[i,,j]   = as.numeric(x$pars[1,])
        SEs[i,,j]  = as.numeric(x$pars[2,])
      }
    }
  }

  # for(i in 1:length(fns)){
  #   cat("\nLoading neuron",i)
  #   fn = paste0(datalib,fns[i])
  #   load(fn)
  #   for(j in 1:10){
  #     tmp = tmpres[[j]]
  #     if(!all(is.na(tmp))){
  #       Pvals[i,j] = tmp$est$ks$p.value
  #       tmp = rbind(tmp$est$P,sqrt(diag(tmp$est$cov)))
  #       Ps[i,,j] = as.numeric(tmp[1,])
  #       SEs[i,,j] = as.numeric(tmp[2,])
  #     }
  #   }
  #   cat(" - loaded!")
  # }

  idx[sort(unique(which(is.na(Pvals), arr.ind = TRUE)[,1]))]


  spks_in_trials
  tmp_spk = spks_in_trials[idx,]



  boxplot(tmp_spk)



  prod(dim(Pvals))
  sum(is.na(Pvals))
  mean(is.na(Pvals))
  tmp_spk[is.na(Pvals)]

  x = seq(1e-8,.1,length=100)
  y = numeric(length(x))
  for(i in 1:length(x)){
    y[i] = na.mean(Pvals>x[i])
  }

  plot(x,y, type='l', log = "x")

  boxplot((Pvals)); abline(h=(c(.01,.05)), lty=3)
  ord = order(apply(Pvals,1,na.median), decreasing = TRUE)
  boxplot(t(Pvals[ord,])); abline(h=(c(.01,.05)), lty=3)

  plot(apply(Pvals>.05, 1, na.mean))

  na.mean(Pvals>.05)
  na.mean(Pvals>.01)
  image(Pvals>.05, col = add.alpha(c('red','dodgerblue'),.75))

  # apply(Pvals>.05)
  par(mfrow=c(2,1), mar=c(3,3,1,1))
  image(Pvals>.05, col = add.alpha(c('red','dodgerblue'),.75))
  image(Pvals>.01, col = add.alpha(c('red','dodgerblue'),.75))

  # plot(spks_in_trials[Pvals<0.05])
  cols = add.alpha(c('red','dodgerblue'),.75)

  layout(1)
  image(spks_in_trials, col=ifelse(Pvals>.05,cols[1],cols[2]))

  na.mean(apply(Pvals>.01, 1, na.mean)>.66)
  na.mean(apply(Pvals>.05, 1, na.mean)>.66)

  cat("Pct of trials with p-val >.01 is",na.mean(Pvals>.01))
  cat("Pct of trials with p-val >.05 is",na.mean(Pvals>.05))

  na.mean(apply(Pvals[,-c(3,8,9)]>.01, 1, na.mean)==1)
  na.mean(apply(Pvals[,-c(3,8,9)]>.05, 1, na.mean)==1)


  ord = order(apply(Pvals>.05, 1, na.mean))

  Pvals.05 = 1*(Pvals>.05)
  Pvals.01 = 1*(Pvals>.01)

  x = 1:10
  y = apply(Pvals, 2, function(x) mean(is.na(x)))
  plot(x,y, type='n')
  text(x,y, labels = as.character(x))


  ord = order(apply(spks_in_trials, 1, median), decreasing=TRUE)
  boxplot(t(spks_in_trials[ord,]))


  spks2 = spks_in_trials[idx,]

  x = as.vector(spks2)
  y = as.vector(Pvals)
  # x = apply(spks_in_trials[idx,], 1, median)
  # y = apply(Pvals>.05,1, na.mean)

  # hist(spks2, breaks=1000, xlim=c(0,100))
  plot(x,y, pch=16, col=add.alpha('red',.75), xlim=c(0,200))
  abline(h=.05, v=10)

  boxplot(t(Pvals[ord,c(3,8,9)]))

  round(100*Pvals)

  isi = lapply(trials, function(x) lapply(x, diff))

  ISI = matrix(nr=length(res), nc=10)
  for(i in 1:length(res)){
    tmpisi = isi[[i]]
    ISI[i,] = unlist(lapply(tmpisi, median))
  }

  x = as.vector(ISI[idx,])
  y = as.vector(Pvals)
  plot(x,jitter(y,amount = .01), pch=16, col=add.alpha('red',.75), xlim=c(0,.2))
  abline(h=c(.01,.05))

  # load(paste0(datalib,fns[1]))
  #
  # Sb = tmpres[[1]]$Sb
  # Sr = tmpres[[1]]$Sr

  # Individual neurons
  Pb = apply(Ps,1:2, na.mean)[,1:nPb]
  Pr = apply(Ps,1:2, na.mean)[,(nPb+1):nP]

  # Individual trials
  dim(Ps)
  Pb = t(apply(Ps,2:3, na.mean))[,1:nPb]
  Pr = t(apply(Ps,2:3, na.mean))[,(nPb+1):nP]

  # Kernels
  Kb = Sb%*%t(Pb[,-1])
  Kr = Sr%*%t(Pr[,-1])
  # Baselines
  Bb = matrix(rep(Pb[,1],nrow(Kb)),nr=nrow(Kb), byrow=TRUE)
  Br = matrix(rep(Pr[,1],nrow(Kr)),nr=nrow(Kr), byrow=TRUE)

  matplot(exp(Kb+Bb), type='l', lty=1, col=add.alpha('black',.25))
  matplot(exp(Kr+Br), type='l', lty=1, col=add.alpha('black',.25))


  Kb = Kb2 = matrix(NA, nr=10*length(res), nc=nrow(Sb))
  Kr = Kr2 = matrix(NA, nr=10*length(res), nc=nrow(Sr))
  Br = Bb = numeric(10*length(res))
  for(i in 1:length(res)){
    for(j in 1:10){
      # i=j=1
      Kb2[10*(i-1)+j,] = Sb%*%Ps[i,2:nPb,j]
      Kr2[10*(i-1)+j,] = Sr%*%Ps[i,(nPb+2):nP,j]

      Br[10*(i-1)+j] = Ps[i,1,j]
      Br[10*(i-1)+j] = Ps[i,nPb+1,j]

      Kb[10*(i-1)+j,] = Sb%*%Ps[i,2:nPb,j] + Ps[i,1,j]
      Kr[10*(i-1)+j,] = Sr%*%Ps[i,(nPb+2):nP,j]   + Ps[i,nPb+1,j]

    }
  }
  Kb = t(Kb)
  Kb2 = t(Kb2)
  Kr = t(Kr)
  Kr2 = t(Kr2)

  layout(1:2)
  matplot(exp(Kb), type='l', lty=1, col=add.alpha('black',.06), ylim=c(0,.02))
  lines(exp(apply(Kb,1,na.mean)), col=add.alpha('red',.75),lwd=2)
  lines(exp(apply(Kb,1,na.median)), col=add.alpha('orange',.75),lwd=2)
  matplot(exp(Kr), type='l', lty=1, col=add.alpha('black',.1), ylim=c(0,.002))
  lines(exp(apply(Kr,1,na.mean)), col=add.alpha('red',.75),lwd=2)
  lines(exp(apply(Kr,1,na.median)), col=add.alpha('orange',.75),lwd=2)




# Find neurons that are close to population average:
  tmp = Kb-matrix(apply(Kb,1,na.mean),nr=nrow(Kb),nc=ncol(Kb))
  tmp = apply(tmp, 2, function(x) na.sum(x^2))

  length(unique(which(is.na(Kb), arr.ind=TRUE)[,2]))
  idxn = which(tmp>0 & tmp < 10) # these have low MSE, MSE=0 means NA kernels...

  tmp = apply(matrix(tmp,nr=10),2,sum) # find total MSE over all trials
  idxn = which(tmp>0 & tmp < 60)

  idxn = idxn[idxn %in% which(apply(is.na(Pvals),1,sum)<3)] # include only neurons with less than 3 missed trials (NA values)
  idxn

  # make index vector of all neurons trials among the 2490 cols:
  idx.tmp = matrix(10*(idxn-1),nr=length(idxn),nc=10)+matrix(1:10,nr=length(idxn),nc=10, byrow = TRUE)
  idx.tmp = sort(as.vector(idx.tmp))

  # Plot the estimated kernels for all the trials of the idxn neurons:
  layout(1)
  matplot(exp(Kb[,idx.tmp]), type='l', lty=1, col=add.alpha('black',.35), ylim=c(0,.02))

  idxn
  which(apply(is.na(Kb[,idx.tmp]),2,all))

  # check p-values of the idxn neurons
  1*round(100*Pvals[idxn,])>5
i=1
layout(1)
plot(0,0, type='n', xlim=c(0,400), ylim=c(0,length(idxn)), bty='n', xlab='ms', axes=FALSE, ylab="neurons")
for(i in 1:length(idxn)){
  for(j in 1:10){
    tmp = trials[[idxn[i]]][[j]]
    points(tmp+(j-1)*40, rep(i, length(tmp)), pch=124)
  }
}
abline(v=(0:10)*40, lty=3)
text((0:9)*40+20, rep(251,10), labels = as.character(1:10), pos = 3)

par(mfrow=c(2,2))
for(i in 1:length(idxn)){
  tmp = 1e3*unlist(isi[idxn[i]])
  hist(tmp, xlim=c(0,250), prob=TRUE, breaks=c(seq(0,260,5),1e6))
}





  layout(1)
  matplot(exp((Kb[,unique(which(exp(Kb2)>1, arr.ind = TRUE)[,2])])), type='l', col=add.alpha('black',.1), lty=1)






  layout(1)
  plot(exp(apply(Kb,1,na.mean)), type='l')
  plot(exp(apply(Kb2,1,na.mean)), type='l')


  isi.all = unlist(isi)
  idxn = which(apply(Pvals>.05,1,na.sum)>0) # neuron index
  layout(1)
  hst = hist(1e3*isi.all, breaks=c(seq(0,150,.1),1e6), xlim=c(0,100), ylim = c(0,250), freq = TRUE, border=NA, col=add.alpha('red',.75))
  x = hst$mids*1e3
  y = hst$density
  y = y[x<100]
  x = x[x<100]
  par(mfrow=c(2,5))
  for(j in 1:10){
    idxt = seq(1,10*249,10)+j-1 # trial index
    matplot(exp(Kb[,idxt[idxn]]), type='l', lty=1, col=add.alpha('black',.06), ylim=c(0,.01), xlim=c(0,50))
    lines(exp(apply(Kb[,idxt[idxn]],1,na.mean)), col=add.alpha('red',.75),lwd=2)
    lines(exp(apply(Kb[,idxt[idxn]],1,na.median)), col=add.alpha('orange',.75),lwd=2)
    mtext(paste("Trial",j), 3, -1.5)
    par(new=TRUE)
    plot(x,y,type='l', xlim=c(0,50), col=add.alpha('dodgerblue',.75), lwd=2, yaxt='n')
  }





  layout(1)



  idxn = 1:length(res)
  spks = trials[idxn]
  Hst = numeric(0)
  tmp.isi = lapply(spks, function(x) 1e3*unlist(lapply(x,diff)))
  for(i in 1:length(res)){
    tmp.hst = hist(tmp.isi[[i]], breaks=c(seq(0,100,.5),1e6), plot=FALSE)
    x = tmp.hst$mids
    n = length(x)
    x = x[-n]
    y = tmp.hst$density[-n]
    Hst = rbind(Hst,y)
  }

  layout(1)
  # matplot(x,t(Hst), type='l', lty=1, col=add.alpha('black',.15), ylim=c(0,.02), xlim=c(0,100))


  boxplot(Hst, xaxt='n', ylim=c(0,.02), xlim=c(0,ncol(Hst)))
  # axis(1, at=1:20, seq(5,100,5))
  Kb = Sb%*%t(Pb[,-1])
  Kb = Kb+matrix(Pb[,1], nr=nrow(Kb), nc=10, byrow=TRUE)
  par(new=TRUE)
  matplot(exp(Kb), type='l', lty=1, col=add.alpha('red',.5), lwd=2, ylim=c(0,.02), xlim=c(0,100))
  # matplot(Kb, type='l', lty=1, col=add.alpha('red',.5), lwd=2, xlim=c(0,100))


  ord = order(apply(Hst,1, median), decreasing=TRUE)
  boxplot(t(Hst[ord,]), xaxt='n')

  Hst[1,]

  par(mfrow=c(4,5))
  idxn = 1:20
  spks = trials[idxn]
  tmp.isi = lapply(spks, function(x) 1e3*unlist(lapply(x,diff)))
  for(i in 1:length(idxn)){
    hist(tmp.isi[[i]], breaks=c(seq(0,100,1),1e6), xlim=c(0,100), prob=TRUE, ylim=c(0,.025))
  }



  # plot(spks,rep(1,length(spks)), pch=124)
  par(mfrow=c(2,5))
  for(j in 1:10){
    tmp.spk = spks[[j]]
    if(is.empty(tmp.spk)){
      plot(0,0, type='n', bty='n', axes=FALSE,ann=FALSE)
    } else{

    }
  }

  # tmp.spk = spks[[7]]
  # plot(1e3*tmp.spk, rep(1, length(tmp.spk)), pch=124)
  # diff(1e3*spks[[7]])



  isi = lapply(trials, function(x) lapply(x, diff))

  idxn1 = which(unlist(lapply(isi, function(y) any(unlist(lapply(y, function(x) x <= 1e-3 )))))) # neurons with ISI's <= 1ms
  idxn2 = setdiff(1:length(trials), idxn1) # neurons with min(ISI's)>1 ms


  # idxn1 = which(unlist(lapply(isi, function(y) any(unlist(lapply(y, function(x) x <= 2e-3 )))))) # neurons with ISI's <= 2ms
  # idxn2 = setdiff(1:length(trials), idxn1) # neurons with min(ISI's)>2 ms

  length(idxn1)/249

  isi.all = 1e3*unlist(isi[idxn2])
  hist(isi.all, xlim=c(0,100), breaks=c(seq(0,150,1),1e6))


  na.mean(Pvals>.01)
  na.mean(Pvals>.05)

  na.mean(Pvals[idxn1,]>.01)
  na.mean(Pvals[idxn2,]>.01)

  Pvals[idxn2,]
  unique(which(is.na(Pvals),arr.ind = TRUE)[,1])
  round(100*Pvals[idxn2,])

  idxnt = which(Pvals[idxn2,]>.01, arr.ind = TRUE) # neurons and trials eligible for parameter aggregation

  prod(dim(Pvals[idxn2,]))
  dim(Ps)

  matplot(Sb, type='l', col='black')
  idxnt

  tmpP = numeric(nP)
  for(i in 1:nP){
    tmpP[i] = na.mean(Ps[idxnt[1],i,idxnt[,2]])
  }

  layout(1:2)
  plot(exp(Sb%*%tmpP[2:nPb]+tmpP[1]), type='l')
  matplot(Sb, type='l', col='black')
tmpP
