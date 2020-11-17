misc::clean_up()
library(SSGLM)
datalib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/plots/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/TeX/SSGLM plos-latex-template/gfx/"
set.seed(1234)

datalib = "/Users/jacob/Data/" #"Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"

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

isi = unlist(lapply(trials, function(x) lapply(x, diff)))


fns = list.files(datalib)
fns = fns[!grepl("turtle2",fns)]
fns = fns[!grepl("_2",fns)]
# fns = fns[grepl("_2",fns)]
fns
# pct = round(100*length(fns)/249,2)
# msg = paste0(pct,"% neurons analyzed at ",tmp)
# cat("\n",msg)

# fns = list.files(datalib)
# fns = fns[grepl("turtle2_n",fns) & grepl("qts3",fns)]
length(fns)
# idx = as.numeric(substr(fns,10,nchar(fns)-9))
# idx

tmp = substr(fns,9,nchar(fns)-4)
if(any(grepl("_",tmp))){
  tmp = substr(tmp, 1, nchar(tmp)-2)
}
idx = as.numeric(tmp)
fns = fns[order(idx)]

# neuron=idx[1]
plot_neuron <- function(neuron){
#neuron=22
  # fn = paste0(datalib,"turtle_n",neuron,".Rda")
  # fn = paste0(datalib,"turtle2_n",neuron,"_qts3.Rda")
  # load(fn)
  fn = paste0(datalib,fns[neuron])
  load(fn)
  totres = tmpres



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
    if(!is.null(res) & !all(is.na(res))){
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
  # i=1

  Btime = numeric(0)
  for(i in 1:10){
    res = totres[[i]]

    if(!all(is.na(res))){


      tmp = res$est$P
      Pest = rbind(Pest,tmp)
      np = length(tmp)
      kb[i,] = exp(res$Sb%*%as.numeric(res$est$P[2:(np-2)])+res$est$P$b0)
      kr[i,] = exp(res$Sr%*%as.numeric(res$est$P[np])+res$est$P$r0)
      Btime = c(Btime,sum(res$X>0))
    }
  }

cat("\nMean bursting kernel area", mean(kb,na.rm = TRUE))
cat("\nTotal bursting time", sum(Btime), "ms\n\n")



  # par(mfrow=c(2,1))

  plot(0,0,  type='n', xlim=c(0,nrow(res$Sb)), ylim=c(0,.05))
  for(i in 1:10){
    lines(kb[i,])
  }
  tmp = apply(Pest,2,mean)
  lines(exp(res$Sb%*%as.numeric(tmp[2:(np-2)])+tmp[1]),lwd=2, col='red')

  plot(0,0,  type='n', xlim=c(0,nrow(res$Sr)), ylim=c(0,0.003))
  for(i in 1:10){
    lines(kr[i,])
  }
  lines(exp(res$Sr%*%as.numeric(tmp[np])+tmp[np-1]),lwd=2, col='red')

  # kts       = c(20,40,60,80,100,120,250)    # add auxiliary knot at max burst 250 ms

  # kts = c(21,36,66,115, 250)
  # paste0("lag",c(0,kts,0,max(kts)))
  # (1-log_logistic(1:(5*250), 250, 10))[250]

  # layout(1)
  # image(round(res$est$cor[1:7,1:7],3))

  est = round(rbind(apply(Pest,2,mean), 2*apply(Pest,2,sd)),3)

  # colnames(est)= paste0("lag",c(0,kts,0,max(kts)))
  return(est)
}

# plot_neuron(idx[4])

idxn = c(22, 156, 177, 179)
plot_neuron(22)
plot_neuron(156)

plot_neuron(177)
plot_neuron(179)

idxn2 = c(95, 135)#,14)
idxn3 = c(12,13)

plot_neuron(135)


id = c(171,103,163)
plot_neuron(171)
plot_neuron(103)
plot_neuron(163)

# for(i in 1:length(idx)){
#   cat("\n-------------------------")
#   cat("\nNumer", i, " Neuron",idx[i],"\n")
#   try(plot_neuron(idx[i]))
#   Sys.sleep(2)
# }


i=j=1
# Pvals = matrix(nr=length(idx),nc=10)
# for(i in 1:length(idx)){
#   fn = paste0(datalib,"turtle_n",idx[i],".Rda")
#   load(fn)
#   for(j in 1:10){
#     tmp = tmpres[[j]]
#     if(!all(is.na(tmp))){
#       Pvals[i,j] = tmp$est$ks$p.value
#     }
#   }
# }

# Pvals2 = matrix(nr=length(idx),nc=10)
# for(i in 1:length(idx)){
#   fn = paste0(datalib,"turtle_n",idx[i],"_qts.Rda")
#   load(fn)
#   for(j in 1:10){
#     tmp = tmpres[[j]]
#     if(!all(is.na(tmp))){
#       Pvals2[i,j] = tmp$est$ks$p.value
#     }
#   }
# }

na.mean = function(x) {
  mean(x, na.rm=TRUE)
}
na.median = function(x) {
  median(x, na.rm=TRUE)
}

Pvals = matrix(nr=length(idx),nc=10)
Ps = SEs = array(dim=c(length(idx),8,10), dimnames = list(paste("n",1:length(idx)),c(paste0("b",0:5),paste0("r",0:1)),paste("trial",1:10)))

for(i in 1:length(idx)){
  cat("\nLoading neuron",i)
  fn = paste0(datalib,"turtle2_n",idx[i],"_qts3.Rda")
  load(fn)
  for(j in 1:10){
    tmp = tmpres[[j]]
    if(!all(is.na(tmp))){
      Pvals[i,j] = tmp$est$ks$p.value
      tmp = rbind(tmp$est$P,sqrt(diag(tmp$est$cov)))
      Ps[i,,j] = as.numeric(tmp[1,])
      SEs[i,,j] = as.numeric(tmp[2,])
    }
  }
  cat(" - loaded!")
}

tmp1 = cbind(1:10,apply(round(100*Pvals,2)>1,2, na.mean),apply(round(100*Pvals,2)>5,2, na.mean))
tmp2 = cbind(1:length(idx),apply(round(100*Pvals,2)>1,1, na.mean),apply(round(100*Pvals,2)>5,1, na.mean))

layout(1)
plot(tmp2[,1:2]);abline(h=c(.01,.05), lty=3)
plot(tmp2[,c(1,3)]);abline(h=c(.01,.05), lty=3)
plot(tmp2[,2:3]);abline(0,1);abline(h=c(.01,.05), v=c(.01,.05), lty=3)
mean(tmp2[,3]>0.5)

round(100*Pvals,2)

boxplot(Pvals);abline(h=c(.01,.05), lty=3)

hist(Pvals[,1], breaks = seq(0,1,.01));abline(v=.05)

ord = order(apply(Pvals,1, na.median), decreasing = TRUE)
boxplot(t(Pvals[ord,]));abline(h=c(.01,.05), lty=3)
2*117/3
round(tmp1,2)
round(tmp2,2)


# kts = c(21,36,66,115, 250)
# paste0("lag",c(0,kts,0,max(kts)))


round(tmp1,2)
round(tmp2,2)

tmpvals = Pvals
tmpvals[tmpvals > 0.01] = NA
image(tmpvals[ord,], col=heat.colors(100,.8), axes=FALSE, ylab="Trials", xlab="Neurons")
axis(2, seq(0,1,length=10), 1:10, tick = FALSE, las=1)


plot(cbind(1:10,apply(spks_in_trials,2,mean)))
plot(cbind(1:10,apply(spks_in_trials,2,median)))
plot(cbind(1:10,apply(spks_in_trials,2,sd)))

hist(spks_in_trials[,1], breaks=100)
hist(spks_in_trials[,8], breaks=100)


tmpvals = Pvals[,-c(3,8,9)]
boxplot(tmpvals);abline(h=c(.01,.05), lty=3)
tmpvals[tmpvals > 0.01] = NA
image(tmpvals[ord,], col=heat.colors(100,.8), axes=FALSE, ylab="Trials", xlab="Neurons")
axis(2, seq(0,1,length=ncol(tmpvals)), 1:ncol(tmpvals), tick = FALSE, las=1)

par(mfrow=c(2,4))
for(i in 1:8){
  boxplot((Ps[,i,]))
}

tmp = apply(spks_in_trials[,-9], 1, mean)
# quantile(tmp)

idx = which(tmp<100 & tmp>25) # roughly 25% and 75% quantiles
trials = trials[idx]

spks_in_trials = spks_in_trials[idx,]

layout(1)
boxplot(spks_in_trials)

# tmp = which(spks_in_trials>500,arr.ind = TRUE)
# tmp[order(tmp[,1]),]
# Pvals[tmp[,1],]


Ps[1,,]
plot(exp(apply(Ps[,7,],2,na.mean)))
# SEs[1,,]
# sum(round(100*Pvals)>1, na.rm = TRUE)
# sum(round(100*Pvals2)>1, na.rm = TRUE)
mean(round(100*Pvals)>1, na.rm = TRUE)

misc::clear_console()

# 1*(round(100*Pvals)>=1)
# 1*(round(100*Pvals2)>=1)
1*(round(100*Pvals3)>=1)





Ps[1,,]
round(apply(Ps,2, na.mean),2)
round(2*apply(Ps,2:3, function(x) sd(x, na.rm =TRUE)),2)


Ptot = apply(Ps,2, na.mean)
Sb   = tmpres[[1]]$Sb
Sr   = tmpres[[1]]$Sr

par(mfrow=c(2,2))
plot(0,0, type='n', bty='n', xlim=c(0,250), ylim=c(0,0.05))
for(i in 1:length(idx)){
  Ptmp = apply(Ps[i,,],1,na.mean)
  ktmp = Sb%*%Ptmp[2:6]
  ktmp = ktmp+Ptmp[1]
  lines(exp(ktmp), col=add.alpha('black',.25))
}
ktmp = Sb%*%Ptot[2:6]
ktmp = ktmp+Ptot[1]
lines(exp(ktmp), col=add.alpha('red',.5), lwd=2)
abline(v=c(21,36,66,115, 250), lty=3)

plot(0,0, type='n', bty='n', xlim=c(0,250), ylim=c(0,1.5))
for(i in 1:length(idx)){
  Ptmp = apply(Ps[i,,],1,na.mean)
  ktmp = Sb%*%Ptmp[2:6]
  # ktmp = ktmp+Ptmp[1]
  lines(exp(ktmp), col=add.alpha('black',.25))
}
ktmp = Sb%*%Ptot[2:6]
# ktmp = ktmp+Ptot[1]
lines(exp(ktmp), col=add.alpha('red',.5), lwd=2)
abline(v=c(21,36,66,115, 250), lty=3)

plot(0,0, type='n', bty='n', xlim=c(0,1250), ylim=c(0,0.005))
for(i in 1:length(idx)){
  Ptmp = apply(Ps[i,,],1,na.mean)
  ktmp = Sr%*%Ptmp[8]
  ktmp = ktmp+Ptmp[7]
  lines(exp(ktmp), col=add.alpha('black',.25))
}
ktmp = Sr%*%Ptot[8]
ktmp = ktmp+Ptot[7]
lines(exp(ktmp), col=add.alpha('red',.5), lwd=2)

plot(0,0, type='n', bty='n', xlim=c(0,1250), ylim=c(0,1))
for(i in 1:length(idx)){
  Ptmp = apply(Ps[i,,],1,na.mean)
  ktmp = Sr%*%Ptmp[8]
  # ktmp = ktmp+Ptmp[7]
  lines(exp(ktmp), col=add.alpha('black',.25))
}
ktmp = Sr%*%Ptot[8]
# ktmp = ktmp+Ptot[7]
lines(exp(ktmp), col=add.alpha('red',.5), lwd=2)




