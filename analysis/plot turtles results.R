library(misc)
library(SSGLM)
library(parallel)
clean_up()
datalib = "/Users/jacob/Data/" #"Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/TeX/SSGLM plos-latex-template/gfx/"

set.seed(1234)
savePDF = FALSE

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

# List files in datalib
fns = list.files(datalib)
fns = fns[grepl("turtle_",fns)]
fns = fns[!grepl("x",fns)]
# fns = fns[grepl("x",fns)]
# fns = fns[grepl("x2",fns)]
tmp = substr(fns,9,nchar(fns)-4)
if(any(grepl("_",tmp))){
  if(any(grepl("x2",tmp))){
    tmp = substr(tmp, 1, nchar(tmp)-3)
  }else{
    tmp = substr(tmp, 1, nchar(tmp)-2)
  }

}

idx = as.numeric(tmp)
fns = fns[order(idx)]


dat = R.matlab::readMat(misc::icloud_lib("GitHub/Source/R/Extracellular triggered EPSP/ExtracellularUnits.mat")) #okay
turtles = list()
for(i in 1:249)
  eval(parse( text=paste0("turtles$n",i," <- dat$units[[1]][[i]][[1]]") ))

ntrains = length(turtles)

# Convert spike times to seconds
tosec   = 2.5e-5
for(i in 1:ntrains)
  turtles[[i]] = turtles[[i]]*tosec

trials = list()
for(i in 1:ntrains){
  trials[[i]] = get_spiketimes(i)
}
names(trials) = paste0("n",1:ntrains)

turtles = turtles[ idx[order(idx)] ]
trials = trials[ idx[order(idx)] ]

spks_in_trials = matrix(unlist(lapply(trials,function(x) unlist(lapply(x,length)))), nr=length(trials),nc=10, byrow=TRUE)

# Load results on p-values and parameter estimates from all 249 spiketrains
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
  res = mclapply(1:length(fns), f, mc.cores = 6)

# Load spline matrices for bursting/resting
  load(paste0(datalib,fns[1]))
  Sb = tmpres[[1]]$Sb
  Sr = tmpres[[1]]$Sr
  rm("tmpres")

# Number of parameters for bursting/resting
  nPb = ncol(Sb)+1
  nPr = ncol(Sr)+1
  nP  = nPb+nPr


# Set parameter estimates into arrays
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
  mean(apply(Pvals>0.05,2,na.mean))
  na.mean(Pvals>.05)
  na.mean(Pvals>.01)

  isi = 1e3*unlist( lapply(trials, function(x) lapply(x,function(y) diff(sort(y)))) )

if(savePDF) pdf(paste0(plotlib,"249neurons_isi.pdf"), width=10, height=4)
        par(mfrow=c(1,2), mar=c(3,3,1,1), oma=c(0,0,0,0))
        hist(isi, breaks=c(seq(0,1e3+1,.5),1e6), xlim=c(0,1e3), border=NA, col=add.alpha('black',.25), main="", xlab="", ylab="", ylim=c(0,.015))
        mtext("A",side=3, line=0, adj=0, family = 'sans', cex=1.25)
        mtext("ms",1, line=2)
        mtext("PDF",2, line=2.5)
        hist(isi, breaks=c(seq(0,1e3+1,.5),1e6), xlim=c(0,100), border=NA, col=add.alpha('black',.5), main="", xlab="", ylab="", ylim=c(0,.015))
        mtext("B",side=3, line=0, adj=0, family = 'sans', cex=1.25)
        mtext("ms",1, line=2)
        mtext("PDF",2, line=2.5)
if(savePDF)  dev.off()

if(savePDF) pdf(paste0(plotlib,"249neurons_spiketrains.pdf"), width=5, height=6)
        par(mfrow=c(1,1), mar=c(3,3,1,1), oma=c(0,0,0,0))
        plot(0,0, type='n', xlim=c(0,400), ylim=c(0,250), bty='n', xlab='ms', axes=FALSE, ylab="neurons")
        for(i in 1:249){
          for(j in 1:10){
            tmp = trials[[i]][[j]]
            points(tmp+(j-1)*40, rep(i, length(tmp)), pch='.')
          }
        }
        abline(v=(0:10)*40, lty=3)
        text((0:9)*40+20, rep(251,10), labels = as.character(1:10), pos = 3)
if(savePDF)  dev.off()



      # Individual neurons
        Pb.neurons = apply(Ps,1:2, na.mean)[,1:nPb]
        Pr.neurons = apply(Ps,1:2, na.mean)[,(nPb+1):nP]

      # Individual trials
        Pb.trials = t(apply(Ps,2:3, na.mean))[,1:nPb]
        Pr.trials = t(apply(Ps,2:3, na.mean))[,(nPb+1):nP]

      # Kernels
        Kb.neurons = Sb%*%t(Pb.neurons[,-1])
        Kr.neurons = Sr%*%t(Pr.neurons[,-1])
        Kb.trials = Sb%*%t(Pb.trials[,-1])
        Kr.trials = Sr%*%t(Pr.trials[,-1])

      # Baselines
        Bb.neurons = matrix(rep(Pb.neurons[,1],nrow(Sb)),nr=nrow(Sb), byrow=TRUE)
        Br.neurons = matrix(rep(Pr.neurons[,1],nrow(Sr)),nr=nrow(Sr), byrow=TRUE)
        Bb.trials = matrix(rep(Pb.trials[,1],nrow(Sb)),nr=nrow(Sb), byrow=TRUE)
        Br.trials = matrix(rep(Pr.trials[,1],nrow(Sr)),nr=nrow(Sr), byrow=TRUE)

      # All kernels
        Kb.all = Kb.all.nobase = matrix(NA, nr=10*length(res), nc=nrow(Sb))
        Kr.all = Kr.all.nobase = matrix(NA, nr=10*length(res), nc=nrow(Sr))
        Br.all = Bb.all = numeric(10*length(res))
        for(i in 1:length(res)){
          for(j in 1:10){
            Kb.all.nobase[10*(i-1)+j,] = Sb%*%Ps[i,2:nPb,j]
            Kr.all.nobase[10*(i-1)+j,] = Sr%*%Ps[i,(nPb+2):nP,j]

            Br.all[10*(i-1)+j] = Ps[i,1,j]
            Br.all[10*(i-1)+j] = Ps[i,nPb+1,j]

            Kb.all[10*(i-1)+j,] = Sb%*%Ps[i,2:nPb,j] + Ps[i,1,j]
            Kr.all[10*(i-1)+j,] = Sr%*%Ps[i,(nPb+2):nP,j]   + Ps[i,nPb+1,j]

          }
        }
        Kb.all = t(Kb.all)
        Kr.all = t(Kr.all)
        Kb.all.nobase = t(Kb.all.nobase)
        Kr.all.nobase = t(Kr.all.nobase)

if(savePDF) pdf(paste0(plotlib,"kernel_overview.pdf"),width=10,height=6)
# if(savePDF) pdf(paste0(plotlib,"kernel_overview_125.pdf"),width=10,height=12)
        par(mfrow=c(2,3), bty='n', mar=c(3,3,1,0), oma=c(0,0,2,0), cex.axis=1.2)
        matplot(exp(Kb.all), type='l', lty=1, col=add.alpha('black',.05), ylim=c(0,.05), xlab="", ylab="")
        mtext("A",side=3, line=0, adj=0, family = 'sans', cex=1.25); mtext("ms", side=1, line=2)
        matplot(exp(Kb.neurons+Bb.neurons), type='l', lty=1, col=add.alpha('black',.1), ylim=c(0,.05), xlab="", ylab="")
        mtext("B",side=3, line=0, adj=0, family = 'sans', cex=1.25); mtext("ms", side=1, line=2)
        matplot(exp(Kb.trials+Bb.trials), type='l', lty=1, col=add.alpha('black',.25), ylim=c(0,.05), xlab="", ylab="")
        mtext("C",side=3, line=0, adj=0, family = 'sans', cex=1.25); mtext("ms", side=1, line=2)
        matplot(exp(Kr.all), type='l', lty=1, col=add.alpha('black',.05), ylim=c(0,.002), xlab="", ylab="")
        mtext("D",side=3, line=0, adj=0, family = 'sans', cex=1.25); mtext("ms", side=1, line=2)
        matplot(exp(Kr.neurons+Br.neurons), type='l', lty=1, col=add.alpha('black',.1), ylim=c(0,.002), xlab="", ylab="")
        mtext("E",side=3, line=0, adj=0, family = 'sans', cex=1.25); mtext("ms", side=1, line=2)
        matplot(exp(Kr.trials+Br.trials), type='l', lty=1, col=add.alpha('black',.25), ylim=c(0,.002), xlab="", ylab="")
        mtext("F",side=3, line=0, adj=0, family = 'sans', cex=1.25); mtext("ms", side=1, line=2)
if(savePDF) dev.off()


# Select two neurons to analyze
      idx1 = apply(Pvals[,-9],1,function(x) sum(is.na(x)) ) <= 1  # should have at most 1 missing trials besides trial 9
      idx2 = apply(Pvals>.05,1,na.mean)==1                        # all non-missing p-values should be >.05
      idxn = which(idx1 & idx2)                                   # neurons eligible for analysis
      idxn = idxn[order(Pb.neurons[which(idx1 & idx2),1])[1:3]]   # pick the ones with the 3 smallest baselines

      spks_in_trials[idxn,] # only few spikes for these neurons trials

      idx3 = apply(spks_in_trials,1,median)>100                   # pick neurons with a large number of spikes
      idx4 = apply(Pvals>.01,1,na.mean)>.5                        # more than half of non-missing p-values should be >.01


      idxn = c(idxn,which(idx3 & idx4))                           # the final 6 neurons we look at
      spks_in_trials[idxn,]

      # idxn = idxn[1:2]
      idxn = c(103,111, 107)  # pick these two, 107 is picked in addition as a poorly fitting neuron

      spks_in_trials[idxn,]
      round(Pvals[idxn,],3) # looks good

      xtable::xtable(round(Pvals[idxn,],4))

      Kb.idx = Kb.neurons[,idxn]+Bb.neurons[,idxn]
      Kr.idx = Kr.neurons[,idxn]+Br.neurons[,idxn]

      # idxx = c(10*(idxn[1]-1)+1:10,10*(idxn[2]-1)+1:10)
      idxx = as.numeric(outer(10*(idxn-1),1:10,"+"))

      Kb.idx.trials = Kb.all[,idxx]
      Kr.idx.trials = Kr.all[,idxx]


if(savePDF) pdf(paste0(plotlib,"3neurons_kernels.pdf"), width=5, height=4)
      # if(savePDF) pdf(paste0(plotlib,"2neurons_kernels_125.pdf"), width=10, height=5)
      par(mfrow=c(1,1), mar=c(3,3,1,1), oma=c(0,0,0,0))
      cols = c('dodgerblue','red', 'springgreen4')
      col1 = add.alpha(cols[1],.75)
      col2 = add.alpha(cols[2],.75)
      col3 = add.alpha(cols[3],.75)
      # matplot(exp(Kb.idx.trials), type='l', lty=3, lwd=1, col=c(rep(col1,10),rep(col2,10),rep(col3,10)), ylim=c(0,.015), xlab="ms", ylab=paste0("Trials for neurons (",paste(idxn, collapse = ","),")"), bty='n')
      # par(new=TRUE)
      matplot(exp(Kb.idx), type='l', lty=1, lwd=2, col=c(col1,col2,col3), ylim=c(0,.015), xlab="", ylab="", bty='n')
      mtext("Bursting kernels", 3, line=0)
      mtext("ms", 1, line=2)
      # mtext("Trial average", 2, line=2)
      legend("topright",paste("Neuron",idxn), fill=cols, bty='n')
if(savePDF) dev.off()

  isi = lapply(trials[idxn], function(x) 1e3*unlist(lapply(x,diff)))


if(savePDF) pdf(paste0(plotlib,"3neurons_histograms.pdf"), width=5, height=4)
  par(mfrow=c(1,1), mar=c(3,3,1,1), oma=c(0,0,0,0))
      col1 = add.alpha('dodgerblue',.75)
      col2 = add.alpha('red',.5)
      col3 = add.alpha('springgreen4',.5)
      hst3 = hist(isi[[3]], breaks=c(seq(0,1100,5),1e6), xlim=c(0,250), border=NA, col=col3, prob=TRUE, ylim=c(0,.01), ann=FALSE, axes=FALSE)
      par(new=TRUE)
      hst2 = hist(isi[[2]], breaks=c(seq(0,1100,5),1e6), xlim=c(0,250), border=NA, col=col2, prob=TRUE, ylim=c(0,.01), ann=FALSE)
      par(new=TRUE)
      hst1 = hist(isi[[1]], breaks=c(seq(0,1100,5),1e6), xlim=c(0,250), border=NA, col=col1, prob=TRUE, ylim=c(0,.01), ann=FALSE, axes=FALSE)

      mtext("ISI distributions", 3, line=0)
      mtext("ms", 1, line=2)
      mtext("density", 2, line=2)
      legend("topright",paste("Neuron",idxn), fill=cols, bty='n')
if(savePDF) dev.off()


    # Load all results for the two specific neurons:
    turtle.res = list()
    for(i in 1:length(idxn)){
      load(paste0(datalib,fns[idxn[i]]))
      turtle.res[[i]] = tmpres
    }




if(savePDF) pdf(paste0(plotlib,"3neurons_trials.pdf"), width=10, height=8)
      # if(savePDF) pdf(paste0(plotlib,"2neurons_trials_125.pdf"), width=10, height=12)
      layout(1)
      par(mar=c(3,3,1,1), mfrow=c(1,1))

      cols = c('dodgerblue','red','springgreen4')
      # layout(matrix(c(1,1,2,3),nc=2))
      plot(0,0, type='n', xlim=c(0*1e3,40*1e3), ylim=c(1,10), yaxt='n', bty='n', ann=FALSE, bty='n')
      axis(side = 2, at = 1:10, labels=10:1, tick = FALSE, las=1)

      offs = c(.2,0,-.2)
      i=j=1
      for(i in 1:length(idxn)){
        axis(side = 2, at = 1:10+offs[i], labels = rep(paste0("(",idxn[i],")"),10), tick = FALSE, las=1, cex.axis=.5, line = -1.25)
        for(j in 1:10){
          tmpres = turtle.res[[i]][[j]]
          if(!is.null(tmpres) & !all(is.na(tmpres))){
            x = tmpres$spk
            tmp = (11-j)*(tmpres$prob>0.5)
            tmp[tmp==0] = NA
            segments(0,11-j+offs[i],4e4,11-j+offs[i], col=add.alpha(cols[i],.15), lwd=5, lend=1)
            lines(tmp+offs[i], col=add.alpha(cols[i],.9), lwd=5, lend=1)
            points(x,rep(11-j+offs[i],length(x)), pch=124, col=add.alpha('black',.45))
          }
        }
      }

      mtext("Trial (neuron)", 2,line=2)
      mtext("Time from stimulus (ms)", 1,line=2)
      # abline(h=1:10, lty=3)
if(savePDF) dev.off()




if(savePDF) pdf(paste0(plotlib,"3neurons.pdf"), width=10, height=10)

      par(mar=c(4,4,1,1), oma=c(0,0,0.5,0), cex.axis=1.2)
      layout(matrix(c(1,2,3,1,2,4), nr=3), heights = c(20,1,10))


      cols = c('dodgerblue','red','springgreen4')
      # layout(matrix(c(1,1,2,3),nc=2))
      plot(0,0, type='n', xlim=c(0*1e3,40*1e3), ylim=c(1,10), yaxt='n', bty='n', ann=FALSE, bty='n')
      axis(side = 2, at = 1:10, labels=10:1, tick = FALSE, las=1, line = -1)


      offs = c(.2,0,-.2)
      i=j=1
      for(i in 1:length(idxn)){
        # axis(side = 2, at = 1:10+offs[i], labels = rep(paste0("(",idxn[i],")"),10), tick = FALSE, las=1, cex.axis=.5, line = -1.25)
        for(j in 1:10){
          tmpres = turtle.res[[i]][[j]]
          if(!is.null(tmpres) & !all(is.na(tmpres))){
            tmp = (11-j)*(tmpres$prob>0.5)
            tmp[tmp==0] = NA
            segments(0,11-j+offs[i],4e4,11-j+offs[i], col=add.alpha(cols[i],.15), lwd=5, lend=1)
            lines(tmp+offs[i], col=add.alpha(cols[i],.9), lwd=5, lend=1)
          }
          # x = tmpres$spk
          x = 1e3*trials[[ idxn[i] ]][[j]]
          points(x,rep(11-j+offs[i],length(x)), pch=124, col=add.alpha('black',.65))

        }
      }
      mtext("A",side=3, line=0, adj=0, family='sans', cex=1.25)
      # mtext("Trial (neuron)", 2,line=2)
      mtext("Trial", 2,line=2)
      mtext("ms", 1,line=2.2)
      # abline(h=1:10, lty=3)
      abline(h=1:9+.5, lty=3)

      par(mar=c(0,0,0,0))
      plot(0,0, type='n', bty='n', axes=FALSE, ann=FALSE)
      legend("top",paste("Neuron",idxn), fill=cols, bty='n', horiz = TRUE, cex = 1.25)

      par(mar=c(5,4,2,1))
      isi3 = lapply(trials[idxn], function(x) 1e3*unlist(lapply(x,diff)))

      col1 = add.alpha('dodgerblue',.75)
      col2 = add.alpha('red',.5)
      col3 = add.alpha('springgreen4',.5)
      hst3 = hist(isi3[[3]], breaks=c(seq(0,1100,5),1e6), xlim=c(0,250), border=NA, col=col3, prob=TRUE, ylim=c(0,.015), ann=FALSE, axes=FALSE)
      par(new=TRUE)
      hst2 = hist(isi3[[2]], breaks=c(seq(0,1100,5),1e6), xlim=c(0,250), border=NA, col=col2, prob=TRUE, ylim=c(0,.015), ann=FALSE)
      par(new=TRUE)
      hst1 = hist(isi3[[1]], breaks=c(seq(0,1100,5),1e6), xlim=c(0,250), border=NA, col=col1, prob=TRUE, ylim=c(0,.015), ann=FALSE, axes=FALSE)
      mtext("B",side=3, line=0, adj=0, family='sans', cex=1.25)
      mtext("ms", 1, line=2)
      mtext("density", 2, line=2)
      # legend("topright",paste("Neuron",idxn), fill=cols, bty='n')

      cols = c('dodgerblue','red', 'springgreen4')
      col1 = add.alpha(cols[1],.75)
      col2 = add.alpha(cols[2],.75)
      col3 = add.alpha(cols[3],.75)
      # matplot(exp(Kb.idx.trials), type='l', lty=3, lwd=1, col=c(rep(col1,10),rep(col2,10),rep(col3,10)), ylim=c(0,.015), xlab="ms", ylab=paste0("Trials for neurons (",paste(idxn, collapse = ","),")"), bty='n')
      # par(new=TRUE)
      matplot(exp(Kb.idx), type='l', lty=1, lwd=2, col=c(col1,col2,col3), ylim=c(0,.015), xlab="", ylab="", bty='n')
      mtext("C",side=3, line=0, adj=0, family='sans', cex=1.25)
      mtext("ms", 1, line=2)
if(savePDF) dev.off()






