library(misc)
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
hist(1000*isi, breaks=1e4, xlim=c(0,200))
1000*range(isi)
100*sum(isi<1e-3)/length(isi)

# apply(spks_in_trials,2,sum)
par(mfrow=c(2,1))
idx = order(apply(spks_in_trials,1,median), decreasing = TRUE)
boxplot(t(spks_in_trials[idx,]))

idx = order(apply(spks_in_trials,2,median), decreasing = TRUE)
boxplot(spks_in_trials[,idx])


par(mar=c(3,3,1,1), oma=c(1,1,1,1), bty='n')

layout(matrix(c(1,2,3,4,4,4), nr=3), width=c(1,2))
plot(apply(spks_in_trials,2,median))
plot(apply(spks_in_trials,2,function(x) diff(quantile(x, c(.25,.75)))))
plot(apply(spks_in_trials,2,sd))


i=j=1
plot(0,0, type='n', xlim=c(0,400), ylim=c(0,250), bty='n', xlab='ms', axes=FALSE, ylab="neurons")
for(i in 1:249){
  for(j in 1:10){
    tmp = trials[[i]][[j]]
    points(tmp+(j-1)*40, rep(i, length(tmp)), pch='.')
  }
}
abline(v=(0:10)*40, lty=3)
text((0:9)*40+20, rep(251,10), labels = as.character(1:10), pos = 3)





plot(apply(spks_in_trials[,-9], 1, mean))
tmp = apply(spks_in_trials[,-9], 1, mean)
quantile(tmp)
# hist(tmp, breaks=100)

# hist(tmp, breaks=c(seq(0,500,10)))

# quantile(tmp,c(.25,.75))

idx = which(tmp<100 & tmp>25)

trials2 = trials[idx]
layout(1)
i=j=1
par(oma=c(0,0,0,0), mar=c(1,1,1,1))
plot(0,0, type='n', xlim=c(0,400), ylim=c(0,255), bty='n', xlab='ms', axes=FALSE, ylab="neurons")
# for(i in 1:length(idx)){
for(i in 1:249){
  for(j in 1:10){
    tmp = trials[[i]][[j]]
    # Col = ifelse(i %in% idx, "red","blue")
    Col = add.alpha('black',.75)
    points(tmp+(j-1)*40, rep(i, length(tmp)), pch='.', col=Col)
  }
}
abline(v=(0:10)*40, lty=3)
text((0:9)*40+20, rep(251,10), labels = as.character(1:10), pos = 3, cex=.75)


diff(round(1000*trials[[1]][[1]]))
