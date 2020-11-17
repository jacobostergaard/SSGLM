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

ids = 1:50

neurons = list()
for(i in ids){
  neurons[[i]] = get_spiketimes(i)
}
names(neurons) = paste0("n",ids)

trs = 1:10
trials = list()
for(i in trs){
  trials[[i]] = lapply(neurons, function(x) x[[i]])
}

names(trials) = paste0("t",trs)

plot_trial <- function(i=1){
  spk = trials[[ids[i]]]
  plot(0,0, xlim=c(0,40), ylim=c(1,length(ids)), type='l', bty='n')
  for(i in 1:length(ids)){
    tmp = spk[[ids[i]]]
    points(tmp, rep(i, length(tmp)), pch='.')
  }
  invisible(spk)
}

spk = plot_trial(1)
isi = numeric(0)
spks = numeric(0)
for(i in trs){
  spk = plot_trial(i)
  spk = lapply(spk, sort)
  spks = c(spks,unlist(spk))
  isi = c(isi,unlist(lapply(spk, diff)))
}
layout(1)

hist(isi, breaks=seq(0,200, .01), xlim=c(0,.1))

hist(isi, breaks=c(seq(0,.005, .0001),200), xlim=c(0,.005))

1000*4*0.000025


ms = 5
sum(isi<=4*min(isi)*ms)

mean(isi<=4*min(isi)*ms)*length(isi)
round(100*mean(isi<=4*min(isi)*ms),2)
length(isi)
min(isi)/2.5e-5


plot(0,0, type='n', xlim=c(0,250), ylim=c(0,1))
for(ms in 1:250){
  points(ms,mean(isi<=4*min(isi)*ms), pch='.')
}


prbs = seq(0,.99,length=100)
prb= .95

plot(quantile(isi, probs = prbs), prbs, type='l')
abline(h=prb, v=quantile(isi, probs = prb), col='red', lty=3)
hist(isi[isi<quantile(isi,probs=prb)], breaks=100)
plot(isi, pch='.', ylim=c(0,quantile(isi,probs=prb)+.1))
plot(isi, pch='.', ylim=c(0,1))
abline(h=quantile(isi,probs=prb), col='red')


hist(isi[isi>quantile(isi,probs=prb)], breaks=100)



i=27
isi_i = unlist(lapply(lapply(neurons[[i]],sort),diff))
hist(isi_i, breaks=100)


bst = burst_info(1000*spks,cutoff = 1000)
mean(bst$burstlength>10*1e3)
mean(bst$no_spikes==1)
# idx = (bst$no_spikes!=1) & (bst$burstlength<10*1e3)
tmp = sort(bst$burstlength[idx])
tmp = sort(bst$burstlength)
range(tmp)/1e3
plot(tmp,(1:length(tmp))/length(tmp), type='l', xlim=c(0,max(tmp)))

z= 1:10000
a = tmppars$pars$bst[1]
b = tmppars$pars$bst[2]
a = 1e9
b = .5
tmpz = Z.cdf(z, a=tmppars$pars$bst[1], b=tmppars$pars$bst[2])
lines(tmpz$z, tmpz$cdf, col='red')


tmppars = find_xpars(1e3*spks, 1)
tmppars$pars
lines(tmppars$bst*1e-3, col='red')
plot(tmppars$bst, type='l')



spk = plot_trial(6)
