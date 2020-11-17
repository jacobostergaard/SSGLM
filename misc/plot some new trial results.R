
par(mfrow=c(5,5), mar=c(2,2,1,1), oma=c(0,0,0,0))

for(i in 1:length(trials)){
  tmp = unlist(lapply(trials[[i]], diff))
  hist(tmp, breaks=c(seq(0,1e4,10),1e5), prob=TRUE, xlim=c(0,5e2), main=i, ylim=c(0,.02))
}

i=3
par(mfrow=c(2,5))
for(j in 1:10){
  tmp = diff(trials[[i]][[j]])
  hist(tmp, breaks=c(seq(0,1e4,10),1e5), prob=TRUE, xlim=c(0,5e2), main=j, ylim=c(0,.02))
}



i=3; j=1
spk = trials[[i]][[j]]

plot(spk, rep(1,length(spk)), pch=124, xlim=c(0,4e4))
abline(v=4e4/3)

y         = as.numeric(SSGLM::make_y(spk,1,N = 40e3, verbose = FALSE)$y)

range(unlist(lapply(trials, function(x) lapply(x, max))))

