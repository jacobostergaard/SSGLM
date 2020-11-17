misc::clean_up()

library(SSGLM)
datalib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"

fn = paste0(datalib,"ssglm_sim.Rda")

load(fn)


layout(1:2)
hist(glm$isi, breaks=1000, border=NA, col=add.alpha('black',.5))

par(mar=c(3,0,0,0))
# Plot Izhikevich spike train
nMax  = min(500000,glm$nsim)
splot = matrix(glm$spk[1:nMax],nr=20,byrow=TRUE)

plot(0,0,type='n', xlim=c(0,nMax/nrow(splot)), ylim=c(1,nrow(splot)), bty='n', ann=FALSE, cex.axis=1.15, yaxt='n')
for(i in 1:nrow(splot)){
  # Plot the spikes
  tmps = which(splot[i,]>0)
  # tmpy = tmps[seq(1,length(tmpy),2)] # Downsample the plotting sequence
  points(tmps, rep(nrow(splot)+1-i-.1,length(tmps)), pch=124, col=add.alpha('black',.75), cex=.8)
}
# mtext("ms",side = 1,line=2, cex=1.15)
segments(0,nrow(splot)+.5,20/glm$dt,nrow(splot)+0.5, lwd=10, lty = 1, lend='butt', col=add.alpha('black',.5))



