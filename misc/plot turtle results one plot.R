misc::clean_up()
library(SSGLM)
datalib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/plots/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/TeX/SSGLM plos-latex-template/gfx/"

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

  trials = list(n1=get_spiketimes(12),n2=get_spiketimes(13))

  rmv = c(3,8,9) # remove these trials: either little to no spikes, or very irregular behavior
  trials$n1 = trials$n1[-rmv]
  trials$n2 = trials$n2[-rmv]

  # subtract initial 10s quiescent period (stimulus at 10s) and convert to milliseconds
  trials$n1 = lapply(trials$n1, function(x) (x-10)*1e3)
  trials$n2 = lapply(trials$n2, function(x) (x-10)*1e3)


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

    trials$n1 = lapply(trials$n1,downsample)
    trials$n2 = lapply(trials$n2,downsample)

    isi = unlist( lapply(trials, function(x) lapply(x,diff)) )

    maxb = 250 # choose max isi for spikes in a burst

# Find all burst lenghts for all trials for both neurons
    bst = numeric(0)
    for(i in 1:7){
      tmp1 = burst_info(trials$n1[[i]],cutoff = maxb, refrac = min(isi))
      tmp2 = burst_info(trials$n1[[i]],cutoff = maxb, refrac = min(isi))
      bst = c(bst, tmp1$burstlength, tmp2$burstlength)
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

    all.res = list(n1=list(),n2=list())
    for(i in 1:2){
      for(j in 1:7){
        fn = paste0(datalib,"turtle_neuron",i,"_trial",j,".Rda")
        load(fn)
        all.res[[i]][[j]] = res
        cat("\nLoaded neuron",i,"trial",j)
      }
    }
    res = all.res

    maxk  = nrow(res$n1[[1]]$S)

    ks.results <- function(neuron,trial, P, k){
      res  = res[[neuron]][[trial]]
      p    = length(P)
      L    = get_intensity(as.numeric(res$y), res$X, B0.bst = P[1], k.bst = k, B0.rst = P[p],k.rst = NULL)
      z    = rescale_waiting_times(res$spk, L, 1)
      ks   = suppressWarnings(ks.test(z, 'pexp'))

      return(list(z=z, pval = ks$p.value))
    }

    col1 = add.alpha('black',.75)
    col2 = add.alpha('red',.75)
    col3 = add.alpha('dodgerblue',.25)
    col4 = add.alpha('dodgerblue',.75)

fn = paste0(plotlib, "turtle_res.pdf")
pdf(fn, height=12, width=9)

    idx = 1:10000
    par(mar=c(3,3,2,2), oma=c(1,1.5,1,1), bty='n', mfrow=c(3,2))

# Plot fit of transition parameters
    plot(0,0, type='n', xlim=c(0,1000), ylim=c(0,1))
        lines(idx,getF(bst.h)$y[idx], col=col3, lty=1, lwd=1)
        points(idx,getF(bst.h)$y[idx], col=col1, cex=.5, pch=16)
        lines(idx,Z.cdf(idx, bst.pars[1],bst.pars[2])$cdf, col=col2, lwd=3)
        legend("topleft", bty='n', legend=c("Observed","Fitted CDF"), pch=c(16,NA), lty=c(NA,1), pt.cex = c(1.25,NA), lwd=c(NA,3), col=c(col1,col2))
        mtext("CDF", side = 2, line=2)
        mtext("Burst length [ms]", side=1, line=2.25)
        mtext("A", side = 3, line=0, adj=0, family = 'sans', cex=1.25)


    plot(0,0, type='n', xlim=c(0,2000), ylim=c(0,1))
        lines(idx,getF(rst.h)$y[idx], col=col3, lty=1, lwd=1)
        points(idx,getF(rst.h)$y[idx], col=col1, cex=.5, pch=16, xlim=c(1000,4000))
        lines(idx,Z.cdf(idx, rst.pars[1],rst.pars[2])$cdf, col=col2, lwd=3)
        legend("topleft", bty='n', legend=c("Observed","Fitted CDF"), pch=c(16,NA), lty=c(NA,1), pt.cex = c(1.25,NA), lwd=c(NA,3), col=c(col1,col2))
        mtext("CDF", side = 2, line=2)
        mtext("Rest length [ms]", side=1, line=2.25)
        mtext("B", side = 3, line=0, adj=0, family = 'sans', cex=1.25)

# Plot the kernel estimates

    P1 = P2 = 0
    plot(0,0, type='n', xlim=c(0,maxk), ylim=c(0,2), bty='n')
        for(i in 1:7){
          tmpres = res$n1[[i]]
          tmpP1  = as.numeric(tmpres$est$P)
          P1     = P1+tmpP1/7
          lines(exp(tmpres$est$Kb), col=add.alpha(col2,.75), lwd=1, lty=3)

          tmpres = res$n2[[i]]
          tmpP2  = as.numeric(tmpres$est$P)
          P2     = P2 + tmpP2/7
          lines(exp(tmpres$est$Kb), col=add.alpha(col4,.75), lwd=1, lty=3)
        }

        S1 = res$n1[[1]]$S
        S2 = res$n2[[1]]$S
        p1 = length(P1)
        p2 = length(P2)

        k1 = S1%*%P1[2:(p1-1)]
        k2 = S2%*%P2[2:(p2-1)]

        lines(exp(k1), col=add.alpha(col2,.85), lwd=2)
        lines(exp(k2), col=add.alpha(col4,.85), lwd=2)
        mtext(expression(e^{S~beta}), side = 2, line=2)
        mtext("History lag [ms]",1, outer=FALSE, line=2.25, cex=1)
        # legend("topright", c("Neuron 1", "Neuron 2", "average", "individual"), col=c(col2,col4,col1,col1), lty=c(NA,NA,1,3), lwd=c(NA,NA,2,1), pch=c(15,15,NA,NA),bty='n')
        legend("topright", c("Neuron 1: Average","Neuron 1: Trials", "Neuron 2: Average","Neuron 2: Trials"), col=c(col2,col2,col4,col4), lty=c(1,3,1,3), lwd=2, bty='n')

        mtext("C", side = 3, line=0, adj=0, family = 'sans', cex=1.25)


    plot(0,0, type='n', xlim=c(0,maxk), ylim=c(0,.03), bty='n')

        for(i in 1:7){
          tmpres = res$n1[[i]]
          lines(exp(tmpres$est$Kb+tmpres$est$B0b), col=add.alpha(col2,.75), lwd=1, lty=3)
          tmpres = res$n2[[i]]
          lines(exp(tmpres$est$Kb+tmpres$est$B0b), col=add.alpha(col4,.75), lwd=1, lty=3)
        }

        lines(exp(k1+P1[1]), col=add.alpha(col2,.85), lwd=2)
        lines(exp(k2+P2[1]), col=add.alpha(col4,.85), lwd=2)
        mtext(expression(e^{S~beta+beta[0]}), side = 2, line=2)
        mtext("History lag [ms]",1, outer=FALSE, line=2.25, cex=1)
        # legend("topright", c("Neuron 1", "Neuron 2", "average", "individual"), col=c(col2,col4,col1,col1), lty=c(NA,NA,1,3), lwd=c(NA,NA,2,1), pch=c(15,15,NA,NA),bty='n')
        legend("topright", c("Neuron 1: Average","Neuron 1: Trials", "Neuron 2: Average","Neuron 2: Trials"), col=c(col2,col2,col4,col4), lty=c(1,3,1,3), lwd=2, bty='n')
        mtext("D", side = 3, line=0, adj=0, family = 'sans', cex=1.25)


        z     = list(n1=numeric(),n2=numeric())
        P     = list(P1, P2)
        k     = list(k1, k2)
        pvals = matrix(nr=7,nc=2)
        for(i in 1:2){
          for(j in 1:7){
            tmp  = ks.results(i,j,P[[i]],k[[i]])
            z[[i]] = c(z[[i]],tmp$z)
            pvals[j,i] = tmp$pval
          }
        }

    # Plot the KS tests
    ks_plot(z$n1, a = .05, col1 = col2, col2 = col3)
        legend("topleft", c("Identity line", "Quantiles","95% CI limit"), lty=c(1,1,1), col=c(col1,col2,col3), lwd=c(2,2,8), bty='n')
        mtext("E", side = 3, line=0, adj=0, family = 'sans', cex=1.25)

    ks_plot(z$n2, a = .05, col1 = col2, col2 = col3)
        legend("topleft", c("Identity line", "Quantiles","95% CI limit"), lty=c(1,1,1), col=c(col1,col2,col3), lwd=c(2,2,8), bty='n')
        mtext("F", side = 3, line=0, adj=0, family = 'sans', cex=1.25)

dev.off()
