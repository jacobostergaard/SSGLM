misc::clean_up()

library(SSGLM)
datalib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"


fn = paste0(datalib,"izh_sim.Rda") # load Izhikevich simulation and mpfpp settings
load(fn)

fn = paste0(datalib,"izh_mpfpp_res.Rda") # load Izhikevich parameter estimates from mpfpp algorithm
load(fn)


# Simulate a SSGLM neuron using estimates bases on Izhikevich simulation
    set.seed(1234)
    sim_dt = izh$dt
    nsim   = izh$nsim
    xpars  = izh$xpars
    # kts = c(5,10)/sim_dt
    # kts[kts> 20/sim_dt] = 20/sim_dt # Knots must be truncated to interval [0,100]
    # kts       = sort(c(0,0,0,0,kts,20/sim_dt))
    kts    = izh$kts
    S      = splines::splineDesign(kts, x=1:max(kts), ord=4, outer.ok=TRUE)
    Bb     = as.numeric(izh.res$P[2:4])
    b0     = as.numeric(izh.res$P[1])
    bk     = as.numeric(S%*%Bb)
    r0     = as.numeric(izh.res$P[5])
    rk     = 0

    sim     = as.data.frame(sim_ssglm(nsim, r0, rk, b0, bk, c(xpars$bst,xpars$rst), sim_dt))

    colnames(sim) = c("y","x","lam")
    spk = which(sim$y>0)
    isi = diff(spk)

    glm = list(sim=sim, spk=spk, isi=isi, simpars = list(B0.bst = b0, B.bst = Bb, B0.rst = r0, B.rst = rk), xpars = xpars, kts = kts, S=S, nsim=nsim, dt=sim_dt)


  fn = paste0(datalib,"ssglm_sim.Rda")
  save(glm,file=fn)

  # Compare simulations from histograms and quantiles
  layout(1:2)
  hist(izh$isi, breaks=seq(0,10000,10), border=NA, col=add.alpha('black',.25), ylim=c(0,0.01), prob=TRUE, ann=FALSE, xlim=c(0,4000))
  par(new=TRUE)
  hist(glm$isi, breaks=seq(0,10000,10), border=NA, col=add.alpha('red',.45), ylim=c(0,0.01), prob=TRUE, ann=FALSE, xlim=c(0,4000))


  x = quantile(izh$isi, probs=seq(0,1,.1))
  y = quantile(glm$isi, probs=seq(0,1,.1))

  plot(x, y, xlab="Izh", ylab="GLM", col=add.alpha('red',.75), pch=16, bty='n')
  abline(0,1, lty=3)
  text(x,y,labels = seq(0,100,10), cex=.5, pos = 3, offset = .25)

  length(izh$spk)
  length(glm$spk)

  table(glm$sim$y)


