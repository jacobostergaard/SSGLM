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


# input settings are simulation values
B0.bst  = 0# glm$simpars$B0.bst
B.bst   = c(-10, 0, 0) #glm$simpars$B.bst
B0.rst  = -50#glm$simpars$B0.rst
B.rst   = NULL # glm$simpars$B.rst
M       = 100
slag    = 2*nrow(glm$S)
q       = 1e-8
pars    = make_mpfpp_par_input(Bb0 = B0.bst, Br0 = B0.rst, S = glm$S, Bb = B.bst, Br = B.rst, xpar = glm$xpars, M, q0 = q, s_lag = slag, dt = glm$dt)
pars$q  = c(1e-9,q,q,q,0)

nobs    = glm$nsim

Sys.time()
system.time({

  y = glm$sim$y[1:nobs]
  x = glm$sim$x[1:nobs]

  # res = mpfpp(y,x,pars, usetrue = TRUE, verbose = TRUE)
  # res$y = glm$sim$y[1:nobs][1:nobs]
  # fn  = paste0(datalib,"glm_mpfpp_truepar_truex.Rda")
  # save(res, file = fn)
  #
  # res = mpfpp(y,x,pars, usetrue = FALSE, verbose = TRUE)
  # res$y = glm$sim$y[1:nobs][1:nobs]
  # fn  = paste0(datalib,"glm_mpfpp_truex.Rda")
  # save(res, file = fn)

  x = sample(c(-1,1),M,replace=TRUE)

  # res = mpfpp(y,x,pars, usetrue = TRUE, verbose = TRUE)
  # res$y = glm$sim$y[1:nobs][1:nobs]
  # fn  = paste0(datalib,"glm_mpfpp_truepar.Rda")
  # save(res, file = fn)

  res = mpfpp(y,x,pars, usetrue = FALSE, verbose = TRUE)
  res$y = glm$sim$y[1:nobs][1:nobs]
  fn  = paste0(datalib,"glm_mpfpp.Rda")
  save(res, file = fn)

  glm.res = get_est(res)
  fn  = paste0(datalib,"glm_mpfpp_res.Rda")
  save(glm.res, file=fn)

})


