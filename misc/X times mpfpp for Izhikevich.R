misc::clean_up()
library(SSGLM)
datalib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/plots/"

# Simulate Izhikevich neuron and fit transition parameters for X and SMC

# Simulate Izhikevich neuron
set.seed(1234)
sim_dt = 0.01

# Simulate Izhikevich neuron
nsim = 5*1000/sim_dt # 20 seconds of observations
izh_par = list(V0 = -70, U0 = 10, sig=20, a = 0.02, b = 0.2, c = -50, d = 2, I = 20, dt = sim_dt, type = "tonic bursting")

izh     = sim_izhikevich(nsim, izh_par)
izh_spk = izh$spk
izh_isi = diff(which(izh_spk>0))

cat("\nIzhikevich:",nsim,"simulations with",sum(izh_spk),"spikes")


izh = list(sim=izh)
izh$spk  = which(izh$sim$spk>0)
izh$isi  = izh_isi
izh$dt   = sim_dt
izh$par  = izh_par
izh$nsim = nsim

spk    = izh$spk
cutoff = 15/izh$dt
refrac = min(izh$isi)
init   = c(20000,2.5,4000,10)
xpars  = SSGLM::find_xpars(spk, cutoff, refrac, init)$pars
izh$xpars = xpars






izh.ks= list()


f <- function(i){

  kts = c(5,10)/sim_dt
  kts[kts> 20/sim_dt] = 20/sim_dt # Knots must be truncated to interval [0,100]
  kts       = sort(c(0,0,0,0,kts,20/sim_dt))
  S         = splines::splineDesign(kts, x=1:max(kts), ord=4, outer.ok=TRUE)
  izh$kts = kts


  B0.bst  = 0
  B.bst   = c(-15,0,0)
  B0.rst  = -50
  B.rst   = NULL
  M       = 100
  slag    = 2*nrow(S)
  q       = 1e-8

  pars       = make_mpfpp_par_input(Bb0 = B0.bst, Br0 = B0.rst, S = S, Bb = B.bst, Br = NULL, xpar = xpars, M, q0 = q, s_lag = slag, dt = izh$dt)

  pars$q = c(1e-9,q,q,q,0)


  nobs = 500e3

  y = izh$sim$spk[1:nobs]
  x = sample(c(-1,1), M, replace = TRUE)

  cat("\nRunning KS iteration ",i)
  system.time({
    res = mpfpp(y,x,pars, usetrue = FALSE, verbose = TRUE)
    res$y = y

    izh.res = get_est(res)

    fn = paste0(datalib, "izh_ks.Rda")
    load(fn)
    nks = length(izh.ks)
    izh.ks[[nks+1]] = izh.res$ks
    save(izh.ks,file=fn)

    })
  return(izh.ks)
  }


nks = 50


iks = list()
for(i in 1:nks){
  Sys.time()
  iks[[i]] = f(i)
}

izh.ks = iks
fn = paste0(datalib, "izh_ks.Rda")
save(izh.ks,file=fn)
# parallel::mclapply(1:2, f, mc.cores = 2)

hist(unlist(lapply(iks, function(x) x$statistic)), breaks=20, border=NA, col=add.alpha('red',.5))
abline(v=ks_crit(.05)/sqrt(length(spk)), lty=3)
mean(unlist(lapply(iks, function(x) x$p.value))<0.05)
