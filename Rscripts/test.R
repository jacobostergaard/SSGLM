library(Matrix)
misc::clean_up()

sparseVector( i = 1:10, length=20)

set.seed(1234)

trials <- function(n){

  S = matrix(ncol=2, nrow=0)
  for(i in 1:n){
    Ns = sample(10,1)
    spk = sort(sample(19,Ns))
    S = rbind(S,cbind(rep(i, length(spk)),spk))
  }

  M = sparseMatrix(i = S[,1], j=S[,2], x=rep(1,nrow(S)))

  return(M)
}

L = 1
y = trials(L)

kts       = c(0,0,0,0,5,10,10)
S         = splines::splineDesign(kts, x=1:max(kts), ord=4, outer.ok=TRUE)

res = mpfpp(y = y, x = matrix(1,1,1),
                  Bb = c(.5,1,2,3), Br = 1, Sb = S, Sr = S,
                  xpars = 1:4, s_lag = 2, Neffpct = 1,
                  M = 6, w0 = .1, q = .05, dt = 1,
                  usetrue = FALSE, verbose = TRUE)

res



