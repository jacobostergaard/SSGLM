misc::clean_up()
library(SSGLM)
datalib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/plots/"
plotlib = "/Users/jacob/Library/Mobile Documents/com~apple~CloudDocs/GitHub/TeX/SSGLM plos-latex-template/gfx/"
set.seed(1234)


fns = list.files(datalib)
fns = fns[grepl("turtle_n",fns) & !grepl("trial",fns) & !grepl("qts",fns)]
idx = as.numeric(substr(fns,9,nchar(fns)-4))

length(idx)

neuron = idx[9]

plot_neuron <- function(neuron, trial=1){

    fn = paste0(datalib,"turtle_n",neuron,".Rda")
    load(fn)

    res = tmpres[[trial]]

    Sb = res$Sb
    Sr = res$Sr
    P = res$est$P
    idxb = grepl("b",names(P))
    idxr = grepl("r",names(P))

    Pb = as.numeric(P[idxb])
    Pr = as.numeric(P[idxr])

    SEb = sqrt(diag(res$est$cov[idxb,idxb]))
    SEr = sqrt(diag(res$est$cov[idxr,idxr]))

    Pb.hi = Pb+1.96*SEb
    Pb.lo = Pb-1.96*SEb

    Pr.hi = Pr+1.96*SEr
    Pr.lo = Pr-1.96*SEr
    Pb = cbind(Pb,Pb.lo, Pb.hi)
    Pr = cbind(Pr,Pr.lo, Pr.hi)

    kb = as.data.frame(Sb%*%Pb[-1,])
    kr = as.data.frame(Sr%*%Pr[-1,])

    kb.est1 = exp(kb[,1]+Pb[1,1])
    kb.lo1  = exp(kb[,2]+Pb[1,2])
    kb.hi1  = exp(kb[,3]+Pb[1,3])

    kb.est2 = exp(kb[,1])
    kb.lo2  = exp(kb[,2])
    kb.hi2  = exp(kb[,3])

    kr.est1 = exp(kr[,1]+Pr[1,1])
    kr.lo1  = exp(kr[,2]+Pr[1,2])
    kr.hi1  = exp(kr[,3]+Pr[1,3])

    kr.est2 = exp(kr[,1])
    kr.lo2  = exp(kr[,2])
    kr.hi2  = exp(kr[,3])

    # par(mar=c(3,3,1,1), mfrow=c(2,2), bty='n')

    plot(kb.est1, ylim=c(0,max(kb.hi1)), type='n')
    polygon(c(1:length(kb.est1),length(kb.est1):1), c(kb.hi1,rev(kb.lo1)), border=NA, col=add.alpha('dodgerblue',.25))
    lines(kb.est1, col=add.alpha('red',.85), lwd=2)

    plot(kb.est2, ylim=c(0,max(kb.hi2)), type='n')
    polygon(c(1:length(kb.est2),length(kb.est2):1), c(kb.hi2,rev(kb.lo2)), border=NA, col=add.alpha('dodgerblue',.25))
    lines(kb.est2, col=add.alpha('red',.85), lwd=2)


    plot(kr.est1, ylim=c(0,max(kr.hi1)), type='n')
    polygon(c(1:length(kr.est1),length(kr.est1):1), c(kr.hi1,rev(kr.lo1)), border=NA, col=add.alpha('dodgerblue',.25))
    lines(kr.est1, col=add.alpha('red',.85), lwd=2)

    plot(kr.est2, ylim=c(0,max(kr.hi2)), type='n')
    polygon(c(1:length(kr.est2),length(kr.est2):1), c(kr.hi2,rev(kr.lo2)), border=NA, col=add.alpha('dodgerblue',.25))
    lines(kr.est2, col=add.alpha('red',.85), lwd=2)

    if(!is.null(res)){
      cat("Trial",trial,"has pval", round(res$est$ks$p.value,3))
    } else{
      cat("Trial",trial,"has pval NULL")
    }
    cat("\n")
}


# for(i in 1:length(idx)){
neuron = idx[6]
  j=2
  for(j in 1:10){
    par(mar=c(3,3,1,1), oma=c(0,0,1,0), mfrow=c(2,2), bty='n')
    plot_neuron(neuron,j)
    txt = paste("Neuron",neuron,"trial",j)
    mtext(txt,3, outer=TRUE, line=-.5)
    Sys.sleep(1)
  }
# }

  neuron = idx[10];j=2
  par(mar=c(3,3,1,1), oma=c(0,0,1,0), mfrow=c(2,2), bty='n')
  plot_neuron(neuron,j)
  txt = paste("Neuron",neuron,"trial",j)
  mtext(txt,3, outer=TRUE, line=-.5)
