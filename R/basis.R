####################################################################################
##
## Create basis functions, either splines, raised cosines or indicator based.
##
####################################################################################

## Create splines matrix
    make_splines <- function(lag, knots = c(5,10,20,50), tension=0, pct=FALSE){
      # knots = qts
      # lag = 100; knots = c(5,10,20,50)

      # input knots are in fractions of the lag lengths
      if(pct){
        if(max(knots)>1) knots = knots/lag
        ctrl_pts = sort(quantile(0:lag, probs = knots))
      } else{
        ctrl_pts = sort(knots)
      }

      # # we need the first knot at 0
      # if(ctrl_pts[1] > 0){
      #   ctrl_pts = trunc(c(0,ctrl_pts))
      # }
      #
      # # we need the last knot at 'lag'
      # if(max(ctrl_pts) < lag){
      #   ctrl_pts = trunc(c(ctrl_pts,lag))
      # }

      # each spline endpoint is a derivative, hence auxillary points are needed in both ends (if not given)
      if(max(ctrl_pts)==lag) ctrl_pts = trunc(c(ctrl_pts,lag+10))
      if(min(ctrl_pts)==0)   ctrl_pts = trunc(c(-10,ctrl_pts))

      s = tension # tension parameter
      S = matrix(0, nr=lag, nc=length(ctrl_pts))

      for(i in 1:lag){
        idx      = max(which(ctrl_pts < i))
        this_pt  = ctrl_pts[idx]
        next_pt  = ctrl_pts[idx+1]
        next2_pt = ctrl_pts[idx+2]

        u = (i-this_pt)/(next_pt-this_pt)
        l = (next2_pt-next_pt)/(next_pt-this_pt)

        tmp = matrix(c(  -s, 2-s/l,   s-2,  s/l,
                        2*s, s/l-3, 3-2*s, -s/l,
                         -s,     0,     s,    0,
                          0,     1,     0,    0), nc=4, byrow=TRUE)

        p = crossprod(c(u^3, u^2, u, 1) ,tmp)

        S[i, (idx-1):(idx+2)] = p
      }

      # Remove columns of zeros!
      f <- function(s) all(s==0)
      if(any(apply(S,2,f))) S = S[,-which(apply(S,2,f))]

      return(S)

    }


####################################################################################

## Create raised cosine basis functions

    make_cosines <- function(lag, num_bases, spacing, t_min=NULL, t_max=NULL){
      nlin <- function(x) log(x+spacing+1e-20)
      t = seq(0,1,length=lag)
      n = length(t)
      nt = nlin(t)
      t0 = ifelse(is.null(t_min), yes = nt[1], no = nlin(t_min))
      tN = ifelse(is.null(t_max), yes = nt[n], no = nlin(t_max))

      db      = (tN - t0) / (num_bases-1)
      a      = pi/db
      phi    = a*seq(t0, tN, db)
      basis    = matrix(nr=n, nc=num_bases)

      for(j in 1:num_bases){
        tmp = a*nt - phi[j]
        idx = which(tmp > pi | tmp < -pi)
        tmp[idx] = pi
        basis[,j] = (cos(tmp) + 1 )*0.5
      }
      # out = list(basis=basis, t=t)
      out = basis
      return(out)
    }


####################################################################################

    ## Create indicator basis functions matrix
    make_indicators <- function(lag, width){

      intervals = ceiling(lag/width)
      out = matrix(0,nr = intervals*width, nc = intervals)
      for(i in 1:intervals){
        out[((i-1)*width+1):(i*width),i] = 1
      }
      out = out[1:lag,]
      return(out)

    }
