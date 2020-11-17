#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace arma;

// Declare functions
arma::vec px_x(arma::vec x, arma::vec pars);        // Transition probability
int rpois(double lam);                              // generate poisson random number from uniform sampling
arma::vec rev_vec(arma::vec x);                     // aux. function: reverses the input vector: x1,...,xn => xn,...,x1
arma::vec sample_x(arma::vec x, arma::vec pars);    // function for sampling latent state x


// [[Rcpp::export]]
arma::mat sim_ssglm(int nsim, double r0, arma::vec rk, double b0, arma::vec bk, arma::vec xpars, double dt){
  /*  Input pars:
   *
   *    nsim    number of simulation steps
   *    r0      resting baseline
   *    rk      resting kernel
   *    b0      bursting baseline
   *    bk      bursting kernel
   *    xpars   transition parameters on the form (alpha_burst, beta_burst, alpha_rest, beta_rest)
   *    dt      timestep
   *
  */

  mat out(nsim,3);

  int lag_r = rk.n_elem;      // length of resting kernel (history dependence)
  int lag_b = bk.n_elem;      // length of bursting kernel (history dependence)
  vec x = -ones<vec>(nsim);   // latent state process
  vec y = zeros<vec>(nsim);   // observed events (spikes)
  vec lam = zeros<vec>(nsim); // intensity process
  int lag = max(lag_r,lag_b); // starting lag value
  // vec U = randu<vec>(nsim);   // uniform random to draw latent states
  vec    y_hst;               // placeholder for history at time t
  double tmp;                 // tmp scalar place holder
  double prb_x;               // placeholder for X transition probabilty
  vec tmpvec;

  // simulation loop
    for(int i=lag; i<nsim; i++){

        tmpvec = x[i-1];
        x[i] = as_scalar(sample_x(tmpvec, xpars));
        // prb_x = as_scalar(px_x(  x.row(i-1) , xpars));

        // if( U[i] > prb_x ){ // stay in mode
        //   x[i] = sign(x[i-1])*(abs(x[i-1])+1); // increase/decrease x by +/- 1
        // } else{ // switch mode
        //   x[i] = -sign(x[i-1]);
        // }

        if(x[i] <= 0){ // in resting mode => resting pars
          y_hst  = y.subvec(i-lag_r,i-1);
          tmp    = as_scalar(rk.t()*rev_vec(y_hst) + r0);
          lam[i] = exp(tmp);

        } else{ // in bursting mode => bursting pars
          y_hst  = y.subvec(i-lag_b,i-1);
          tmp    = as_scalar(bk.t()*rev_vec(y_hst) + b0);
          lam[i] = exp(tmp);

        }
        y[i] = rpois(lam[i]*dt);
    }

    out.col(0) = y;
    out.col(1) = x;
    out.col(2) = lam;

  return out;
}
