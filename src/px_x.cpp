// Transition probability for x[k] given x[k-1]

#include <RcppArmadillo.h>
#include <stdio.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace arma;


// Declare functions
arma::vec log_logistic(arma::vec x, double a, double b);


// [[Rcpp::export]]
arma::vec px_x(arma::vec x, arma::vec pars) {

    /*  x:    input vector of states
     *  pars: input vector of parameters
     */

    // Output vector of length equal to input x

    int M = x.n_elem;           // number of particles
    uvec idx;                   // generic index vector
    vec out(M);
    vec tmp(1);

    double alpha1 = pars(0);
    double beta1  = pars(1);
    double alpha2 = pars(2);
    double beta2  = pars(3);

    // Output probabilites based on log-logistic function, varies depending on the sign of x
    idx = find(  x > 0 );  // find particles that switch mode...
    out(idx) = log_logistic(x(idx),alpha1,beta1);

    idx = find(  x <= 0 );  // find particles that switch mode...
    out(idx) = log_logistic(-x(idx),alpha2,beta2);

    // for(int i=0;i<x.n_elem;i++){
    //   tmp = x(i);
    //
    //   if(x(i)>0){ // x is positive (bursting) hence use parameters a1,b1
    //     out(i) = as_scalar(log_logistic(tmp,a1,b1));
    //   } else{     // x is negativ (resting) hence use parameters a2,b2
    //     out(i) = as_scalar(log_logistic(-tmp,a2,b2));
    //   }
    // }
    return out;
}

