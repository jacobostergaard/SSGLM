#include <RcppArmadillo.h>
#include <stdio.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace arma;

// [[Rcpp::export]]
arma::vec log_logistic(arma::vec x, double a, double b){
  // log-logistic distribution function
  // https://en.wikipedia.org/wiki/Log-logistic_distribution
  vec out(x.n_elem);

  x = x/a; // scale x
  out = pow(x,-b);
  out = pow(1+out,-1);

  return out;
}
