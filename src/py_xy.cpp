#include <RcppArmadillo.h>
#include <stdio.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace arma;

// Declare functions
//  Find Poisson probility given vector input values and scalar intensity
arma::vec dpois(arma::vec x, double lambda);


// [[Rcpp::export]]
arma::mat py_xy(arma::vec x, arma::vec y, double l0, double l1) {

  int nx = x.n_elem;
  int ny = y.n_elem;

  uvec idx; // generic index vector

  // Output matrix of size x-by-y
  mat out(ny,nx);

  // For each value of x, find the corresponding probabilities and insert in out matrix
  idx = find(  x <= 0 );  // find particles that are in resting mode...
  out.each_col(idx) = dpois(y,l0);   // copy Po(y;l1) to columns matched by idx

  idx = find(  x > 0 );  // find particles that are in bursting mode...
  out.each_col(idx) = dpois(y,l1);   // copy Po(y;l0) to columns matched by idx

  return out;
}

