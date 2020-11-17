#include <RcppArmadillo.h>
#include <stdio.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;

// Resampling (size) of input vector x according to weights w

// [[Rcpp::export]]
arma::uvec resample(arma::vec w){

  int n = w.n_elem;

  if(sum(w)==0){ // if weight are all zero set them equal to 1/n, where n is the number of elements in x
    w = ones<vec>(n)/n;
  } else{ // normalize weights
    w = w/sum(w);
  }
  w = cumsum(w);

  uvec ind = cumsum(ones<uvec>(n))-1;
  vec rnd = randu(n);
  unsigned int idx;
  uvec out(n);

  for(int i=0;i<n;i++){
      idx = min(find(rnd(i)<w));

      // cout << find(rnd(i)<w).t() << endl;
      out(i) = ind(idx);
  }

  return out;
}
