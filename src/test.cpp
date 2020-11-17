#include <RcppArmadillo.h>
#include <stdio.h>
#include <Rcpp.h>
#include <ctime>
// #include <conio.h>
// #include <iostream.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;

// Declare functions
// resample indices from input weights w
arma::uvec resample(arma::vec w);


void display_progress(double pct);



// [[Rcpp::export]]
arma::cube test(int p){ // generate random index

  field<cube> F(p,1);

  cube A(2,2,3);

  for(int i=0;i<p;i++){
    F(i,0) = A;
    cout << F(i,0) << endl;
  }

  cout << F.n_elem << endl;
  cout << F.n_rows << endl;
  cout << F.n_cols << endl;
  cout << F.n_slices << endl;


  return A;
}


// [[Rcpp::export]]
arma::mat test_input(arma::mat y){ // look at sparse input

  vec x = vectorise(y);

  int nr = y.n_rows;
  int nc = y.n_cols;

  mat z = reshape( x, nr, nc );


  cout << y << endl;
  cout << x << endl;
  cout << z << endl;


  return z;
}

