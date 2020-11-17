#include <RcppArmadillo.h>
#include <stdio.h>
#include <Rcpp.h>
#include <ctime>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;


// [[Rcpp::export]]
arma::vec rev_vec(arma::vec x){ // reverse the input vector x

  int n = x.n_elem;
  uvec idx = conv_to<uvec>::from(cumsum(ones(n))-1);
  idx = sort(idx,"descend");
  x = x(idx);

  // vec out(n);
  // for(int i=0;i<n;i++){
  //   out(i) = x(n-i-1);
  // }
  // return out;
  return x;
}

//  Find Poisson probility given vector input values and scalar intensity
arma::vec dpois(arma::vec x, double lambda){
  // length of input vector x
  int n = x.n_elem;
  vec lout(n);
  vec y;

  for(int i=0;i<n;i++){
    // intermediary vector of integers: 1,...,x
    y = cumsum(ones<vec>(x(i)));
    // log-Poisson probabilities
    lout(i) =-lambda+x(i)*log(lambda)-sum(log(y));
  }

  // Return probabilites (not log probabilites)
  return(exp(lout));
}

// arma::vec dpois2(int x, arma::vec lambda){
//
//   // length of input vector x
//   int n = lambda.n_elem;
//   vec lout(n);
//   vec y;
//
//   for(int i=0;i<n;i++){
//     // intermediary vector of integers: 1,...,x
//     y = cumsum(ones<vec>(x));
//     // log-Poisson probabilities
//     lout(i) =-lambda(i)+x*log(lambda(i))-sum(log(y));
//   }
//
//   // Return probabilites (not log probabilites)
//   return(exp(lout));
// }


double ess(arma::vec w){ // Effective sample size sum of: 1/w^2 for each weight w
  return 1/as_scalar(sum(pow(w,2)));
}

// create a vector v=(0,1,2,...,n)'
arma::uvec idx_vec(int n){
  // return linspace<uvec>(0, n-1);
  return conv_to<uvec>::from(cumsum(ones(n))-1);
}


// print progress in percent
void display_progress(double pct){

  // change to percent
  pct = round(pct*100);

  if(pct==100){
    cout << "\rProgress: " << pct << "% completed";
  } else if(pct >= 10){
    cout << "\rProgress:  " << pct << "% completed";
  } else{
    cout << "\rProgress:   " << pct << "% completed";
  }

}



int rpois(double lam){ // generate poisson random number from uniform sampling
  // From Knuths algorithm

  double  L = exp(-lam);
  double  p = 1;
  int     k = 0;
  int     res;

  if(lam < 1e-10){ // lambda is practically zero, hence
    res = 0;
  } else{

    while(p>L){
      k = k+1;
      p = p*as_scalar(randu<vec>(1));
      res = k-1;
    }
  }


  return res;
}

