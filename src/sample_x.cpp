// Sample particles for  marginalized particle filter

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace arma;

// Declare functions
arma::vec px_x(arma::vec x, arma::vec pars);      // Transition probability for x state

// [[Rcpp::export]]
arma::vec sample_x(arma::vec x, arma::vec pars) {

  /*
   * Inputs:
   * x:     vector of particles to sample
   * pars:  parameters for the transition probabilites
   */

  vec tmpvec;                   // generic vector
  uvec idx;                     // generic index vector
  int M = x.n_elem;             // number of particles
  vec randvec = randu<vec>(M);  // for sampling particles
// cout << randvec << endl;
  vec x_prob = px_x(x,pars);    // probability of switching given current states
  vec signvec = sign(x);        // determines current mode: resting=-1, bursting=1

  if(any(signvec==0)){          // exception handling of any x=0... (should be impossible) then treat as sign(x)=-1
    idx = find(signvec == 0);
    signvec(idx) = -ones<vec>(idx.n_elem);
  }

  idx = find( randvec < x_prob == 1 );  // find particles that switch mode...
  x(idx)   = -signvec(idx);             // and switch
  idx = find( randvec < x_prob == 0 );  // find particles that stay in same mode...
  x(idx) = x(idx)+signvec(idx);         // and increase time in that mode

  return x;
}
