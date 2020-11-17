#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace arma;


// [[Rcpp::export]]
arma::mat izhikevich(arma::mat V, arma::mat U, arma::mat I, arma::vec pars) {

  // Transpose to work with columns vectors in each simulation step
  V = V.t();
  U = U.t();
  I = I.t();

  // Get n and number of neurons from the dimensions of V
  int n_sim = V.n_cols;
  int n_neu = V.n_rows;

  // Parameters for simulation
  double dt = pars[0];
  double a  = pars[1];
  double b  = pars[2];
  double c  = pars[3];
  double d  = pars[4];

  // Store values for V and U
  arma::vec Vt = zeros(n_neu);
  arma::vec Ut = zeros(n_neu);
  arma::vec dV = zeros(n_neu);
  arma::vec dU = zeros(n_neu);
  // indices at which V_t > 30
  uvec idx;


  // Simulate Izhikevich neurons (starting at second index)
  for(int i=1; i<n_sim; i++){

    Vt = V.col(i-1);
    Ut = U.col(i-1);

    if( sum(Vt >= 30) > 0 ){
      // U first, otherwise it messes up!
        idx = find(Vt > 30);
        Ut.elem(idx) = Ut.elem(idx)+d;
        Vt.elem(idx) = ones(size(idx))*c;
    }

    dV = 0.04*pow(Vt,2)+5*Vt+140-Ut+I.col(i);
    dU = a*(b*Vt-Ut);

    V.col(i) = Vt+dV*dt;
    U.col(i) = Ut+dU*dt;

  }

  // Transpose back to input layout
  arma::mat out = join_cols(V.t(),U.t());
  // arma::mat out = V.t();

  return out;
}

