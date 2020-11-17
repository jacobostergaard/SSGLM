// /* This is an implementation of a marginal particle filter for a bursting neuron with a possible history dependent structure.
// * The model takes an observed spike train as input, as well as intial values of parameters. The iteration uses a particle filter to estimate
// * the latent state, determining whether the neuron is in a bursting regime or not. Parameters controlling the intensity process are baseline+history structure,
// * where the latter is constructed using spline basis functions, hence the spline parameters (weights) determine the history dependent structure.
// *
// * Model equations:
// *
// *  Observed process  Y[n] ~ Pois(lambda[n]) where lambda[n] = 1{X[n]<=0}(b0+k0[n])+1{X[n]>0}(b1+k1[n])
// *                    where b0, b1 are the baselines for resting/bursting respectively and
// *                    k0,k1 are history kernels where each is determined by the spline basis functions and the past history of Y[n]:
// *                    k[n] = B'*(Y[n-1],...,Y[n-l])', where l is the max lag that is included in the history and B is a l-by-1 vector
// *                    B=S*w, where S is the l-by-p design matrix with spline basis functions and w is a p-by-1 vector of spline weights.
// *
// *
// *  Latent process    X[n] ~ P(X[n] |X[n-1],Y[n-1]) depends only on Y[n-1] if X[n-1]=0. In this case
// *                    P(X[n]=0|X[n-1]=0,Y[n-1]=0) = 0 = 1-P(X[n]=1|X[n-1]=0,Y[n-1]=1), otherwise
// *                    P(X[n]=X[n-1]+1|X[n-1]>0) = 1-(X[n-1]/k)^alpha
// *
// * August 2017, Jacob Østergaard
// */
//
// #include <RcppArmadillo.h>
// // #include <Rcpp.h>
// #include <stdio.h>
// #include <math.h>
//
// // [[Rcpp::depends(RcppArmadillo)]]
// using namespace std;
// using namespace arma;
// using namespace Rcpp;
//
// // Declare functions
// arma::vec px_x(arma::vec x, arma::vec pars);                              // Transition probability
// arma::mat py_xy(arma::vec x, arma::vec y, double l0, double l1);          // likelihood of y, given latent state X and parameters l0, l1
// // arma::vec dpois(arma::vec x, double lambda);                              // Poisson probility given vector input values and scalar intensity
// // arma::vec dpois2(double x, arma::vec lambda);                              // Poisson probility given vector input values and scalar intensity
// arma::uvec resample(arma::vec w);                                         // resample indices from input weights w
// arma::uvec idx_vec(int n);                                                // create a vector of indices 1,...,n which in c++ is 0,...,n-1
// arma::vec rev_vec(arma::vec x);                                           // aux. function: reverses the input vector: x1,...,xn => xn,...,x1
// arma::vec sample_x(arma::vec x, arma::vec pars);                         // function for sampling particles
// double ess(arma::vec w);                                                  // aux. function: effective sample size
// void display_progress(double pct);                                        // aux. function: will display the progress of the main loop over observations
//
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
//
//
//
// // [[Rcpp::export]]
// arma::mat mpfpp_cpp(arma::sp_mat y, arma::vec x, arma::vec Bb, arma::vec Br, arma::mat Sb, arma::mat Sr, arma::vec xpars, int s_lag_in, double Neffpct, int M, double w0, arma::vec q, double dt, bool usetrue, bool verbose) {
//
//
//   /* This is the main function, it returns a matrix with N (number of observations) columns, with parameter estimates, covariance matrices (vectorized), particles and estimated intensities
//   *
//   *  Inputs
//   *    y:         spike train (0-1 observations)
//   *    x:         latent states
//   *    Bb:        inital values of bursting parameters: baseline and kernel pars, if input is exactly zero this is interpreted as leaving out this parameter
//   *    Sb:        bursting spline basis functions (a matrix) if some columns are zeros, then effectively the parameter is removed (such as when resting kernel is 0)
//   *    Br:        inital values of resting parameters: baseline and kernel pars, if input is exactly zero this is interpreted as leaving out this parameter
//   *    Sr:        resting spline basis functions (a matrix) if some columns are zeros, then effectively the parameter is removed (such as when resting kernel is 0)
//   *               note that is the number of (non-zero) columns in S differ from the number of entries in B-2, then it is assumed that a baseline value is present
//   *    xpars:     Parameters for the distributions determining latent state dynamics, format: (alpha1,beta1,alpha2,beta2)
//   *    s_lag:     Fixed lag smoother
//   *    w0:        inital value of variance
//   *    q:         variances for RW propagation model
//   *    a1,a2:     initial values of alpha's for log-logistic funtions
//   *    b1,b2:     initial values of beta's for log-logistic funtions
//   *    M:         number of particles
//   *    dt:        time step
//   *    usetrue:   boolean that determines whether parameters should be estimated (otherwise the initial guess is used throughout)
//   *    verbose:   boolean whether to show the progress status while running the main loop (over observations)
//   */
//
//   /*  Timeline:
//    *
//    *    2018-10-30: Restarted package structure with the intend to make it publicly available.
//    *    2018-11-06: Currently the mpfpp only works with a bursting kernel: hence no resting kernel is included and only one baseline parameter.
//    *    2018-12-03: Working on time scale separation for latent and observed states...
//    *    2018-12-14: Need to use sparse matrices (vectors) to deal with memory overload issues...
//    *    2019-04-25: Introduced non-zero resting baseline paramter. TODO: non-zero resting kernel...
//    *    2020-07-01: Updated and restructured code to simplify readability
//    */
//
//     // if X trajectory is given, no need to filter => no need for the fixed lag smoother
//
//
//         int N    = y.n_elem; // number of observations
//
//         if(x.n_elem == N){
//           M = 1;
//         }
//
//         if(y.n_cols > 1){
//           y = y.t();        // make sure y is interpreted as a column vector.
//         }
//
//         if(x.n_cols > 1){
//           x = x.t();        // make sure y is interpreted as a column vector.
//         }
//
//     // Number of pars equal number of columns in spline matrix + baseline(s), i.e. p = Sb.n_cols + 2
//         int pb = Bb.n_elem;   // no. of bursting parameters
//         int pr = Br.n_elem;   // no. of resting parameters
//         int p  = pb+pr;       // total no. of parameters
//         int p2 = p*p;         // number of parameters in covariance matrices
//
//     // Declaring overall parameters and containers
//         int lag  = Sb.n_rows;                           // length of the bursting kernel
//         int s_lag;                                      // s_lag used at time n, only needed when n < s_lag_in (everything is replaced then!)
//         // int lagr  = Sr.n_rows;                          // length of the resting kernel
//         // int lag   = max(Sb.n_rows,Sr.n_rows);           // overall lag length
//
//
//         mat  X  = ones<mat>(M,N);                       // container for M particles through N timesteps
//         mat  Y  = zeros<mat>(M,N);                      // container for M intensities through N timesteps
//         cube P  = zeros<cube>(M,N,p);                   // container for p parameters, for M particles and N timesteps
//         cube W  = zeros<cube>(M,N,p2);                  // container for p^2 covariance parameters, for M particles and N timesteps
//
//         vec  Xn = zeros<vec>(M);                        // working particles (size M)
//         vec  Yn = zeros<vec>(M);                        // working intensities (size M)
//         mat  Pn = zeros<mat>(M,p);                      // working parameters (p pars, M particles)
//         cube Wn = zeros<cube>(M,p,p);                   // working covariance parameters (p^2 pars, M particles)
//
//         rowvec Xout = zeros<rowvec>(N);                 // aggregated output of particles at each timestep
//         rowvec Yout = zeros<rowvec>(N);                 // aggregated output of intensities  at each timestep
//         rowvec Zout = zeros<rowvec>(N);                 // probability of bursting at each timestep
//         mat Pout    = zeros<mat>(p,N);                  // aggregated parameters at each timestep
//         mat Wout    = zeros<mat>(p2,N);                 // aggregated covariance parameters at each timestep
//
//         mat w    = ones<mat>(M,N)/M;                    // vector of previous step particle weights
//         vec wn  = zeros<vec>(M);                        // working weights
//
//         uvec idx(M);                        // index vector used for resampling
//         uvec idx_h;                         // index vector used for current history steps
//         double Neff;                        // efficient sample size
//         vec y_hst      = zeros<vec>(lag);   // trailing past observations used in weighing and estimation
//         mat dloglam    = zeros<mat>(M,p);   // container of derivatives for observations update step
//
//         double tmp;                         // generic scalar
//         vec tmpvec;                         // generic vector
//         vec tmpvec2;                        // generic vector
//         mat tmpmat;                         // generic matrix
//         mat tmpmat2;                        // generic matrix
//         mat tmpmat3;                        // generic matrix
//         cube tmpcube;                       // generic cube
//
//
//         if(x.n_elem == N){
//           // use true values of X, if given.
//           X.row(0) = x.t();
//         } else if(x.n_elem == M){
//           X.col(0) = x;
//         } else if(x.n_elem == 1){
//           X = x(0)*X;
//         }
//
//
//         // if(x.n_elem == M){
//         //   tmpvec = ones<vec>(M);
//         //   idx = find(x<0);
//         //   tmpvec(idx) = -ones<vec>(idx.n_elem);
//         //   for(int i=1; i<lag+1; i++){
//         //     X.col(i) = X.col(i)+tmpvec;
//         //   }
//         // } else if(x.n_elem==1){
//           for(int i=1; i<lag+1; i++){
//             X.col(i) = X.col(i)+x(0);
//           }
//         // }
//
//         if(x.n_elem!=N){
//           for(int n=1; n<lag+1; n++){
//             X.col(n) = X.col(0)*n;
//           }
//         }
//         for(int n=0; n<lag+1; n++){
//           Xout(n) = sum(X.col(n))/M;
//         }
//
//     // Prediction model for parameters is P[k+1] = F*P[k]+e, where e~N(0,Q)
//         mat F = eye<mat>(p,p);      // F matrix is the (p x p) identity matrix
//         mat Qmat = eye<mat>(p,p);   // Q matrix is a scaled (p x p) identity matrix
//
//
//           if(q.n_elem>1){
//             // if predicition model variances are given for individual parameters, set diagonal
//             Qmat.diag() = q;
//           } else{
//             Qmat =q(0)*Qmat;  // Q matrix is diagonal with equal entries, scaled by input q
//           }
//           tmpvec = vectorise(Qmat);         // vectorize to collapse dimensions
//           mat Q = repmat(tmpvec.t(),M,1);   // used in time update of variance
//
//
//     // define initial values
//         vec P0 = join_cols(Bb,Br);         // initial kernel/baseline weights
//         mat W0 = w0*eye<mat>(p,p)+w0/100;  // initial variance matrix and ensure that covariances aren't initialized at zero...
//
//         for(int i=0; i<p; i++){
//           P.slice(i) = P0(i)*ones<mat>(M,N);      // all parameters are set to initial values
//         }
//
//         tmpvec= vectorise(W0);
//         for(int i=0; i<p2; i++){
//           W.slice(i) = tmpvec(i)*ones<mat>(M,N);  // all covariance parameters are set to initial values
//         }
//
//
//     /*                    END OF SETUP                      */
//
//
//
//
//
//
//
// // *********************************************    main loop over observations    *********************************************
//
//     for(int n=lag;n<N;n++){ // loop over N observations, omitting the first 'lag' observations and 's_lag' (due to fixed lag smoothing)
//
//             if(n < s_lag_in){
//               s_lag = n;
//             } else{
//               s_lag = s_lag_in;
//             }
//
//             // 1: sample particles according to the px_xy distribution
//
//                   if(x.n_elem!=N){
//                     // sample particles
//                     X.col(n) = sample_x(X.col(n-1), xpars);
//                   }
//
//
//             // 2: predict parameters and set intensities
//
//                   tmpmat = P.col(n-1);
//                   P.col(n) = tmpmat*F;    // time prediction of parameters
//
//                   tmpmat   = W.col(n-1);
//                   W.col(n) = tmpmat + Q;  // time prediction of covariances
//
//
//             // 3: set intensities and weights by likelihood function
//
//                       Pn = P.col(n);      // set working parameters
//                       Xn = X.col(n);      // set working covariance parameters
//
//                       if(pb>1){
//                         // history dependency in the bursting kernel, ie spline based regression, otherwise only baseline is included
//                         y_hst  = y.rows(n-lag,n-1);
//
//                         tmpmat = Pn.cols(1,pb-1)*Sb.t();
//                         tmpvec = tmpmat*rev_vec(y_hst)+Pn.col(0);
//
//                       } else{
//                         // only baseline firing is included when bursting, hence no contribution other than baseline intensity
//                         tmpvec = Pn.col(0);
//                       }
//
//                       Yn = exp(tmpvec);  // set current intensity using bursting modulation
//
//                       if(any(Xn<=0)){ // detect resting particles
//                         idx    = find(Xn <= 0);        // find particles in rest mode
//                         tmpvec = Pn.col(pb);   // extract resting parameters
//                         Yn(idx) = exp(tmpvec(idx));  // set resting intensity level
//                       }
//
//                       Y.col(n) = Yn;
//
//
//                       tmpvec = dpois2(y(n), Yn*dt);
//                       tmpvec = log(tmpvec)+log(w.col(n-1));
//                       tmpvec = exp(tmpvec);
//
//                   // normalize weights
//                       if(sum(tmpvec)>0){
//                         tmpvec = tmpvec/sum(tmpvec);   // normalize weights if not all are not zero...
//                       } else{
//                         tmpvec = ones<vec>(M)/M; // ... else set all weights equal to 1/M
//                       }
//
//                       wn = tmpvec;
//
//
//
//
//             // 4: resample particles and reset weights
//
//                       Neff = as_scalar(pow(sum(pow(wn,2)),-1));
//
//
//
//                       if(Neff < Neffpct*M){
//                         // Particles are resampled, should technically be from 0 to n... but it slows the runtime significantly down!
//                         // Hence, resample up to s_lag lagged values, i.e. a fixed lag smoother.
//
//                           idx   = resample(wn);        // resample indices based on weights
//                           idx_h = regspace<uvec>(n-s_lag, n);
//
//
//                           // tmpmat  = X.cols(n-s_lag,n);
//                           // X.cols(n-s_lag,n) = tmpmat.rows(idx);
//                           X.cols(idx_h) = X.submat(idx,idx_h);
//                           Y.cols(idx_h) = Y.submat(idx,idx_h);
//
//                           tmpcube = P.cols(n-s_lag,n);
//                           for(int j=0; j<p; j++){
//                             tmpmat           = tmpcube.slice(j);
//                             tmpcube.slice(j) = tmpmat.rows(idx);
//                           }
//                           P.cols(n-s_lag,n) = tmpcube;
//
//                           tmpcube = W.cols(n-s_lag,n);
//                           for(int j=0; j<p2; j++){
//                             tmpmat           = tmpcube.slice(j);
//                             tmpcube.slice(j) = tmpmat.rows(idx);
//                           }
//                           W.cols(n-s_lag,n) = tmpcube;
//
//                           // reset weights
//                           wn = ones<vec>(M)/M;
//                       }
//
//                       w.col(n) = wn;
//
//                       Xn = X.col(n);
//                       Yn = Y.col(n);
//                       Pn = P.col(n);
//
//                       for(int m=0; m<M; m++){
//                         tmpmat = W.tube(m,n);
//                         tmpmat.reshape(p,p);
//                         Wn.row(m) = tmpmat;
//                       }
//
//
//
//             // 5: part i) evaluate particles and set derivative of log-intensity
//
//                       // reset dloglam (derivative of log-intensity) before filling in the blanks
//                       dloglam = 0*dloglam;
//
//                       if( any(Xn>0) ){
//                         idx = find(Xn > 0);
//                         tmpmat = dloglam.rows(idx);
//                         tmpvec = zeros<vec>(p);
//                         tmpvec(0) = 1;
//                         tmpvec.subvec(1,pb-1) = Sb.t()*rev_vec(y_hst);
//                         tmpmat.each_row() = tmpvec.t();
//                         dloglam.rows(idx) = tmpmat;
//                       }
//
//                       if( any(Xn<=0) ){
//                         idx = find(Xn <= 0);
//                         tmpmat = dloglam.rows(idx);
//                         tmpvec = zeros<vec>(p);
//                         tmpvec(pb) = 1;
//                         tmpmat.each_row() = tmpvec.t();
//                         dloglam.rows(idx) = tmpmat;
//                       }
//
//
//
//                       // set derivative of loglam to 0 for fixed parameters => no update on information
//                       if(any(Qmat.diag()==0)){
//                         idx = find(Qmat.diag() == 0);
//                         dloglam.cols(idx) = 0*dloglam.cols(idx);
//                       }
//
//
//             // 5: part ii) update linear posteriors for each particle based on current observation
//
//                       // find weighted estimates of parameters
//
//                       for(int m=0;m<M;m++){
//                             // update variance estimate for particle m
//                             tmpmat    = dloglam.row(m).t()*dloglam.row(m)*Yn(m)*dt; // the second derivative part is zero, hence only part using first derivatives is needed
//                             tmpmat2   = Wn.row(m);
//                             tmpmat    = inv( tmpmat2.i()+tmpmat);
//                             Wn.row(m) = tmpmat;
//                             W.tube(m,n) = vectorise(tmpmat);
//
//                             // update parameter estimate for particle m
//                             tmpvec    = dloglam.row(m).t()*( y(n)-Yn(m)*dt);
//                             tmpvec    = tmpmat*tmpvec;
//                             if(any(Qmat.diag()==0)){
//                               idx = find(Qmat.diag() == 0);
//                               tmpvec(idx) = 0*tmpvec(idx);
//                             }
//
//                             Pn.row(m) = Pn.row(m) + tmpvec.t();
//                       }
//
//                       if(usetrue){
//                         Pn.each_row() = P0.t();
//                       }
//                       P.col(n) = Pn;
//
//                     if(verbose){ // display current progress
//                       display_progress((double) (n-lag)/(N-lag));
//                     }
//       } // end of observation loop
//
//
//     // aggregate parameter estimates
//     for(int n=lag;n<N;n++){ // loop over N observations, omitting the first 'lag' observations and 's_lag' (due to fixed lag smoothing)
//
//           wn = w.col(n);
//           Pn = P.col(n);
//
//           tmpvec = Pn.t()*wn;
//           Pout.col(n) = tmpvec;
//
//
//           tmpmat = Pn-repmat(tmpvec.t(),M,1);
//
//           Wn = W.col(n);
//
//           tmpmat3 = zeros<mat>(p,p);
//
//           for(int m=0;m<M;m++){
//                 tmpmat2  = Wn.row(m);
//                 tmpmat2.reshape(p,p);
//                 tmpmat3 += wn(m)*( tmpmat2+tmpmat.row(m).t()*tmpmat.row(m) );
//               }
//           Wout.col(n) = vectorise(tmpmat3);
//
//           Xn  = X.col(n);
//           Yn  = Y.col(n);
//
//           Xout(n) = as_scalar(Xn.t()*wn);
//           Yout(n) = as_scalar(Yn.t()*wn);
//
//           idx = find(Xn>0);             // bursting states
//           // tmp = idx.n_elem/M;           // proportion of bursting particles
//           tmp = sum( wn(idx) )/sum(wn); // estimated probability of bursting in current state by weights of bursting particles
//
//           Zout(n) = tmp;
//     }
//
//
//     cout << "\n" << endl;
//
//     mat out;
//     out = join_cols(Yout,Xout);
//     out = join_cols(out,Zout);
//     out = join_cols(out,Pout);
//     out = join_cols(out,Wout);
//
//     return out;
// }
//
