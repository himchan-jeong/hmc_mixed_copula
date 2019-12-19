# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// In all functions, H is (converted to) a cube, with 
// dim 0: sample no, dim 1: phi from likelihood, dim 2: phi from sample
// Thus, H[[j]][i,k] = H[i,k,j]

//Convert RCpp::NumericVector object, which comes from R 3d array, to arma::cube
arma::cube vect2cube(Rcpp::NumericVector H_){
  Rcpp::IntegerVector x_dims = H_.attr("dim");
  arma::cube H(H_.begin(), x_dims[0], x_dims[1], x_dims[2], false);
  return(H);
}

// [[Rcpp::export]]
double log_sum_exp(Rcpp::NumericVector x){
  // compute log(exp(x[0]) + exp(x[1]) + ... )
  double m = max(x);
  int n = x.size();
  double s = 0.0;
  for(int i=0; i<n; ++i){
    s += exp(x[i] - m);
  }
  return(m + log(s));
}

// // Copula double derivative
// // [[Rcpp::export]]
// double dCphi_du1du2(double u1, double u2, double phi){
//   double b;
//   b = pow(((1-exp(-phi)) - (1-exp(-phi*u1))*(1-exp(-phi*u2))),-2);
//   return(phi*(1-exp(-phi))*exp(-phi*(u1+u2))*b);
// }

// Copula double derivative, with option for different copulas
// 1=normal, 2=clayton, 3=gumbel, 4=frank
// [[Rcpp::export]]
double dCphi_du1du2(double u1, double u2, double phi, int copula) {
  double c;
  // temporary 1d vectors, hack to use qnorm of Rcpp
  Rcpp::NumericVector a1(1);Rcpp::NumericVector a2(1);
  a1[0] = u1; a2[0] = u2;
  if (copula==4) { // Frank
    c = phi * 
        (1-exp(-phi)) * 
        exp(-phi*(u1+u2)) * 
        pow((1-exp(-phi)) - (1-exp(-phi*u1))*(1-exp(-phi*u2)), -2);
  } else if (copula==3) { //Gumbel
    if (phi < 1) {
      throw std::invalid_argument("phi must be >= 1 for Gumbel copula");
    }
    c = exp(-pow((pow(-log(u1), phi) + pow(-log(u2), phi)), 1/phi)) * 
        (1/(u1*u2)) * 
        pow(pow(-log(u1), phi) + pow(-log(u2), phi), (-2 + 2/phi)) * 
        pow(log(u1)*log(u2), phi-1) * 
        (1 + (phi-1)*pow(pow(-log(u1), phi) + pow(-log(u2), phi), -1/phi));
  } else if (copula==2) { //Clayton
    if (phi <= 0) {
      throw std::invalid_argument("phi must be > 0 for Clayton copula");
    }
    c = (1 + phi) * pow(u1*u2, -1 - phi) / 
        pow(pow(u1, -phi) + pow(u2, -phi) - 1, 2 + 1/phi);
  } else if (copula==1) { //Gaussian
    // phi = rho, the correlation
    if (phi <= -1 || phi >= 1) {
      throw std::invalid_argument("phi must be between -1 and 1 for Gaussian copula");
    }
    c = exp(-(pow(phi, 2)*(pow(qnorm(a1, 0.0, 1.0, false)[0], 2) + 
                           pow(qnorm(a2, 0.0, 1.0, false)[0], 2)) - 
                2*phi*qnorm(a1, 0.0, 1.0, false)[0]*qnorm(a2, 0.0, 1.0, false)[0]) / 
            (2*(1 - pow(phi, 2)))
            ) / sqrt(1 - pow(phi, 2));
  } else {
    throw std::invalid_argument("copula parameter must be an integer between 1 and 4.");
  }
  return c;
}

// Copula single derivative (h function), with option for different copulas
// 1=normal, 2=clayton, 3=gumbel, 4=frank
// [[Rcpp::export]]
double dCphi_du2(double u1, double u2, double phi, int copula) {
  double c;
  // temporary 1d vectors, hack to use qnorm of Rcpp
  // Rcpp::NumericVector a1(1);Rcpp::NumericVector a2(1);
  // a1[0] = u1; a2[0] = u2;
  if (copula==4) { // Frank
    c = -exp(phi) * (exp(phi*u1) - 1) / 
      (exp(phi*u2 + phi*u1) - exp(phi*u2 + phi) - exp(phi*u1 + phi) + exp(phi));
  } else if (copula==3) { //Gumbel
    if (phi < 1) {
      throw std::invalid_argument("phi must be >= 1 for Gumbel copula");
    }
    c = -exp(-pow((pow(-log(u1), phi) + pow(-log(u2), phi)), 1/phi)) * 
      pow(pow(-log(u1), phi) + pow(-log(u2), phi), 1/phi - 1) *
      pow(-log(u2), phi) /
        (u2 * log(u2));
  } else if (copula==2) { //Clayton
    if (phi <= 0) {
      throw std::invalid_argument("phi must be > 0 for Clayton copula");
    }
    c = pow(u2, -phi-1) * pow(pow(u1, -phi) + pow(u2, -phi) - 1, -1 -1/phi);
  } else if (copula==1) { //Gaussian
    // phi = rho, the correlation
    if (phi <= -1 || phi >= 1) {
      throw std::invalid_argument("phi must be between -1 and 1 for Gaussian copula");
    }
    c = Rcpp::stats::pnorm_0((Rcpp::stats::qnorm_0(u1, 1, 0) - phi*Rcpp::stats::qnorm_0(u2, 1, 0)) / sqrt(1 - pow(phi, 2)), 1, 0);
  } else {
    throw std::invalid_argument("copula parameter must be an integer between 1 and 4.");
  }
  return c;
}

// /////////////////////////////////////////////////////////////////////////////////////////
// ///                                  LATENT FORMULATION FUNCTIONS                    ////
// /////////////////////////////////////////////////////////////////////////////////////////
// 
// // Copula pdf, latent structure. Marginals excluded as they cancel out in eta/BF functions
// // [[Rcpp::export]]
// double copulapdf_log(Rcpp::NumericVector y1, Rcpp::NumericVector z,
//                         Rcpp::NumericVector mu1, Rcpp::NumericVector mu2,
//                         double sigma1, double phi, int copula){
//   double logl = 0.0;
//   double s;
//   int n = y1.size();
//   // temporary 1d vectors, hack to use d/p functions of RCpp
//   Rcpp::NumericVector a1(1);Rcpp::NumericVector a2(1);
//   for(int i=0; i<n; ++i){
//     a1[0] = y1[i]; a2[0] = z[i];
//     s = log(dCphi_du1du2(pnorm(a1,mu1[i],sigma1,false)[0], 
//                          plogis(a2,mu2[i],1.0,false)[0], 
//                          phi,
//                          copula));
//     logl += s;
//   }
//   return(logl);
// }
// 
// // [[Rcpp::export]]
// double copulapdf_wrapcov_log(Rcpp::NumericVector y1, Rcpp::NumericVector z, 
//                                arma::mat X1, arma::mat X2,
//                                arma::mat b1, arma::mat b2, double sigma1, 
//                                double phi, int copula){
//   arma::mat t1 = X1 * b1;
//   arma::mat t2 = X2 * b2;
//   Rcpp::NumericVector m1 = Rcpp::NumericVector(t1.begin(), t1.end());
//   Rcpp::NumericVector m2 = Rcpp::NumericVector(t2.begin(), t2.end());
//   return(copulapdf_log(y1, z, m1, m2, sigma1, phi, copula));
// }
// 
// // log(p_j(x, eta)) from Geyer
// // j will be on 0-indexing (add 1 to compare to paper formulas)
// // [[Rcpp::export]]
// double p_j_log(Rcpp::NumericVector eta, Rcpp::NumericVector skelpts, int j, 
//                Rcpp::NumericVector y1, Rcpp::NumericVector z,
//                arma::mat b1, arma::mat b2, double sigma1,
//                arma::mat X1, arma::mat X2,
//                int copula) {
//   // create vector to feed into log_sum_exp;
//   int M = skelpts.size();
//   Rcpp::NumericVector lse(M);
//   double log_h_j = copulapdf_wrapcov_log(y1, z, X1, X2, b1, b2, sigma1, 
//                                          skelpts[j], copula);
//   for(int k=0; k<M; ++k){
//     lse[k] = copulapdf_wrapcov_log(y1, z, X1, X2, b1, b2, sigma1, 
//                                    skelpts[k], copula) - 
//       log_h_j + eta[k] - eta[j];
//   }
//   return(-log_sum_exp(lse));
// }
// 
// // function to maximize for estimating r
// // Here, the parameter matrices and z have dimensions (NM x p), where p is p1 or p2 or n etc.
// // They should have been supplied as row concatenated scalars/vectors from the MCMC samples
// // The structure of concatenation is N, N, N,...N (M times).
// // For maximization, provide 2nd element onward for eta. First is assumed zero wlog.
// // [[Rcpp::export]]
// double etafn(Rcpp::NumericVector eta_reduced, Rcpp::NumericVector skelpts,
//              Rcpp::NumericVector y1, Rcpp::NumericMatrix z,
//              arma::mat b1, arma::mat b2,
//              Rcpp::NumericVector sigma1,
//              arma::mat X1, arma::mat X2,
//              int copula){
//   // augment the eta. wlog first element of eta is zero.
//   int n_eta = eta_reduced.size();
//   Rcpp::NumericVector eta(n_eta+1); //defaults of zero
//   for(int i=0; i<n_eta; ++i) eta[i+1] = eta_reduced[i];
//   
//   int M = skelpts.size();
//   if(M != n_eta + 1) throw std::invalid_argument("In etafn(), length of eta_reduced should be number of skeleton points - 1");
//   int N = sigma1.size()/M;
//   double s = 0.0;
//   for(int j=0; j<M; ++j){
//     for(int i=0; i<N; ++i){
//       int ij = (N*j + i);
//       Rcpp::NumericVector z_ij = z(ij, Rcpp::_);
//       arma::mat b1_ij = b1.row(ij).t();
//       arma::mat b2_ij = b2.row(ij).t();
//       double sigma1_ij = sigma1[ij];
//       s += p_j_log(eta, skelpts, j, 
//                    y1, z_ij,
//                    b1_ij, b2_ij, sigma1_ij,
//                    X1, X2, 
//                    copula);
//     }
//   }
//   return(s);
// }
// 
// // gradient of etafn, return vector of length length(eta_reduced).
// // note that we omit the first element of the gradient since we assume wlog that eta[0] = constant (0).
// // [[Rcpp::export]]
// Rcpp::NumericVector etafn_gr(Rcpp::NumericVector eta_reduced, Rcpp::NumericVector skelpts,
//                              Rcpp::NumericVector y1, Rcpp::NumericMatrix z,
//                              arma::mat b1, arma::mat b2,
//                              Rcpp::NumericVector sigma1,
//                              arma::mat X1, arma::mat X2,
//                              int copula){
//   // augment the eta. wlog first element of eta is zero.
//   int n_eta = eta_reduced.size();
//   Rcpp::NumericVector eta(n_eta+1); //defaults of zero
//   for(int i=0; i<n_eta; ++i) eta[i+1] = eta_reduced[i];
// 
//   int M = skelpts.size();
//   if(M != n_eta + 1) throw std::invalid_argument("In etafn(), length of eta_reduced should be number of skeleton points - 1");
//   int N = sigma1.size()/M;
//   Rcpp::NumericVector gr(n_eta);
// 
//   for(int r=1; r<n_eta+1; ++r){
//     double s = 0.0;
//     for(int j=0; j<M; ++j){
//       for(int i=0; i<N; ++i){
//         int ij = (N*j + i);
//         Rcpp::NumericVector z_ij = z(ij, Rcpp::_);
//         arma::mat b1_ij = b1.row(ij).t();
//         arma::mat b2_ij = b2.row(ij).t();
//         double sigma1_ij = sigma1[ij];
//         s += exp(p_j_log(eta, skelpts, r,
//                      y1, z_ij,
//                      b1_ij, b2_ij, sigma1_ij,
//                      X1, X2,
//                      copula));
//       }
//     }
//     gr[r-1] = N - s;
//   }
//   return(gr);
// }
// 
// // Component in sum of Bayes factor estimator
// // Does not have N_k factor in lse as in Roy, because we assume N_k is constant
// // ie same number of iterations per skel pt, hence it can be taken outside and is
// // irrelevant for maximization of Bfn
// // [[Rcpp::export]]
// double p_phi_log(Rcpp::NumericVector eta, Rcpp::NumericVector skelpts, double phi, 
//                Rcpp::NumericVector y1, Rcpp::NumericVector z,
//                arma::mat b1, arma::mat b2, double sigma1,
//                arma::mat X1, arma::mat X2,
//                int copula){
//   // create vector to feed into log_sum_exp;
//   int M = skelpts.size();
//   Rcpp::NumericVector lse(M);
//   double log_h_phi = copulapdf_wrapcov_log(y1, z, X1, X2, b1, b2, sigma1, phi,
//                                            copula);
//   for(int k=0; k<M; ++k){
//     lse[k] = copulapdf_wrapcov_log(y1, z, X1, X2, b1, b2, sigma1, skelpts[k],
//                                    copula) - 
//       log_h_phi + eta[k] - eta[0];
//   }
//   return(-log_sum_exp(lse));
// }
// 
// 
// // Bayes factor estimator
// // [[Rcpp::export]]
// double Bfn(double phi, Rcpp::NumericVector eta, Rcpp::NumericVector skelpts,
//              Rcpp::NumericVector y1, Rcpp::NumericMatrix z,
//              arma::mat b1, arma::mat b2,
//              Rcpp::NumericVector sigma1,
//              arma::mat X1, arma::mat X2,
//              int copula){
//   int M = skelpts.size();
//   int N = sigma1.size()/M;
//   double s = 0.0;
//   for(int j=0; j<M; ++j){
//     for(int i=0; i<N; ++i){
//       int ij = (N*j + i);
//       Rcpp::NumericVector z_ij = z(ij, Rcpp::_);
//       arma::mat b1_ij = b1.row(ij).t();
//       arma::mat b2_ij = b2.row(ij).t();
//       double sigma1_ij = sigma1[ij];
//       s += exp(p_phi_log(eta, skelpts, phi, 
//                    y1, z_ij,
//                    b1_ij, b2_ij, sigma1_ij,
//                    X1, X2, 
//                    copula));
//     }
//   }
//   return(s);
// }
// // Bayes factor estimator, return a vector of BFs for a vector of phi
// // [[Rcpp::export]]
// Rcpp::NumericVector Bfn_vec(Rcpp::NumericVector phi, Rcpp::NumericVector eta, Rcpp::NumericVector skelpts,
//            Rcpp::NumericVector y1, Rcpp::NumericMatrix z,
//            arma::mat b1, arma::mat b2,
//            Rcpp::NumericVector sigma1,
//            arma::mat X1, arma::mat X2,
//            int copula){
//   int n_phi = phi.size();
//   Rcpp::NumericVector bf(n_phi);
//   for(int i=0; i<n_phi; ++i){
//     bf[i] = Bfn(phi[i], eta, skelpts,
//                 y1, z,
//                 b1, b2,
//                 sigma1,
//                 X1, X2, 
//                 copula);
//     printf("Finished BF: %d\n", i);
//   }
//   return(bf);
// }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Direct (non latent) functions                                      //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// // full pdf, direct, no latents
// // [[Rcpp::export]]
// double copulapdf_log_obs(Rcpp::NumericVector y1, Rcpp::NumericVector y2, 
//                   Rcpp::NumericVector mu1, Rcpp::NumericVector mu2, 
//                   double sigma1, double phi, int copula){
//   double logl = 0.0;
//   double z;
//   double p_b_oneminus;
//   double t1;
//   double t2;
//   double t3;
//   int n = y1.size();
//   for (int i=0; i<n; ++i){
//     z = (y1[i] - mu1[i])/sigma1;
//     p_b_oneminus = 1 / (1 + exp(mu2[i]));
//     t1 = -log(sigma1);
//     t2 = Rcpp::stats::dnorm_0(z, 1);
//     if (y2[i] == 0) {
//       t3 = log(dCphi_du2(p_b_oneminus, Rcpp::stats::pnorm_0(z, 1, 0), phi, copula));
//     } else {
//       t3 = log(1-dCphi_du2(p_b_oneminus, Rcpp::stats::pnorm_0(z, 1, 0), phi, copula));
//     }
//     logl += t1 + t2 + t3;
//   }
//   return logl;
// }

// [[Rcpp::export]]
double copulapdf_log_obs(Rcpp::NumericVector y1, Rcpp::NumericVector y2,
                         arma::vec mu1, arma::vec mu2,
                         double sigma1, double phi, int copula){
  double logl = 0.0;
  double z;
  double p_b_oneminus;
  double t1;
  double t2;
  double t3;
  int n = y1.size();
  for (int i=0; i<n; ++i){
    z = (y1[i] - mu1(i))/sigma1;
    p_b_oneminus = 1 / (1 + exp(mu2(i)));
    t1 = -log(sigma1);
    t2 = Rcpp::stats::dnorm_0(z, 1);
    if (y2[i] == 0) {
      t3 = log(dCphi_du2(p_b_oneminus, Rcpp::stats::pnorm_0(z, 1, 0), phi, copula));
    } else {
      t3 = log(1-dCphi_du2(p_b_oneminus, Rcpp::stats::pnorm_0(z, 1, 0), phi, copula));
    }
    logl += t1 + t2 + t3;
  }
  return logl;
}

// [[Rcpp::export]]
double copulapdf_wrapcov_log_obs(Rcpp::NumericVector y1, Rcpp::NumericVector y2, 
                          arma::mat X1, arma::mat X2,
                          arma::mat b1, arma::mat b2, 
                          double sigma1, double phi, int copula){
  // IMP: arma::vec saves much more memory than Rcpp::NumericVector
  arma::vec m1 = X1*b1;
  arma::vec m2 = X2*b2;
  return(copulapdf_log_obs(y1, y2, m1, m2, sigma1, phi, copula));
}

// log(p_j(x, eta)) from Geyer
// j will be on 0-indexing (add 1 to compare to paper formulas)
// [[Rcpp::export]]
double p_j_log_obs(Rcpp::NumericVector eta, Rcpp::NumericVector skelpts, int j, 
               Rcpp::NumericVector y1, Rcpp::NumericVector y2,
               arma::mat b1, arma::mat b2, double sigma1,
               arma::mat X1, arma::mat X2,
               int copula) {
  // create vector to feed into log_sum_exp;
  int M = skelpts.size();
  Rcpp::NumericVector lse(M);
  double log_h_j = copulapdf_wrapcov_log_obs(y1, y2, X1, X2, b1, b2, sigma1,
                                         skelpts[j], copula);
  for(int k=0; k<M; ++k){
    lse[k] = copulapdf_wrapcov_log_obs(y1, y2, X1, X2, b1, b2, sigma1,
                                   skelpts[k], copula) -
                                     log_h_j + eta[k] - eta[j];
  }
  return(-log_sum_exp(lse));
}

// function to maximize for estimating r
// Here, the parameter matrices and z have dimensions (NM x p), where p is p1 or p2 or n etc.
// They should have been supplied as row concatenated scalars/vectors from the MCMC samples
// The structure of concatenation is N, N, N,...N (M times).
// For maximization, provide 2nd element onward for eta. First is assumed zero wlog.
// [[Rcpp::export]]
double etafn_obs(Rcpp::NumericVector eta_reduced, Rcpp::NumericVector skelpts,
             Rcpp::NumericVector y1, Rcpp::NumericVector y2,
             arma::mat b1, arma::mat b2,
             Rcpp::NumericVector sigma1,
             arma::mat X1, arma::mat X2,
             int copula){
  // augment the eta. wlog first element of eta is zero.
  int n_eta = eta_reduced.size();
  Rcpp::NumericVector eta(n_eta+1); //defaults of zero
  for(int i=0; i<n_eta; ++i) eta[i+1] = eta_reduced[i];
  
  int M = skelpts.size();
  if(M != n_eta + 1) throw std::invalid_argument("In etafn(), length of eta_reduced should be number of skeleton points - 1");
  int N = sigma1.size()/M;
  double s = 0.0;
  for(int j=0; j<M; ++j){
    for(int i=0; i<N; ++i){
      int ij = (N*j + i);
      arma::mat b1_ij = b1.row(ij).t();
      arma::mat b2_ij = b2.row(ij).t();
      double sigma1_ij = sigma1[ij];
      s += p_j_log_obs(eta, skelpts, j, 
                   y1, y2,
                   b1_ij, b2_ij, sigma1_ij,
                   X1, X2, 
                   copula);
    }
  }
  return(s);
}

// gradient of etafn, return vector of length length(eta_reduced).
// note that we omit the first element of the gradient since we assume wlog that eta[0] = constant (0).
// [[Rcpp::export]]
Rcpp::NumericVector etafn_gr_obs(Rcpp::NumericVector eta_reduced, Rcpp::NumericVector skelpts,
                             Rcpp::NumericVector y1, Rcpp::NumericVector y2,
                             arma::mat b1, arma::mat b2,
                             Rcpp::NumericVector sigma1,
                             arma::mat X1, arma::mat X2,
                             int copula){
  // augment the eta. wlog first element of eta is zero.
  int n_eta = eta_reduced.size();
  Rcpp::NumericVector eta(n_eta+1); //defaults of zero
  for(int i=0; i<n_eta; ++i) eta[i+1] = eta_reduced[i];
  
  int M = skelpts.size();
  if(M != n_eta + 1) throw std::invalid_argument("In etafn(), length of eta_reduced should be number of skeleton points - 1");
  int N = sigma1.size()/M;
  Rcpp::NumericVector gr(n_eta);
  
  for(int r=1; r<n_eta+1; ++r){
    double s = 0.0;
    for(int j=0; j<M; ++j){
      for(int i=0; i<N; ++i){
        int ij = (N*j + i);
        arma::mat b1_ij = b1.row(ij).t();
        arma::mat b2_ij = b2.row(ij).t();
        double sigma1_ij = sigma1[ij];
        s += exp(p_j_log_obs(eta, skelpts, r,
                             y1, y2,
                             b1_ij, b2_ij, sigma1_ij,
                             X1, X2,
                             copula));
      }
    }
    gr[r-1] = N-s;
  }
  return(gr);
}

// Component in sum of Bayes factor estimator
// Does not have N_k factor in lse as in Roy, because we assume N_k is constant
// ie same number of iterations per skel pt, hence it can be taken outside and is
// irrelevant for maximization of Bfn
// [[Rcpp::export]]
double p_phi_log_obs(Rcpp::NumericVector eta, Rcpp::NumericVector skelpts, double phi, 
                 Rcpp::NumericVector y1, Rcpp::NumericVector y2,
                 arma::mat b1, arma::mat b2, double sigma1,
                 arma::mat X1, arma::mat X2,
                 int copula){
  // create vector to feed into log_sum_exp;
  int M = skelpts.size();
  Rcpp::NumericVector lse(M);
  double log_h_phi = copulapdf_wrapcov_log_obs(y1, y2, X1, X2, b1, b2, sigma1, phi,
                                           copula);
  for(int k=0; k<M; ++k){
    lse[k] = copulapdf_wrapcov_log_obs(y1, y2, X1, X2, b1, b2, sigma1, skelpts[k],
                                   copula) - 
                                     log_h_phi + eta[k] - eta[0];
  }
  return(-log_sum_exp(lse));
}


// Bayes factor estimator
// [[Rcpp::export]]
double Bfn_obs(double phi, Rcpp::NumericVector eta, Rcpp::NumericVector skelpts,
           Rcpp::NumericVector y1, Rcpp::NumericVector y2,
           arma::mat b1, arma::mat b2,
           Rcpp::NumericVector sigma1,
           arma::mat X1, arma::mat X2,
           int copula){
  int M = skelpts.size();
  int N = sigma1.size()/M;
  double s = 0.0;
  for(int j=0; j<M; ++j){
    for(int i=0; i<N; ++i){
      int ij = (N*j + i);
      arma::mat b1_ij = b1.row(ij).t();
      arma::mat b2_ij = b2.row(ij).t();
      double sigma1_ij = sigma1[ij];
      s += exp(p_phi_log_obs(eta, skelpts, phi, 
                         y1, y2,
                         b1_ij, b2_ij, sigma1_ij,
                         X1, X2, 
                         copula));
    }
  }
  return(s);
}
// Bayes factor estimator, return a vector of BFs for a vector of phi
// [[Rcpp::export]]
Rcpp::NumericVector Bfn_vec_obs(Rcpp::NumericVector phi, Rcpp::NumericVector eta, Rcpp::NumericVector skelpts,
                            Rcpp::NumericVector y1, Rcpp::NumericVector y2,
                            arma::mat b1, arma::mat b2,
                            Rcpp::NumericVector sigma1,
                            arma::mat X1, arma::mat X2,
                            int copula){
  int n_phi = phi.size();
  Rcpp::NumericVector bf(n_phi);
  for(int i=0; i<n_phi; ++i){
    bf[i] = Bfn_obs(phi[i], eta, skelpts,
                y1, y2,
                b1, b2,
                sigma1,
                X1, X2, 
                copula);
    printf("Finished BF: %d\n", i);
  }
  return(bf);
}

////////////////////////////////////////////////////////////////////////////////
//
//                   functions for standard error calculations
//
////////////////////////////////////////////////////////////////////////////////
// Following the notation from Roy's standard error paper
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
//                     Comment/uncomment from below to end.
//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// // Matrix of r's (ratios of normalizing contants)
// // [[Rcpp::export]]
// arma::mat D(arma::vec eta){
//   int M = eta.n_elem;
//   arma::mat D(M,M-1);
//   arma::vec d = arma::exp(-eta.tail(M-1)+eta(0));
//   D.head_rows(1) = d.t();
//   D.tail_rows(M-1) = -arma::diagmat(d);
//   return(D);
// }
// 
// // This function get's all the Y's and stores as a cube
// // Indices: 0 (i): MCMC iterations, 1 (r): skelpt in likelihood, 2 (l): skelpt in MCMC
// // Inputs are similar to those for etafn
// // Unlike etafn, we feed in the full eta, not reduced (first element is probably zero)
// // [[Rcpp::export]]
// arma::cube Y(Rcpp::NumericVector eta, Rcpp::NumericVector skelpts,
//                  Rcpp::NumericVector y1, Rcpp::NumericMatrix z,
//                  arma::mat b1, arma::mat b2,
//                  Rcpp::NumericVector sigma1,
//                  arma::mat X1, arma::mat X2,
//                  int copula) {
//   int M = skelpts.size();
//   int N = sigma1.size()/M; // no of MCMC iterations per skeleton point
//   arma::cube Y(N,M,M);
//   for(int r=0; r<M; ++r){
//     for(int l=0; l<M; ++l){
//       for(int i=0; i<N; ++i){
//         int il = (N*l + i);
//         Rcpp::NumericVector z_il = z(il, Rcpp::_);
//         arma::mat b1_il = b1.row(il).t();
//         arma::mat b2_il = b2.row(il).t();
//         double sigma1_il = sigma1[il];
//         Y(i,r,l) = exp(p_j_log(eta, skelpts, r,
//                        y1, z_il,
//                        b1_il, b2_il, sigma1_il,
//                        X1, X2,
//                        copula));
//       }
//     }
//   }
//   return(Y);
// }
// 
// // Get W's and store in matrix
// // Indices: r, l correspond to Y(i,r,l)
// // [[Rcpp::export]]
// arma::mat W(arma::cube Y){
//   arma::mat W = mean(Y, 0);
//   return(W);
// }
// 
// // Get B's and store in matrix
// // Indices: B(r,r) and B(s,s). Both r and s refer to skelpt in likelihood
// // [[Rcpp::export]]
// arma::mat B(arma::cube Y){
//   int N = Y.n_rows;
//   int M = Y.n_cols;
//   int NM = N*M;
//   arma::mat B(M,M);
//   for(int r=0; r<M; ++r){
//     for(int s=0; s<M; ++s){
//       if(r==s){
//         arma::mat P_r = Y(arma::span(), arma::span(r), arma::span());
//         B(r,r) = arma::accu(P_r % (1-P_r))/NM;
//       }
//       else{
//         arma::mat P_r = Y(arma::span(), arma::span(r), arma::span());
//         arma::mat P_s = Y(arma::span(), arma::span(s), arma::span());
//         B(r,s) = -arma::accu(P_r % P_s)/NM;
//       }
//     }
//   }
//   return(B);
// }
// 
// // Get Y_bar, the cube of batch means Y_bar(m, r, l)
// // Indices 0 (m): batch mean number (1,...,e) 2 (r) 3(l) have usual meanings
// // Don't assume that N is a perfect square
// // If N is not a perfect square, the last batch of size < e is discarded
// // [[Rcpp::export]]
// arma::cube Y_bar(arma::cube Y){
//   int N = Y.n_rows;
//   int M = Y.n_cols;
//   int e = floor(sqrt(N));
//   arma::cube Y_bar(e, M, M);
//   for(int m=0; m<e; ++m){
//     Y_bar(arma::span(m,m), arma::span(), arma::span()) =
//       mean(Y(arma::span(m*e,((m+1)*e-1)), arma::span(), arma::span()), 0);
//   }
//   return(Y_bar);
// }
// 
// // Sigma matrix, for diagnostics
// // Put l from 0,...,M-1
// // [[Rcpp::export]]
// arma::mat Sigma_l(arma::cube Y_bar, arma::mat W, int l){
//   int e = Y_bar.n_rows;
//   int M = Y_bar.n_cols;
//   arma::mat Sigma_l(M,M,arma::fill::zeros);
//   //for(int l=0; l<M; ++l){
//     for(int m=0; m<e; ++m){
//       arma::vec d = arma::vectorise(Y_bar(arma::span(m), arma::span(), arma::span(l))) -
//         arma::vectorise(W(arma::span(), arma::span(l)));
//       Sigma_l += d * d.t();
//     }
//   //}
//   Sigma_l *= e/((e-1));
//   return(Sigma_l);
// }
// 
// // Omega matrix. 
// // Omega = \sum{Sigma_1,...,Sigma_M} / M
// // [[Rcpp::export]]
// arma::mat Omega(arma::cube Y_bar, arma::mat W){
//   int M = Y_bar.n_cols;
//   arma::mat Omega(M,M,arma::fill::zeros);
//   for(int l=0; l<M; ++l){
//     arma::mat Sigma_l_ = Sigma_l(Y_bar, W, l);
//     Omega += Sigma_l_;
//   }
//   Omega /= M;
//   return(Omega);
// }
// 
// // Lambda matrix
// // [[Rcpp::export]]
// arma::mat Lambda(arma::mat D, arma::mat B, arma::mat Omega){
//   if(arma::approx_equal(Omega, arma::zeros(arma::size(Omega)), "absdiff", 0.0001)){
//     return(arma::zeros<arma::mat>(D.n_cols, D.n_cols));
//   }
//   else{
//     return(D.t()*solve(B, Omega)*solve(B, D));
//   }
// }
// 
// // Store V(i,l) in matrix V
// // function of unknown phi
// // z, b1, b2, sigma1 below must be from the NEW chain, ie not the chain
// // used to estimate eta
// // [[Rcpp::export]]
// arma::mat V(Rcpp::NumericVector eta, Rcpp::NumericVector skelpts,
//             Rcpp::NumericVector y1, Rcpp::NumericMatrix z,
//             arma::mat b1, arma::mat b2,
//             Rcpp::NumericVector sigma1,
//             arma::mat X1, arma::mat X2,
//             double phi, int copula){
//   int M = skelpts.size();
//   int N = sigma1.size()/M; // no of MCMC iterations per skeleton point
//   arma::mat V(N,M);
//   for(int l=0; l<M; ++l){
//     for(int i=0; i<N; ++i){
//       int il = (N*l + i);
//       Rcpp::NumericVector z_il = z(il, Rcpp::_);
//       arma::mat b1_il = b1.row(il).t();
//       arma::mat b2_il = b2.row(il).t();
//       double sigma1_il = sigma1[il];
//       V(i,l) = M*exp(p_phi_log(eta, skelpts, phi,
//                      y1, z_il,
//                      b1_il, b2_il, sigma1_il,
//                      X1, X2,
//                      copula));
//     }
//   }
//   return(V);
// }
// 
// // V_bar batch means of V
// // Matrix of dim e x M
// // Indices V_bar(m, l)
// // [[Rcpp::export]]
// arma::mat V_bar(arma::mat V){
//   int N = V.n_rows;
//   int M = V.n_cols;
//   int e = floor(sqrt(N));
//   arma::mat V_bar(e, M);
//   for(int m=0; m<e; ++m){
//     V_bar(arma::span(m,m), arma::span()) =
//       mean(V(arma::span(m*e,((m+1)*e-1)), arma::span()), 0);
//   }
//   return(V_bar);
// }
// 
// // tau (actually tau^2)
// // [[Rcpp::export]]
// double tausq(arma::mat V_bar){
//   int e = V_bar.n_rows;
//   int M = V_bar.n_cols;
//   arma::vec tau_vec(M);
//   arma::vec V_bar_bar = mean(V_bar, 0).t(); // grand mean
//   // Compute tau[l] for each l
//   for(int l=0; l<M; ++l){
//     tau_vec(l) = (e/(e-1))*arma::accu( arma::square(V_bar.col(l) - V_bar_bar(l)) );
//   }
//   return(mean(tau_vec));
// }
// 
// // c vector
// // [[Rcpp::export]]
// arma::vec c_vec(Rcpp::NumericVector eta, Rcpp::NumericVector skelpts,
//                 Rcpp::NumericVector y1, Rcpp::NumericMatrix z,
//                 arma::mat b1, arma::mat b2,
//                 Rcpp::NumericVector sigma1,
//                 arma::mat X1, arma::mat X2,
//                 double phi, int copula){
//   int M = skelpts.size();
//   if(M != eta.size()) throw std::invalid_argument("In etafn(), eta.size() should be equal to number of skeleton points");
//   int N = sigma1.size()/M; // no of MCMC iterations per skeleton point
//   arma::vec c_vec(M-1, arma::fill::zeros);
//   for(int r=1; r<M; ++r){ // r here is j in Roy's formula. Also start from 1, not 0
//     double d_r = exp(eta[0] - eta[r]);
//     for(int l=0; l<M; ++l){
//       for(int i=0; i<N; ++i){
//         int il = (N*l + i);
//         Rcpp::NumericVector z_il = z(il, Rcpp::_);
//         arma::mat b1_il = b1.row(il).t();
//         arma::mat b2_il = b2.row(il).t();
//         double sigma1_il = sigma1[il];
//         c_vec(r-1) += exp( p_j_log(eta, skelpts, r,
//                                    y1, z_il,
//                                    b1_il, b2_il, sigma1_il,
//                                    X1, X2, copula) +
//                            p_phi_log(eta, skelpts, phi,
//                                      y1, z_il,
//                                      b1_il, b2_il, sigma1_il,
//                                      X1, X2, copula) );
//       }
//     }
//     c_vec(r-1) /= (N*d_r*d_r);
//   }
//   return(c_vec);
// }
// 
// // Final standard error of Bayes Factor
// // [[Rcpp::export]]
// double stderr_BF(int N_new_total, double q, arma::vec c_vec, arma::mat Lambda,
//                  double tausq){
//   return(sqrt((q*(arma::as_scalar(c_vec.t()*Lambda*c_vec)) + tausq) / N_new_total));
// }
// 
// // main function to compute standard error at each of a vector of phis
// // [[Rcpp::export]]
// Rcpp::NumericVector se_BF_wrap(Rcpp::NumericVector phi,
//                                Rcpp::NumericVector eta, Rcpp::NumericVector skelpts,
//                                Rcpp::NumericVector y1,
//                                Rcpp::NumericMatrix z_old, Rcpp::NumericMatrix z_new,
//                                arma::mat b1_old, arma::mat b1_new,
//                                arma::mat b2_old, arma::mat b2_new,
//                                Rcpp::NumericVector sigma1_old, Rcpp::NumericVector sigma1_new,
//                                arma::mat X1, arma::mat X2,
//                                int copula){
//   // Compute old and new N's, the number of iterations for each skep pts MCMC
//   int M = skelpts.size();
//   if(M != eta.size()) throw std::invalid_argument("In etafn(), eta.size() should be equal to number of skeleton points");
//   int N_old = sigma1_old.size()/M;
//   int N_new = sigma1_new.size()/M;
//   double q = (float) N_new / (float) N_old;
//   Rcpp::NumericVector se_BF(phi.size());
//   for(int g=0; g<phi.size(); ++g){
//     //--Lambda
//     //----D
//     //----B
//     //------Y
//     //----Omega
//     //------Y_bar
//     //--------Y
//     //------W
//     //--------Y
//     arma::mat DD = D(eta);
//     // DD.print();
//     // printf("\nFinished D\n");
//     arma::cube YY = Y(eta, skelpts,
//                       y1, z_old,
//                       b1_old, b2_old,
//                       sigma1_old,
//                       X1, X2,
//                       copula);
//     // printf("\nFinished Y\n");
//     arma::mat BB = B(YY);
//     // BB.print();
//     // printf("\nFinished B\n");
//     arma::cube Ybar = Y_bar(YY);
//     // printf("\nFinished Y_bar\n");
//     arma::mat WW = W(YY);
//     // WW.print();
//     // printf("\nFinished W\n");
//     arma::mat OO = Omega(Ybar, WW);
//     // OO.print();
//     // printf("\nFinished Omega\n");
//     arma::mat L = Lambda(DD, BB, OO);
//     // L.print();
//     // printf("Finished Lambda\n");
//     //--tausq
//     //----V_bar
//     //----V
//     arma::mat Vmat = V(eta, skelpts,
//                     y1, z_new,
//                     b1_new, b2_new,
//                     sigma1_new,
//                     X1, X2,
//                     phi[g], copula);
//     // printf("\nFinished V\n");
//     arma::mat Vbar = V_bar(Vmat);
//     //printf("\nFinished V_bar\n");
//     double tsq = tausq(Vbar);
//     //printf("\ntausq: %f\n", tsq);
//     //printf("\nFinished tau\n");
//     //--c_vec
//     arma::vec cvec = c_vec(eta, skelpts,
//                            y1, z_new,
//                            b1_new, b2_new,
//                            sigma1_new,
//                            X1, X2,
//                            phi[g], copula);
//     //cvec.print();
//     //printf("\nFinished c_vec\n");
//     //se
//     se_BF[g] = stderr_BF(N_new*M, q, cvec, L, tsq);
//     printf("Finished se(BF): %d\n", g);
//   }
//   return(se_BF);
// }
// 
// // First term of varinace of Bayes Factor (involving q) for diagnostics
// // [[Rcpp::export]]
// double stderr_BF_term1(int N_new_total, double q, arma::vec c_vec, arma::mat Lambda){
//   return((q*(arma::as_scalar(c_vec.t()*Lambda*c_vec))) / N_new_total);
// }
// 
// // Second term of variance of Bayes Factor (not involving q) for diagnostics
// // [[Rcpp::export]]
// double stderr_BF_term2(int N_new_total, double tausq){
//   return( tausq / N_new_total );
// }
// 
// // diagnostic function to return separate parts of the variance of BF
// // [[Rcpp::export]]
// arma::mat se_BF_twoterms(Rcpp::NumericVector phi,
//                                Rcpp::NumericVector eta, Rcpp::NumericVector skelpts,
//                                Rcpp::NumericVector y1,
//                                Rcpp::NumericMatrix z_old, Rcpp::NumericMatrix z_new,
//                                arma::mat b1_old, arma::mat b1_new,
//                                arma::mat b2_old, arma::mat b2_new,
//                                Rcpp::NumericVector sigma1_old, Rcpp::NumericVector sigma1_new,
//                                arma::mat X1, arma::mat X2,
//                                int copula){
//   // Compute old and new N's, the number of iterations for each skep pts MCMC
//   int M = skelpts.size();
//   if(M != eta.size()) throw std::invalid_argument("In etafn(), eta.size() should be equal to number of skeleton points");
//   int N_old = sigma1_old.size()/M;
//   int N_new = sigma1_new.size()/M;
//   double q = (float) N_new / (float) N_old;
//   arma::mat se_BF(phi.size(), 2);
//   for(int g=0; g<phi.size(); ++g){
//     //--Lambda
//     //----D
//     //----B
//     //------Y
//     //----Omega
//     //------Y_bar
//     //--------Y
//     //------W
//     //--------Y
//     arma::mat DD = D(eta);
//     // DD.print();
//     // printf("\nFinished D\n");
//     arma::cube YY = Y(eta, skelpts,
//                       y1, z_old,
//                       b1_old, b2_old,
//                       sigma1_old,
//                       X1, X2,
//                       copula);
//     // printf("\nFinished Y\n");
//     arma::mat BB = B(YY);
//     // BB.print();
//     // printf("\nFinished B\n");
//     arma::cube Ybar = Y_bar(YY);
//     // printf("\nFinished Y_bar\n");
//     arma::mat WW = W(YY);
//     // WW.print();
//     // printf("\nFinished W\n");
//     arma::mat OO = Omega(Ybar, WW);
//     // OO.print();
//     // printf("\nFinished Omega\n");
//     arma::mat L = Lambda(DD, BB, OO);
//     // L.print();
//     // printf("Finished Lambda\n");
//     //--tausq
//     //----V_bar
//     //----V
//     arma::mat Vmat = V(eta, skelpts,
//                        y1, z_new,
//                        b1_new, b2_new,
//                        sigma1_new,
//                        X1, X2,
//                        phi[g], copula);
//     // printf("\nFinished V\n");
//     arma::mat Vbar = V_bar(Vmat);
//     //printf("\nFinished V_bar\n");
//     double tsq = tausq(Vbar);
//     //printf("\ntausq: %f\n", tsq);
//     //printf("\nFinished tau\n");
//     //--c_vec
//     arma::vec cvec = c_vec(eta, skelpts,
//                            y1, z_new,
//                            b1_new, b2_new,
//                            sigma1_new,
//                            X1, X2,
//                            phi[g], copula);
//     //cvec.print();
//     //printf("\nFinished c_vec\n");
//     //se
//     se_BF(g, 0) = stderr_BF_term1(N_new*M, q, cvec, L);
//     se_BF(g, 1) = stderr_BF_term2(N_new*M, tsq);
//     printf("Finished se(BF) diagnostics: %d\n", g);
//   }
//   return(se_BF);
// }

////////////////////////////////////////////////////////////////////////////////
//
//                   functions for model selection
//
////////////////////////////////////////////////////////////////////////////////

// full pdf, direct, no latents, for one observation
// [[Rcpp::export]]
double copulapdf_log_obs_single(double y1, double y2, 
                         double mu1, double mu2, 
                         double sigma1, double phi, int copula){
  double z = (y1 - mu1)/sigma1;
  double p_b_oneminus = 1 / (1 + exp(mu2));
  if (y2 == 0) {
    return (-log(sigma1) + Rcpp::stats::dnorm_0(z, 1) + 
            log(dCphi_du2(p_b_oneminus, Rcpp::stats::pnorm_0(z, 1, 0), phi, copula)));
  } else {
    return (-log(sigma1) + Rcpp::stats::dnorm_0(z, 1) + 
            log(1-dCphi_du2(p_b_oneminus, Rcpp::stats::pnorm_0(z, 1, 0), phi, copula)));
  }
}

// lpml
// [[Rcpp::export]]
double lpml_Cpp(Rcpp::NumericVector y1, Rcpp::NumericVector y2,
                arma::mat mu1, arma::mat mu2,
                Rcpp::NumericVector sigma1, Rcpp::NumericVector phi,
                int copula){
  int n = y1.size();
  int B = phi.size();
  double LPML = 0.0;
  for (int i=0; i<n; i++) {
    double s = 0.0;
    for (int b=0; b<B; b++) {
      s += exp(-copulapdf_log_obs_single(y1[i], y2[i], mu1(b,i), mu2(b,i), 
                                  sigma1[b], phi[b], copula));
    }
    LPML += log(s);
  }
  return(n*log(B) - LPML);
}