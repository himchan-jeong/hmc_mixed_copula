functions {
// copula
// 1: Gaussian
// 2: Clayton
// 3: Gumbel
// 4: Frank
real dCphi_du1du2_s(real u1, real u2, real phi, int copula) {
  real c;
  if (copula==4) { // Frank
    c = phi*(1-exp(-phi))*exp(-phi*(u1+u2))*((1-exp(-phi)) - (1-exp(-phi*u1))*(1-exp(-phi*u2)))^(-2);
  } else if (copula==3) { //Gumbel
    if (phi < 1) {
      reject("phi must be >= 1 for Gumbel copula");
    }
    c = exp(-(((-log(u1))^phi)+((-log(u2))^phi))^(1/phi)) * (1/(u1*u2)) * ((((-log(u1))^phi)+((-log(u2))^phi))^(-2+2/phi)) * ((log(u1)*log(u2))^(phi-1)) * (1 + (phi-1)*((((-log(u1))^phi)+((-log(u2))^phi))^(-1/phi)));
  } else if (copula==2) { //Clayton
    if (phi <= 0) {
      reject("phi must be > 0 for Clayton copula");
    }
    c = ((1+phi)*((u1*u2)^(-1-phi))) / ((u1^(-phi) + u2^(-phi) - 1)^(2+1/phi));
  } else if (copula==1) { //Gaussian
    // phi = rho, the correlation
    if (phi <= -1 || phi >= 1) {
      reject("phi must be between -1 and 1 for Gaussian copula");
    }
    c = exp(-(phi^2*(inv_Phi(u1)^2 + inv_Phi(u2)^2) - 2*phi*inv_Phi(u1)*inv_Phi(u2)) / (2*(1-(phi^2)))) / sqrt(1-phi^2);
  } else {
    reject("copula parameter must be an integer between 1 and 4.");
  }
  return c;
}

real copulapdf_log(real[] y1, real[] z, vector mu1, vector mu2, real sigma1, real phi, int n, int copula){
  real logl;
  real s;
  logl = 0.0;
  for (i in 1:n){
    s = log(dCphi_du1du2_s(normal_cdf(y1[i] , mu1[i], sigma1), 
                           logistic_cdf(z[i] , mu2[i], 1), phi,
                           copula)) + 
        normal_log(y1[i], mu1[i], sigma1) + 
        logistic_log(z[i], mu2[i], 1);
    logl = logl + s;
  }
  return logl;
}
real dCphi_du2(real u1, real u2, real phi, int copula) {
  real c;
  if (copula == 1) {  # Gaussian
    c = normal_cdf((inv_Phi(u1) - phi*inv_Phi(u2))/sqrt(1-pow(phi, 2)), 0, 1);
  } else if (copula == 2) {  # Clayton
    c = pow(u2, -phi-1) * pow(pow(u1, -phi)+pow(u2, -phi)-1, -1-1/phi);
  } else if (copula == 3) {
    c = -exp(-pow(pow(-log(u1), phi)+pow(-log(u2), phi), 1/phi)) * pow(pow(-log(u1), phi)+pow(-log(u2), phi), 1/phi-1) * pow(-log(u2), phi) / (u2 * log(u2));
  } else if (copula == 4) {
    c = -exp(phi) * (exp(phi*u1)-1) / (exp(phi*u2 + phi*u1) - exp(phi*u2 + phi) - exp(phi*u1 + phi) + exp(phi));
  } else {
    reject("copula parameter must be an integer between 1 and 4.");
  }
  return c;
}

real obslik_log(real[] y1, real[] y2, vector mu1, vector mu2, real sigma1, real phi, int n, int copula){
  real logl;
  real s;
  real z;
  real p_b_oneminus;
  real t1;
  real t2;
  real t3;
  logl = 0.0;
  for (i in 1:n){
    z = (y1[i] - mu1[i])/sigma1;
    p_b_oneminus = 1 / (1 + exp(mu2[i]));
    t1 = -log(sigma1);
    t2 = normal_lpdf(z | 0, 1);
    if (y2[i] == 0) {
      t3 = log(dCphi_du2(p_b_oneminus, normal_cdf(z, 0, 1), phi, copula));
    } else {
      t3 = log(1-dCphi_du2(p_b_oneminus, normal_cdf(z, 0, 1), phi, copula));
    }
    s = t1 + t2 + t3;
    logl = logl + s;
  }
  return logl;
}
}

// data {
//   int<lower=0> n; // number of subjects
//   int<lower=0> k1; // number of predictors for y1
//   int<lower=0> k2; // number of predictors for y2
//   real y1[n]; // continuous data
//   real y2[n]; // 0/1 binary data
//   matrix[n, k1] x1; // predictor variables for y1
//   matrix[n, k2] x2; // predictor variables for y2
//   int<lower=1, upper=4> copula; // copula type, integer from 1-4.
//   real phi; //copula parameter, estimated from EB method
// }
// 
// transformed data{
//   int<lower=-1, upper=1> sign[n];
//   for (i in 1:n) {
//     if (y2[i]==1) 
//       sign[i] = 1;
//     else
//       sign[i] = -1;
//   }
// }
// 
// parameters {
//   vector[k1] b1; // beta coefficients for y1
//   vector[k2] b2; // beta coefficients for y2
//   real<lower=0> abs_z[n]; // abs value of latent variable
//   real<lower=0> sigma1; // sd for y1's normal distribution
// }
// 
// transformed parameters {
//   real z[n];
//   vector[n] mu1; // location for y1
//   vector[n] mu2; // location for z
//   for (i in 1:n) {
//     z[i] = sign[i] * abs_z[i];
//   }
//   mu1 = x1 * b1;
//   mu2 = x2 * b2;
// }
// 
// model {
//   b1 ~ normal(0, 100);
//   b2 ~ normal(0, 100);
//   sigma1 ~ cauchy(0, 5); //half-Cauchy implicit
//   increment_log_prob(copulapdf_log(y1, z, mu1, mu2, sigma1, phi, n, copula));
// }

// ************************************ Direct model ******************************************
data {
  int<lower=0> n; // number of subjects
  int<lower=0> k1; // number of predictors for y1
  int<lower=0> k2; // number of predictors for y2
  real y1[n]; // continuous data
  real y2[n]; // 0/1 binary data
  matrix[n, k1] x1; // predictor variables for y1
  matrix[n, k2] x2; // predictor variables for y2
  int<lower=1, upper=4> copula; // copula type, integer from 1-4.
  real phi; //copula parameter, estimated from EB method
}

parameters {
  vector[k1] b1; // beta coefficients for y1
  vector[k2] b2; // beta coefficients for y2
  real<lower=0> sigma1; // sd for y1's normal distribution
}

transformed parameters {
  vector[n] mu1; // location for y1
  vector[n] mu2; // location for z
  mu1 = x1 * b1;
  mu2 = x2 * b2;
}

model {
  b1 ~ normal(0, 100);
  b2 ~ normal(0, 100);
  sigma1 ~ cauchy(0, 5); //half-Cauchy implicit
  increment_log_prob(obslik_log(y1, y2, mu1, mu2, sigma1, phi, n, copula));
}
// ********************************************************************************************
