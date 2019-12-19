# library(rstan)
# library(copula)
# library(VineCopula)
# sourceCpp("Cppfns.cpp")

# observed data log likelihood
llik_obs <- function(y_c, y_b, eta_c, eta_b, sigma_c, phi, cop.int) {
  n <- length(y_c)
  z_c <- (y_c - eta_c) / sigma_c
  p_b_oneminus <- 1 / (1 + exp(eta_b))
  if (cop.int == 1) {  # gaussian
    fam <- 1
  } else if (cop.int == 2) {  # clayton
    fam <- 3
  } else if (cop.int == 3) {  # gumbel
    fam <- 4
  } else if (cop.int == 4) {  # frank
    fam <- 5
  } else stop("cop.int must be either 1, 2, 3 or 4")
  t1 <- -log(sigma_c)  # times n factor will happen with sum()
  t2 <- dnorm(z_c, log = TRUE)
  t3 <- ifelse(y_b == 0,
               log(BiCopHfunc2(p_b_oneminus, pnorm(z_c),
                               fam, phi)),
               log(1 - BiCopHfunc2(p_b_oneminus, pnorm(z_c),
                                   fam, phi))
               )
  sum(t1 + t2 + t3)
}

# compute DIC1, based on observed likelihood
dic1 <- function(samps, dat, cop.int) {
  # samps: extract(fit)
  # dat: data list
  # cop.int: 1=Gaussian, 2=Clayton, 3=Gumbel, 4=Frank
  # compute average conditional log likelihood across draws
  mean_llik <- numeric(length(samps$phi))
  for (i in seq_along(mean_llik)) {
    mean_llik[i] <- llik_obs(dat$y1, dat$y2, samps$mu1[i,], samps$mu2[i,], 
                             samps$sigma1[i], samps$phi[i], cop.int)
  }
  mean_llik <- mean(mean_llik)
  
  # compute conditional log likelihood at posterior means
  llik_at_means <- llik_obs(dat$y1, 
                            dat$y2, 
                            colMeans(samps$mu1), 
                            colMeans(samps$mu2), 
                            mean(samps$sigma1),
                            mean(samps$phi), 
                            cop.int)
  
  # DIC1
  dic1 <- -4*mean_llik + 2*llik_at_means
  pd <- -2*mean_llik + 2*llik_at_means  # effective number of parameters
  list(dic1=dic1, pd=pd)
}

# compute DIC7, based on conditional likelihood (w/latent)
# -4*average_likelihood + 2*likelihood_at_average_parameters
# means refer to etas, not actual means
dic7 <- function(samps, dat, cop.int) {
  # samps: extract(fit)
  # dat: data list
  # cop.int: 1=Gaussian, 2=Clayton, 3=Gumbel, 4=Frank
  
  # compute average conditional log likelihood across draws
  mean_llik <- numeric(length(samps$phi))
  for (i in seq_along(mean_llik)) {
    mean_llik[i] <- copulapdf_log(dat$y1, samps$v[i,], samps$mu1[i,], samps$mu2[i,], 
                                  samps$sigma1[i], samps$phi[i], 
                                  dat$n, cop.int)
  }
  mean_llik <- mean(mean_llik)
  
  # compute conditional log likelihood at posterior means
  llik_at_means <- copulapdf_log(dat$y1, 
                                 colMeans(samps$v), 
                                 colMeans(samps$mu1), 
                                 colMeans(samps$mu2), 
                                 mean(samps$sigma1),
                                 mean(samps$phi), 
                                 dat$n, cop.int)
  
  # DIC7
  -4*mean_llik + 2*llik_at_means
}

# # LPML
# # LPML = n*log(B) - \sum\limits_{i=1}^n log(\sum\limits_{b=1}^B exp(-loglik(i,b)))
 lpml <- function(samps, dat, cop.int) {
   n <- length(dat$y1)  # number of observations
   B <- length(samps$phi)  # number of MCMC iterations
   LPML <- 0
   for (i in 1:n) {
   neg_llik_i <- numeric(B)  # vector of negative log likelihoods, the bth component is llik(i,b)
     for (b in 1:B) {
       neg_llik_i[b] <- -llik_obs(dat$y1[i], dat$y2[i], samps$mu1[b,i], samps$mu2[b,i], 
                     samps$sigma1[b], samps$phi[b], cop.int)
     }
   LPML <- LPML + log_sum_exp(neg_llik_i)
   }
   LPML <- n*log(B) - LPML
 }
 
# LPML = n*log(B) - \sum\limits_{i=1}^n log(\sum\limits_{b=1}^B exp(-loglik(i,b)))
# lpml2 <- function(samps, dat, cop.int) {
#   n <- length(dat$y1)  # number of observations
#   B <- length(samps$phi)  # number of MCMC iterations
#   LPML <- 0
#   for (i in 1:n) {
#     s <- 0  # vector of negative log likelihoods, the bth component is llik(i,b)
#     for (b in 1:B) {
#       s <- s + exp(-llik_obs(dat$y1[i], dat$y2[i], samps$mu1[b,i], samps$mu2[b,i], 
#                                  samps$sigma1[b], samps$phi[b], cop.int))
#     }
#     LPML <- LPML + log(s)
#   }
#   LPML <- n*log(B) - LPML
# }