#!/home/statsadmin/R/bin/Rscript

rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)

source('simulate_data.R') # source code for simulation of the covariates and bivariate response
source('corrmat.R') # incorporating possible correlation in the covariates
source('decodeCombo.R') # Decode combination code to get simulation settings
source('copulalinks.R') # Different transform of copula parameter for each copula family
library(copula)
library(rstan)
library(optimx)
library(Rcpp)
sourceCpp('Cppfns.cpp') # copula functions written in C++ for effective optimization

# Get simulation settings
SEED <- 123
extCombo <- 111 # n=500, tau=0.1, Gaussian copula
# (N)(tau)(cop.int)
# N: 1=500, 2=1000
# tau: tau/10 .1, .3, .6
# cop.int: 1=normal, 2=clayton, 3=gumbel, 4=frank


settings <- decodeCombo(extCombo)
N <- settings$N
tau <- settings$tau
cop.int <- settings$cop.int
cop.str <- switch(cop.int, "normal", "clayton", "gumbel", "frank")
#### DATA GENERATION #######################################################
sd <- SEED
print(paste0("Seed is ", sd))
set.seed(sd)
p <- 2 # No of covariates
m <- 2 # Cluster size (equal across clusters)
Btn <- c(1.2, 2.0)
Btb <- c(1.2, .8)
Bt <- c(Btn, Btb)

# set phi according to kendall's tau
tmp <- switch(cop.int, normalCopula(), claytonCopula(), gumbelCopula(), 
              frankCopula())
phi.true <- iTau(tmp, tau)

datlist <- gen_mixed_data(Bt, N, rho=phi.true, intercept=TRUE, cov="same", 
                          coptype=cop.str)
X <- datlist$X
Z <- datlist$Z
y <- datlist$y
y_c <- y[seq(1,m*N,2)]; y_b <- y[seq(2,m*N,2)]
rm(y);rm(datlist);gc()

#### Full Bayes ###########################
print(paste0("Starting Full Bayes"))
fb_start <- Sys.time()  # record starting time
dat <- list(y1=y_c, y2=y_b, x1=X, x2=Z, k1=p, k2=p, n=N, copula=cop.int)
fit <- stan(file = 'FB.stan', data=dat, iter=2000, chains=2, warmup=1000)
rm(list = c("dat"));gc()
fb_end <- Sys.time()  # record ending time
fb_time <- fb_end - fb_start
# extract samples for b1, b2, z, sigma
pars_ <- c('b1','b2','sigma1','phi')
samps <- extract(fit, pars_)
rm(list = c("fit"));gc()
results.fb <- list( samples=samps,
                    fb_time=fb_time,
                    n=N,
                    p=p,
                    y_c=y_c,
                    y_b=y_b,
                    X_c=X,
                    X_b=Z,
                    phi_true=phi.true,
                    beta_c=Btn,
                    beta_b=Btb,
                    seed=sd)
save(results.fb, file = paste0("FILE-FB-", extCombo, ".RData"))

# free up some memory
rm(list = c("samps", "results.fb")); gc()


#### Empirical Bayes ###############
print(paste0("Starting Empirical Bayes, direct"))

eb_start <- Sys.time()

# Assume skeleton points
# transform to real scale, then take skeleton, then untransform
phi.true.trans <- copula_link(phi.true, cop.int)
skel.trans <- phi.true.trans*c(1.0, 0.5, 0.8, 1.2, 1.5)
skel <- copula_invlink(skel.trans, cop.int)
nskel <- length(skel)

## Estimate r's ##
# Generate MCMC chains for each skeleton point
print(paste0("Starting chains for eta estimation"))
# r_MCMC <- list()
r_samps <- list()
iter_ <- 1500
thin_ <- 1
chains_ <- 2
warmup_ <- 500
dat <- list(y1=y_c, y2=y_b, x1=X, x2=Z, k1=p, k2=p, n=N, phi=NA,
            copula=cop.int)
# extract and stack coefficients.
# Parameters should be (iter_*nskel) x p1 (or p2 or N)
# pars_ <- c('b1','b2','sigma1')
pars.nrow <- (iter_ - warmup_)*chains_
r_samps[["b1"]] <- matrix(, nrow = pars.nrow*nskel, ncol = p)
r_samps[["b2"]] <- matrix(, nrow = pars.nrow*nskel, ncol = p)
r_samps[["sigma1"]] <- numeric(pars.nrow*nskel)
for(i in seq_along(skel)){
  dat$phi <- skel[i]
  r_MCMC <- stan(file = 'EB.stan', data=dat, iter=iter_, chains=chains_,
                 thin=thin_, warmup=warmup_)
  rowids <- ((i-1)*pars.nrow+1):(i*pars.nrow)
  r_samps[["b1"]][rowids,] <- extract(r_MCMC, "b1")[[1]]
  r_samps[["b2"]][rowids,] <- extract(r_MCMC, "b2")[[1]]
  r_samps[["sigma1"]][rowids] <- extract(r_MCMC, "sigma1")[[1]]
}

# free up some memory
rm(list = c("r_MCMC", "dat")); gc()

print(paste0("Starting optimization for eta estimation"))

eb_s <- Sys.time()

# feed into etafn and optimize
eta_red_init <- rep(0, nskel-1)
eta_opt <- optimx(par=eta_red_init, fn=etafn_obs, gr=etafn_gr_obs, method='BFGS',
                  control=list(maximize=TRUE),
                  skelpts=skel, y1=y_c, y2=y_b,
                  b1=r_samps['b1'][[1]], b2=r_samps['b2'][[1]],
                  sigma1=r_samps['sigma1'][[1]],
                  X1=X, X2=Z, copula=cop.int)

eb_e <- Sys.time()

eb_e - eb_s

# extract eta
eta <- unlist(eta_opt[1:(nskel-1)])
eta <- c(0, eta)
Sys.time()
# free up some memory
rm(list = c("r_samps")); gc()
## B function, estimate phi 
print(paste0("Starting chains for phi estimation"))
# get fresh MCMC chains for B function
# r_MCMC <- list()
r_samps <- list()
iter_ <- 1500
chains_ <- 2
warmup_ <- 1000
dat <- list(y1=y_c, y2=y_b, x1=X, x2=Z, k1=p, k2=p, n=N, phi=NA, copula=cop.int)
# extract and stack coefficients.
# Parameters should be (iter_*nskel) x p1 (or p2 or N)
# pars_ <- c('b1','b2','sigma1')
pars.nrow <- (iter_ - warmup_)*chains_
r_samps[["b1"]] <- matrix(, nrow = pars.nrow*nskel, ncol = p)
r_samps[["b2"]] <- matrix(, nrow = pars.nrow*nskel, ncol = p)
r_samps[["sigma1"]] <- numeric(pars.nrow*nskel)
for(i in seq_along(skel)){
  dat$phi <- skel[i]
  r_MCMC <- stan(file = 'EB.stan', data=dat, iter=iter_, chains=chains_,
                 thin=thin_, warmup=warmup_)
  rowids <- ((i-1)*pars.nrow+1):(i*pars.nrow)
  r_samps[["b1"]][rowids,] <- extract(r_MCMC, "b1")[[1]]
  r_samps[["b2"]][rowids,] <- extract(r_MCMC, "b2")[[1]]
  r_samps[["sigma1"]][rowids] <- extract(r_MCMC, "sigma1")[[1]]
}

# free up some memory
rm(list = c("r_MCMC", "dat")); gc()

# now maximize Bfn to obtain EBayes estimate of phi
# search interval for phi depends on range of phi
# search on range of Tau (0, 0.8)
# 1=normal, 2=clayton, 3=gumbel, 4=frank
print(paste0("Starting optimization for phi estimation"))
search.interval <- switch(cop.int, c(-.99, .99),
                          c(0.001, 8),
                          c(1.001, 5),
                          c(0.001, 20))
phi_opt <- optimize(Bfn_obs, search.interval, maximum=TRUE,
                    eta=eta, skelpts=skel, y1=y_c, y2=y_b,
                    b1=r_samps['b1'][[1]], b2=r_samps['b2'][[1]],
                    sigma1=r_samps['sigma1'][[1]],
                    X1=X, X2=Z, copula=cop.int)
phi <- phi_opt$maximum

# free up some memory
rm(list = c("r_samps")); gc()
## Estimate remaining parameters ##
print(paste0("Starting final MCMC for empirical Bayes estimation"))
# Fit with Stan
dat <- list(y1=y_c, y2=y_b, x1=X, x2=Z, k1=p, k2=p, n=N, phi=phi, copula=cop.int)
fit <- stan(file = 'EB.stan', data=dat, iter=2000, chains=2, thin=1, warmup=1000,
            control = list(adapt_delta = 0.99))
eb_end <- Sys.time()
eb_time <- eb_end - eb_start

# extract samples for b1, b2, sigma1
samps <- extract(fit, c("b1","b2","sigma1"))
# free up some memory
rm(list = c("fit")); gc()
results.eb <- list(samples=samps, ## final MCMC samples for marginal parameters
                   eb_time=eb_time,
                   phi_est=phi,
                   eta=eta,
                   skelpts=skel,
                   eta_opt=eta_opt)   ## for convergence code
save(results.eb, file = paste0("FILE-EB-", extCombo, ".RData"))



