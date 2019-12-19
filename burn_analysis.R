
library(data.table)
library(dplyr)
library(copula)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(optimx)
library(Rcpp)
library(bayesplot)
library(VineCopula)
library(gridExtra)

source('simulate_data.R') # source code for simulation of the covariates and bivariate response
source('corrmat.R') # incorporating possible correlation in the covariates
source('decodeCombo.R') # Decode combination code to get simulation settings
source('copulalinks.R') # Different transform of copula parameter for each copula family

sourceCpp('Cppfns.cpp') # copula functions written in C++ for effective optimization
source('model_selection.R') # Calculate DIC and LPML for model selection

# read in the data and generate the outcomes
burn <- read.table(url("http://www-personal.umich.edu/~pxsong/Burn.dat.txt"))[,c(2,10,30)]
names(burn) <- c("age", "total_burn_area", "survival_status")
burn <- burn %>% mutate(ln_burn_area = log(total_burn_area+1))
burn$survival_status[burn$survival_status == 2] <- 0
# create standardized version of age, better for MCMC
burn <- burn %>% mutate(age_std = scale(age))

# fit models with Gaussian, Frank, Gumbel and Clayton copulas
dat <- with(burn, 
            list(y1=ln_burn_area, 
                 y2=survival_status, 
                 x1=cbind(1, age_std), 
                 x2=cbind(1, age_std), 
                 k1=2, 
                 k2=2, 
                 n=length(age), 
                 copula=NULL))
# Gaussian
gauss.start <- Sys.time()
dat$copula <- 1
gaussian <- stan(file = 'FB.stan', data=dat, iter=2000, chains=2)
# stan(file = 'FB.stan', data=dat, iter=2000, chains=2)
gauss.end <- Sys.time()

# Clayton
clay.start <- Sys.time()
dat$copula <- 2
clayton <- stan(file = 'FB.stan', data=dat, iter=2000, chains=2)
# stan(file = 'FB.stan', data=dat, iter=2000, chains=2)
clay.end <- Sys.time()

# Gumbel
gumb.start <- Sys.time()
dat$copula <- 3
gumbel <- stan(file = 'FB.stan', data=dat, iter=2000, chains=2)
# stan(file = 'FB.stan', data=dat, iter=2000, chains=2)
gumb.end <- Sys.time()

# Frank
frank.start <- Sys.time()
dat$copula <- 4
frank <- stan(file = 'FB.stan', data=dat, iter=2000, chains=2)
# stan(file = 'FB.stan', data=dat, iter=2000, chains=2)
frank.end <- Sys.time()

# model selection
dic.gaussian <- dic1(extract(gaussian), dat, 1)  # 3620.975
dic.clayton <- dic1(extract(clayton), dat, 2)    # 3689.076
dic.gumbel <- dic1(extract(gumbel), dat, 3)      # 3607.364
dic.frank <- dic1(extract(frank), dat, 4)        # 3653.4
lpml.gaussian <- lpml(extract(gaussian), dat, 1) # -1811.103
lpml.clayton <- lpml(extract(clayton), dat, 2)   # -1845.921
lpml.gumbel <- lpml(extract(gumbel), dat, 3)     # -1804.144
lpml.frank <- lpml(extract(frank), dat, 4)       # -1827.48

# both in DIC and LPML, Gumbel is the best.

summary(gumbel, pars = c("b1[1]","b1[2]","b2[1]","b2[2]","sigma1", "phi"))

# save results
results <- list(dic.gaussian = dic.gaussian,
                dic.clayton = dic.clayton,
                dic.gumbel = dic.gumbel,
                dic.frank = dic.frank,
                lpml.gaussian = lpml.gaussian,
                lpml.clayton = lpml.clayton,
                lpml.gumbel = lpml.gumbel,
                lpml.frank = lpml.frank,
                fit = gumbel)
save(results, file = 'burn_analysis.RData')

### diagnostics
fit <- gumbel
rhats <- rhat(fit)
summary(rhats)  # > 1.1 is worrying
neffs <- neff_ratio(fit)
summary(neffs)  # < 0.1 is worrying
# needed for further diagnostics
lp <- log_posterior(fit)
np <- nuts_params(fit)
# traceplot
color_scheme_set("mix-brightblue-gray")
mcmc_nuts_divergence(np, lp)
mcmc_nuts_energy(np, merge_chains = FALSE)

pairs(fit, pars = c("b1[1]","b1[2]","b2[1]","b2[2]","sigma1", "phi", "lp__"))

summary(fit, pars = c("b1[1]","b1[2]","b2[1]","b2[2]","sigma1", "phi"))

mcmc_trace(as.array(fit), 
           pars = c("b1[1]","b1[2]","b2[1]","b2[2]","sigma1","phi")) 

theme_update(plot.title = element_text(hjust = 0.5))

t1 <- mcmc_trace(as.array(fit), pars = "b1[1]") + theme(axis.title.y=element_blank(), axis.title.x=element_blank(),legend.position = "none") +
                                labs(title = expression(beta[c0]))
t2 <- mcmc_trace(as.array(fit), pars = "b1[2]") + theme(axis.title.y=element_blank(), axis.title.x=element_blank(),legend.position = "none") +
                                labs(title = expression(beta[c1]))
t3 <- mcmc_trace(as.array(fit), pars = "b2[1]") + theme(axis.title.y=element_blank(), axis.title.x=element_blank(),legend.position = "none") +
                                labs(title = expression(beta[b0]))
t4 <- mcmc_trace(as.array(fit), pars = "b2[2]") + theme(axis.title.y=element_blank(), axis.title.x=element_blank(),legend.position = "none") +
                                labs(title = expression(beta[b1]))
t5 <- mcmc_trace(as.array(fit), pars = "sigma1")+ theme(axis.title.y=element_blank(), axis.title.x=element_blank(),legend.position = "none") +
                                labs(title = expression(sigma[c]))
t6 <- mcmc_trace(as.array(fit), pars = "phi")+ theme(axis.title.y=element_blank(), axis.title.x=element_blank(),legend.position = "none") +
                                labs(title = expression(phi))

d1 <- mcmc_hist(fit, pars = "b1[1]") + theme(axis.title.x=element_blank()) + labs(title = expression(beta[c0]))
d2 <- mcmc_hist(fit, pars = "b1[2]") + theme(axis.title.x=element_blank()) + labs(title = expression(beta[c1]))
d3 <- mcmc_hist(fit, pars = "b2[1]") + theme(axis.title.x=element_blank()) + labs(title = expression(beta[b0]))
d4 <- mcmc_hist(fit, pars = "b2[2]") + theme(axis.title.x=element_blank()) + labs(title = expression(beta[b1]))
d5 <- mcmc_hist(fit, pars = "sigma1")+ theme(axis.title.x=element_blank()) + labs(title = expression(sigma[c]))
d6 <- mcmc_hist(fit, pars = "phi")+ theme(axis.title.x=element_blank()) + labs(title = expression(phi))

grid.arrange(d1, d2, d3, d4, d5, d6, ncol=3)
grid.arrange(t1, t2, t3, t4, t5, t6, ncol=3)



# convert chains for betas to unstandardized scale, then do summary
samps <- extract(gumbel, 
                 pars = c("b1[1]","b1[2]","b2[1]","b2[2]","sigma1", "phi"))
samps$`b1[2]_us`  <- samps$`b1[2]`/attr((burn$age_std), "scaled:scale")
samps$`b1[1]_us`  <- samps$`b1[1]` - 
  samps$`b1[2]_us`*attr((burn$age_std), "scaled:center")
samps$`b2[2]_us`  <- samps$`b2[2]`/attr((burn$age_std), "scaled:scale")
samps$`b2[1]_us`  <- samps$`b2[1]` - 
  samps$`b2[2]_us`*attr((burn$age_std), "scaled:center")
# also get Kendall's tau
samps$tau <- numeric(length(samps$phi))
library(copula)
for(i in seq_along(samps$tau)) samps$tau[i] <- tau(gumbelCopula(samps$phi[i]))
library(coda)
samps_mat <- do.call(cbind, samps)
summary.burn <- summary(mcmc(samps_mat))
save(summary.burn, file = "burn_analysis_MCMC_summary.RData")
