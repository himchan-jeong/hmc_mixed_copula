library(mvtnorm)
library(methods)
library(Matrix)
# library(bindata)
library(copula)
library(compiler)
library(sn)

# coptype should be one of "normal", "frank", "clayton" or "gumbel"
gen_mixed_data <- function(Beta.true, N, rho, intercept, cov, xcor=0.25, coptype) {
  
  # Beta.true = c(B.true.n, B.true.b)
  # Assume equal length since same covariates
  stopifnot(length(Beta.true) %% 2 == 0)
  p <- length(Beta.true)/2
  B.true.n <- Beta.true[1 : p]
  B.true.b <- Beta.true[(p+1) : (2*p)]
  
  if(intercept) {
    W1 <- rmvnorm(N, mean=rep(0,p-1), sigma=corrmat(xcor, p-1, "AR1"))
    W2 <- rmvnorm(N, mean=rep(0,p-1), sigma=corrmat(xcor, p-1, "AR1"))
    # w3 <- rbinom(N, 1, 0.3)
    w3 <- rnorm(N)
    if(cov=="separate") {
      # W1 <- scale(W1);W2 <- scale(W2)
      X <- cbind(1, W1)
      Z <- cbind(1, W2)
    } else if (cov=="same") {
      X <- Z <- cbind(1, W2)
    } else if (cov=="shared") {
      X <- Z <- cbind(1, w3, W1[,1])
      X <- cbind(X, W1[,2:(p-2)]) # X = [1 w3 w1(1) w1(2)...w1(p-2)] 
      Z <- cbind(Z, W2[,1:(p-3)]) # Z = [1 w3 w1(1) w2(1)...w2(p-3)]
    } else stop("Error in gen_mixed_data: Bad value for cov specified")
  } else {
    X <- matrix(runif(N*p), nrow=N)
    # X <- scale(X)
    if(sepcov) {
      Z <- matrix(rnorm(N*p), nrow=N)
      # Z <- scale(Z)
    } else Z <- X
  }

  # Method: Generate correlated uniform samples from Gaussian copula
  #   Then apply normal to one, logistic to other. Cutoff logistic at 0
  y <- numeric(N*2)
  q <- numeric(N)
  cop <- switch(coptype, 
                normal = normalCopula(rho),
                frank = frankCopula(param=rho),
                clayton = claytonCopula(param=rho),
                gumbel = gumbelCopula(param=rho))
  for (i in 1:N) {
    copdata <- rCopula(1, cop)
    y[2*i - 1] <- qnorm(copdata[1], mean = crossprod(X[i,], B.true.n))
    # y[2*i - 1] <- qt(copdata[1], df=1) + crossprod(X[i,], B.true.n)
    # y[2*i - 1] <- qsn(copdata[1], xi=crossprod(X[i,], B.true.n), alpha=4)
    # y[2*i] <- ( (qlogis(copdata[2], location=crossprod(Z[i,], B.true.b)) > 0) * 2) - 1 # +1 -1 coding
    q[i]   <- qlogis(copdata[2], location = crossprod(Z[i,], B.true.b)) #0 1 coding
    y[2*i] <- (q[i] > 0)*1 #0 1 coding
    #y[2*i] <- qnorm(copdata[2], mean = crossprod(X[i,], B.true.b))
  }
  # Return generated data as list
  clustid <- rep(1:N, each=2)
  l <- list(clustid=clustid, y=y, q=q, X=X, Z=Z, N=N, m=2)
}

