source('requirements.R')
source('utils.R')

n <- 200
p <- 10

iter <- 500
length <- coverage <- rep(0, p)
estVar <- est <- matrix(0, nrow = iter, ncol = p)
lambdas.qr <- seq(10, 60, 2)
pars <- pars_init(n, p)
Sigma0 <- pars$Sigma0
beta <- pars$beta


pb <- txtProgressBar(0, iter, style = 3)

set.seed(1)

for(i in 1:iter){
  setTxtProgressBar(pb, i)
  
  # Generate data:
  # Errors with t-distribution with 3 degrees of freedom
  nu <- 3
  eps <- rt(n,nu) / sqrt(nu/(nu-2)) 
  X <- mvrnorm(n, rep(0,p), Sigma0)
  Y <- X %*% beta + eps
  
  ## Compute de-biased LAD estimator 
  out <- debiasedLAD(n, p, X, Y, lambdas.qr, nu)
  UE = out$UE
  LE = out$LE
  beta.ds = out$beta.ds
  sigma = out$sigma
  Theta = out$Theta
  f0 = out$f0
  
  # Compute average coverage and length
  length <-  length + 2 * (UE - LE)
  coverage <- coverage + as.numeric((beta <= UE)*(beta >= LE) == 1)
  
  est[i,] <- beta.ds
  estVar[i,] <- sigma
  
}

pb.close()

# Print results
round(length/iter,3) 
round(coverage/iter,3) 

# Plot histograms
expr <- c(expression(beta [1]),expression(beta [2]),expression(beta [3]),expression(beta [4]))
par(mar=c(1,1,1,1)*2)
par(mfrow=c(2,2))

for(j in 1:4){
  his(n, est, beta, Theta/4/f0^2, expr[j], j, 12, "De-biased LAD")
}

