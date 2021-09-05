
source('requirements.R')
source('utils.R')

n <- 200
p <- 10

diag.value <- 1
rho <- 0.3
lower <- rbind(rep(0,p), (diag(p))[-p,])
upper <- t(lower)
Theta0 <- diag.value * diag(p) + rho * (upper + lower)
Sigma0 <- solve(Theta0)

sqrSigma <- sqrtm(Sigma0)

# the regression vector
beta <- rep(0,p)
beta[1:3] <- c(1, 1, 1)

iter <- 500

length <- coverage <- rep(0,p)
estVar <- est <- matrix(0, nrow = iter, ncol = p)
sig <- 1

lambdas.qr <- seq(10,60,2)

pb <- txtProgressBar(0, iter, style = 3)

set.seed(1)

for(i in 1:iter){
  setTxtProgressBar(pb, i)
  # Generate data 
  nu <- 3
  # errors with t-distribution with 3 degrees of freedom
  eps <- rt(n,nu) / sqrt(nu/(nu-2)) 
  X <- mvrnorm(n, rep(0,p), Sigma0)
  Y <- X %*% beta + eps
  
  ## Compute de-biased LAD estimator 
  Theta <- nodewise(X, p, n, sqrt(log(p)/n))
  
  # LAD estimator
  lam <- lambdas.qr[which.min(cv.qr(X,Y,lambdas.qr,n))]
  
  f <- rq.fit.lasso(X, Y, tau = 0.5,beta = .9995, eps = 1e-06)
  beta.hat <- f$coefficients
  lamZ <-  t(X) %*% (sign(Y-X%*% beta.hat))  /n

  f0 <- dt(0,nu)*sqrt(nu/(nu-2))
  beta.ds <- beta.hat +  1/2/ f0 *Theta  %*% (lamZ)  
  
  # Compute confidence intervals 
  sigma <-  1/4/f0^2 *(diag(Theta)) # diag(Theta)
  rad <- qnorm(0.975)*sqrt(sigma) / sqrt(n)
  LE <- beta.ds - rad 
  UE <- beta.ds + rad
  
  # Compute average coverage and length
  length <-  length + 2*rad
  coverage <- coverage + as.numeric((beta <= UE)*(beta >= LE)==1)
  
  est[i,] <- beta.ds
  estVar[i,] <- sigma
  
}

pb.close()

# Print results
average_lengths <- length/iter
average_coverages <- coverage/iter
round(average_lengths,3) 
round(average_coverages,3) 

# Plot histograms
expr <- c(expression(beta [1]),expression(beta [2]),expression(beta [3]),expression(beta [4]))
par(mar=c(1,1,1,1)*2)
par(mfrow=c(2,2))

for(j in 1:4){
  his(n,odhady,beta,Theta/4/f0^2,expr[j],j,12,"De-biased LAD")
}

