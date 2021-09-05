## Cross-validation for penalized quantile regression
cv.qr <- function(X,Y,lambdas.qr,n){
  it <- length(lambdas.qr)
  err <- rep(0,it)
  folds <- floor(n/5)
  fm <- matrix(sample(1:n),nrow=folds)
  for(j in 1:it){
    for(i in folds){
      fit <- rq.fit.lasso(X[-fm[i,],], Y[-fm[i,]], tau = 0.5, lambda = lambdas.qr[j], beta = .9995, eps = 1e-06)
      err[j] <- err[j] +  sum(abs(Y[fm[i,]] - X[fm[i,],]%*% fit$coefficients))/n #+ lambdas.qr[j]*sum(abs(fit$coefficients))
    }
  }
  lambdas.qr[which.min(err)]
}

## Approximate inverse of the matrix t(X)%*%X/n 
nodewise <- function(x, p, n,lambdaj){
  C <- matrix(rep(0,p*p),nrow=p)
  for(j in 1:p){
    xj <- x[,j]
    xmj <- x[,-j]
    res2 <- glmnet(xmj,xj,lambda=lambdaj)
    gamaj <- res2$beta
    tau2j <- as.numeric((x[,j] %*% (x[,j] - predict(res2,x[,-j],s=lambdaj)))/n)
    Cj <- gamaj/tau2j;      
    if(j==1){ prva <- numeric(0);}else prva <- Cj[1:(j-1)];
    if(j==p){ druha <- numeric(0);} else druha <- Cj[j:(p-1)];
    C[j,] <- c(-prva,1/tau2j,-druha);
  }
  C
}

## Histogram plotting
his <- function(n,odhady,beta,vari,expr,j,nc,labe){
  
  hist(sqrt(n)*( ( odhady[,j]-beta[j] )/sqrt(vari[j,j]) ), 
       nclass=nc, 
       probability=TRUE,
       main = "", xlab=expr, ylab="Density",
       xlim=c(-3.2,3.2),ylim=c(0,0.5),col="grey",yaxt="n")
  
  lines(seq(-4,4,by=0.1), dnorm(seq(-4,4,by=0.1)),col="red")
  
}