##################################################################
# R code for Armillotta, Luati and Lupparelli 
# Observation driven models for discrete-valued time series
##################################################################

# load required libraries

library(stats)
library(xtable)
library(pracma)

##################################################################
# Simulation MLE for Bernoulli GLARMA model
##################################################################

set.seed(12346)  
n1 <- 200         #sample size
n2 <- 500 
n3 <- 1000
p <- 1           # select lag 1 model
q <- 1
s <- 1000        # simulation iterations

delt1  <- c(0.5, -0.4, 0.8)  #true values parameters
delt3 <- c(0.5, 0.4, 0.2)  
delt5 <- c(0.5, 0.4, 1.2)


sim<-function(n,p,q,s,delt){
  
  r <- max(p,q)
  k <- p+q+1             #number of parameters
  d0 <- rep(0, k)        #initial values for estimation
  deltaB <- matrix(nrow=s, ncol=k)
  stdB <- matrix(nrow=s, ncol=k)
  biaB <- matrix(nrow=s, ncol=k)
  vb <- matrix(nrow=s, ncol=k)
  logB <- matrix(nrow=s, ncol=1)
  aicB <- matrix(nrow=s, ncol=1)
  bicB <- matrix(nrow=s, ncol=1)
  meanB <- vector()
  biasB <- vector()
  sdB <- vector()
  mseB <- vector()
  inf <- vector()
  sup <- vector()
  simulation <- matrix(0, n, s)
  prob <- matrix(0, n, s)
  
  # loglikelihood for GLARMA model with Bernoulli distribution
  
  sloglgl11 <- function(delta){
    
    g <- 0
    
    beta <- delta[1]
    gamma1 <- delta[2]
    theta1 <-  delta[3]
    
    for (i in 2:n){
      g <- g + y[i]*(beta + gamma1*eta[i-1] + theta1*(y[i-1] - lam[i-1])/var[i-1]) -
        log(1 + exp(beta + gamma1*eta[i-1] + theta1*(y[i-1] - lam[i-1])/var[i-1])) 
    }
    g <- -1 * g
    return(g)
  }
  
  slogl <- sloglgl11
  
  bet <- delt[1]
  gamm1 <- delt[2]
  thet1 <- delt[3]
  
  for(j in 1:s){
    
    #generate samples
    
    y <- vector()
    x <-vector()
    eta <- vector()
    lam <- vector()
    var <- vector()
    pr <- vector()
    c <- 0.1
    
    y[1] <- 1
    lam[1]<- 0.5
    eta[1] <- 0
    var[1] <- 1
    pr[1] <- 0.5
    
    for(i in 2:n) {
      eta[i] <- bet + gamm1*eta[i-1] + thet1*(y[i-1] - lam[i-1])/var[i-1]
      pr[i] <- exp(eta[i])/(1 + exp(eta[i]))
      lam[i] <- pr[i]
      var[i] <- 1
      y[i] <- rbinom(1, 1, pr[i])
      
    }
    
    simulation[,j] <- y
    prob[,j] <- pr
    
    # MLE estimation
    
    d0[1] <- mean(y)
    mleb <- optim(d0, slogl, method="BFGS", hessian=TRUE)
    varb <- solve(mleb$hessian)
    VB <- diag(varb)
    VB <- t(VB)
    par <- mleb$par
    l <- -mleb$value
    
    #storing results
    
    deltaB[j,] <- par
    
    biaB[j,] <- (par - delt)^2
    
    vb[j,] <- VB
    
    logB[j] <- l
    aicB[j] <- 2*k - 2*logB[j] 
    bicB[j] <- log(n)*k - 2*logB[j]
    
  }
  
  stdB <- sqrt(vb)
  
  #summary statistics
  
  meanB <- colMeans(deltaB)
  
  biasB <- meanB - delt
  
  sdB <- apply(deltaB, 2, sd)
  
  inf <- meanB - 1.96*sdB/sqrt(s)
  sup <- meanB + 1.96*sdB/sqrt(s)
  
  ciB <- cbind(inf,sup)
  
  mseB <- colMeans(biaB)
  
  #mean information criteria
  
  aic_ <- mean(aicB)
  bic_ <- mean(bicB)
  
  #Kolmogorov-Smirnov test (asymptotic normality)
  
  zB1<-(deltaB[,1]-meanB[1])/sdB[1]
  zB2<-(deltaB[,2]-meanB[2])/sdB[2]
  zB3<-(deltaB[,3]-meanB[3])/sdB[3]
  ks1<-ks.test(zB1,pnorm)
  ks2<-ks.test(zB2,pnorm)
  ks3<-ks.test(zB3,pnorm)
  ks<-c(ks1$p.value, ks2$p.value, ks3$p.value)
  
  #results
  
  resultsB <- list(prob=prob, simulation=simulation, par.simu=deltaB, true=delt, mean=meanB, sd=sdB, conf_low=inf, conf_up=sup, bias=biasB,
                   mse=mseB, aic=aic_, bic=bic_, ks_test=ks)
  return(resultsB)
  
}

GLsim11<-sim(n1,p,q,s,delt1)
GLsim21<-sim(n2,p,q,s,delt1)
GLsim31<-sim(n3,p,q,s,delt1)

GLsim12<-sim(n1,p,q,s,delt3)
GLsim22<-sim(n2,p,q,s,delt3)
GLsim32<-sim(n3,p,q,s,delt3)

GLsim13<-sim(n1,p,q,s,delt5)
GLsim23<-sim(n2,p,q,s,delt5)
GLsim33<-sim(n3,p,q,s,delt5)




##################################################################
# Simulation QMLE for Poisson GARMA model (data from geometric)
##################################################################

set.seed(12346)  
n1 <- 200         #sample size
n2 <- 500 
n3 <- 1000
p <- 1           #lag 1 model
q <- 1
s <- 1000        # simulation iterations

delt1  <- c(0.5, -0.4, 0.8)  #true values parameters
delt3 <- c(0.5, 0.4, 0.2) 
delt5 <- c(0.5, 0.4, 1.2)


sim<-function(n,p,q,s,delt){
  
  r <- max(p,q)
  k <- p+q+1             #number of parameters
  d0 <- rep(0, k)        #initial values for estimation
  deltaB <- matrix(nrow=s, ncol=k)
  stdB <- matrix(nrow=s, ncol=k)
  biaB <- matrix(nrow=s, ncol=k)
  vb <- matrix(nrow=s, ncol=k)
  logB <- matrix(nrow=s, ncol=1)
  aicB <- matrix(nrow=s, ncol=1)
  bicB <- matrix(nrow=s, ncol=1)
  meanB <- vector()
  biasB <- vector()
  sdB <- vector()
  mseB <- vector()
  inf <- vector()
  sup <- vector()
  simulation <- matrix(0, n, s)
  prob <- matrix(0, n, s)
  
  # loglikelihood for GARMA model with Poisson distribution
  
  sloglg11p <- function(delta){
    
    g <- 0
    
    beta <- delta[1]
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    for (i in 2:n){
      g <- g + y[i]*(beta + phi1*gy[i-1] + theta1*(gy[i-1] - eta[i-1])) -
        exp(beta + phi1*gy[i-1] + theta1*(gy[i-1] - eta[i-1])) 
    }
    g <- -1 * g
    return(g)
  }
  slogl <- sloglg11p
  
  bet <- delt[1]
  ph1 <- delt[2]
  thet1 <- delt[3]
  
  for(j in 1:s){
    
    #generate samples (form Geometric ditribution)
    
    y <- vector()
    y_ <- vector()
    gy <- vector()
    x <-vector()
    eta <- vector()
    lam <- vector()
    pr <- vector()
    c <- 0.1
    
    y[1] <- 1
    y_[1] <- 1-c
    gy[1] <- 1-c
    eta[1] <- 0
    pr[1] <- 0.5
    
    for(i in 2:n) {
      eta[i] <- bet + ph1*gy[i-1] + thet1*(gy[i-1] - eta[i-1])
      lam[i] <- exp(eta[i])
      pr[i] <- 1/(1 + lam[i])
      y[i] <- rgeom(1, pr[i])
      y_[i] <- max(y[i], c)
      gy[i] <- log(y_[i])
    }
    
    simulation[,j] <- y
    prob[,j] <- pr
    
    # QMLE estimation
    
    d0[1] <- mean(y)
    mleb <- optim(d0, slogl, method="BFGS", hessian=TRUE)
    varb <- solve(mleb$hessian)
    VB <- diag(varb)
    VB <- t(VB)
    par <- mleb$par
    l <- -mleb$value
    
    #storing results
    
    deltaB[j,] <- par
    
    biaB[j,] <- (par - delt)^2
    
    vb[j,] <- VB
    
    logB[j] <- l
    aicB[j] <- 2*k - 2*logB[j] 
    bicB[j] <- log(n)*k - 2*logB[j]
    
  }
  
  stdB <- sqrt(vb)
  
  #summary statistics
  
  meanB <- colMeans(deltaB)
  
  biasB <- meanB - delt
  
  sdB <- apply(deltaB, 2, sd)
  
  inf <- meanB - 1.96*sdB/sqrt(s)
  sup <- meanB + 1.96*sdB/sqrt(s)
  
  ciB <- cbind(inf,sup)
  
  mseB <- colMeans(biaB)
  
  #mean information criteria
  
  aic_ <- mean(aicB)
  bic_ <- mean(bicB)
  
  #Kolmogorov-Smirnov test (asymptotic normality)
  
  zB1<-(deltaB[,1]-meanB[1])/sdB[1]
  zB2<-(deltaB[,2]-meanB[2])/sdB[2]
  zB3<-(deltaB[,3]-meanB[3])/sdB[3]
  ks1<-ks.test(zB1,pnorm)
  ks2<-ks.test(zB2,pnorm)
  ks3<-ks.test(zB3,pnorm)
  ks<-c(ks1$p.value, ks2$p.value, ks3$p.value)
  
  #results
  
  resultsB <- list(prob=prob, simulation=simulation, par.simu=deltaB, true=delt, mean=meanB, sd=sdB, conf_low=inf, conf_up=sup, bias=biasB,
                   mse=mseB, aic=aic_, bic=bic_, ks_test=ks)
  return(resultsB)
  
}

GPsim11<-sim(n1,p,q,s,delt1)
GPsim21<-sim(n2,p,q,s,delt1)
GPsim31<-sim(n3,p,q,s,delt1)

GPsim12<-sim(n1,p,q,s,delt3)
GPsim22<-sim(n2,p,q,s,delt3)
GPsim32<-sim(n3,p,q,s,delt3)

GPsim13<-sim(n1,p,q,s,delt5)
GPsim23<-sim(n2,p,q,s,delt5)
GPsim33<-sim(n3,p,q,s,delt5)



##################################################################
# Simulation QMLE for Poisson Log-ar model (data from geometric)
##################################################################

set.seed(12346)  
n1 <- 200         #sample size
n2 <- 500 
n3 <- 1000
p <- 1           #lag 1 model
q <- 1
s <- 1000        # simulation iterations
delt1  <- c(0.5, -0.4, 0.8)  #true values parameters
delt3 <- c(0.5, 0.4, 0.2)  


sim<-function(n,p,q,s,delt){
  
  r <- max(p,q)
  k <- p+q+1             #numer of parameters
  d0 <- rep(0, k)        #initial values estimation
  deltaB <- matrix(nrow=s, ncol=k)
  stdB <- matrix(nrow=s, ncol=k)
  biaB <- matrix(nrow=s, ncol=k)
  vb <- matrix(nrow=s, ncol=k)
  logB <- matrix(nrow=s, ncol=1)
  aicB <- matrix(nrow=s, ncol=1)
  bicB <- matrix(nrow=s, ncol=1)
  meanB <- vector()
  biasB <- vector()
  sdB <- vector()
  mseB <- vector()
  inf <- vector()
  sup <- vector()
  simulation <- matrix(0, n, s)
  prob <- matrix(0, n, s)
  
  
  # loglikelihood for log-ar model with Poisson distribution
  
  sloglpar <- function(delta){
    
    g <- 0
    
    beta <- delta[1]
    phi1 <- delta[2]
    gamma1 <-  delta[3]
    
    for (i in 2:n){
      g <- g + y[i]*(beta + phi1*hy[i-1] + gamma1*eta[i-1]) -
        exp(beta + phi1*hy[i-1] + gamma1*eta[i-1]) 
    }
    g <- -1 * g
    return(g)
  }
  slogl <- sloglpar
  
  bet <- delt[1]
  ph1 <- delt[2]
  gamm1 <- delt[3]
  
  for(j in 1:s){
    
    #generate samples (from geometric ditribution)
    
    y <- vector()
    y_ <- vector()
    hy <- vector()
    x <-vector()
    eta <- vector()
    lam <- vector()
    pr <- vector()
    
    y[1] <- 1
    hy[1] <- log(2)
    eta[1] <- 0
    pr[1] <- 0.5
    
    for(i in 2:n) {
      eta[i] <- bet + ph1*hy[i-1] + gamm1*eta[i-1]
      lam[i] <- exp(eta[i])
      pr[i] <- 1/(1 + lam[i])
      y[i] <- rgeom(1, pr[i])
      hy[i] <- log(y[i]+1)
    }
    
    simulation[,j] <- y
    prob[,j] <- pr
    
    # QMLE estimation
    
    d0[1] <- mean(y)
    mleb <- optim(d0, slogl, method="BFGS", hessian=TRUE)
    varb <- solve(mleb$hessian)
    VB <- diag(varb)
    VB <- t(VB)
    par <- mleb$par
    l <- -mleb$value
    
    #storing results
    
    deltaB[j,] <- par
    
    biaB[j,] <- (par - delt)^2
    
    vb[j,] <- VB
    
    logB[j] <- l
    aicB[j] <- 2*k - 2*logB[j] 
    bicB[j] <- log(n)*k - 2*logB[j]
    
  }
  
  stdB <- sqrt(vb)
  
  #summary statistics
  
  meanB <- colMeans(deltaB)
  
  biasB <- meanB - delt
  
  sdB <- apply(deltaB, 2, sd)
  
  inf <- meanB - 1.96*sdB/sqrt(s)
  sup <- meanB + 1.96*sdB/sqrt(s)
  
  ciB <- cbind(inf,sup)
  
  mseB <- colMeans(biaB)
  
  #mean information criteria
  
  aic_ <- mean(aicB)
  bic_ <- mean(bicB)
  
  #Kolmogorov-Smirnov test (asymptotic normality)
  
  zB1<-(deltaB[,1]-meanB[1])/sdB[1]
  zB2<-(deltaB[,2]-meanB[2])/sdB[2]
  zB3<-(deltaB[,3]-meanB[3])/sdB[3]
  ks1<-ks.test(zB1,pnorm)
  ks2<-ks.test(zB2,pnorm)
  ks3<-ks.test(zB3,pnorm)
  ks<-c(ks1$p.value, ks2$p.value, ks3$p.value)
  
  #results
  
  resultsB <- list(prob=prob, simulation=simulation, par.simu=deltaB, true=delt, mean=meanB, sd=sdB, conf_low=inf, conf_up=sup, bias=biasB,
                   mse=mseB, aic=aic_, bic=bic_, ks_test=ks)
  return(resultsB)
  
}

Psim11<-sim(n1,p,q,s,delt1)
Psim21<-sim(n2,p,q,s,delt1)
Psim31<-sim(n3,p,q,s,delt1)

Psim12<-sim(n1,p,q,s,delt3)
Psim22<-sim(n2,p,q,s,delt3)
Psim32<-sim(n3,p,q,s,delt3)

