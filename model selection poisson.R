##################################################################
# R code for Armillotta, Luati and Lupparelli 
# Observation driven models for discrete-valued time series
##################################################################

####################################################################
## Simulation study for model selection, Poisson distribution
###################################################################

# load required libraries

library(stats)
library(xtable)
library(numDeriv)
library(nloptr)
library(ggpubr)


###################################################################
# Function for simulations of model selection
###################################################################

simm_p<-function(n, s, delt){
  
  s_1p <- matrix(0, n, s)
  la_1p <- matrix(0, n, s)
  e_1p <- matrix(0, n, s)
  
  x_1pp <- matrix(0, n, s)
  la_1pp <- matrix(0, n, s)
  x_1pg <- matrix(0, n, s)
  la_1pg <- matrix(0, n, s)
  x_1pgl <- matrix(0, n, s)
  la_1pgl <- matrix(0, n, s)
  x_1pl <- matrix(0, n, s)
  la_1pl <- matrix(0, n, s)
  
  d_1pp <- matrix(nrow=s, ncol=3)
  d_1pg <- matrix(nrow=s, ncol=3)
  d_1pgl <- matrix(nrow=s, ncol=3)
  
  t_1pp <- matrix(nrow=s, ncol=3)
  t_1pg <- matrix(nrow=s, ncol=3)
  t_1pgl <- matrix(nrow=s, ncol=3)
  
  l_1p <- matrix(nrow=s, ncol=3)
  aic_1p <- matrix(nrow=s, ncol=3)
  bic_1p <- matrix(nrow=s, ncol=3)
  
  mean_1pp <- matrix(nrow=3, ncol=3)
  mean_1pg <- matrix(nrow=3, ncol=3)
  mean_1pgl <- matrix(nrow=3, ncol=3)
  
  icp <- matrix(nrow=2, ncol=3)
  
  #######################################################################
  # functions for generating the data
  #######################################################################
  
  g_logar.p <- function(del){
    
    y <- vector()
    x <-vector()
    mu <- vector()
    y[1] <- 0
    x[1] <- 0
    mu[1] <- 1
    
    be <- del[1]
    ph <- del[2]
    the <- del[3]
    
    
    for(i in 2:n) {
      x[i] <- be + ph*log(y[i-1]+1) + the*x[i-1]
      mu[i] <- exp(x[i])
      y[i] <- rpois(1, mu[i])
    }
    
    return(list(x=x, mu=mu, y=y))
  }
  
  
  
  g_garma.p <- function(del){
    
    y <- vector()
    x <-vector()
    mu <- vector()
    y[1] <- 0
    x[1] <- 0
    mu[1] <- 1
    
    be <- del[1]
    ph <- del[2]
    the <- del[3]
    
    y_ <- vector()
    gy <- vector()
    c <- 0.1
    
    y_[1] <- max(y[1], c)
    gy[1] <- log(y_[1])
    
    
    for(i in 2:n) {
      x[i] <- be + ph*gy[i-1] + the*(gy[i-1] - x[i-1])
      mu[i] <- exp(x[i])
      y[i] <- rpois(1, mu[i])
      y_[i] <- max(y[i], c)
      gy[i] <- log(y_[i])
    }
    
    return(list(x=x, mu=mu, y=y))
  }
  
  
  
  g_glarma.p <- function(del){
    
    y <- vector()
    x <-vector()
    mu <- vector()
    sd <- vector()
    y[1] <- 0
    x[1] <- 0
    mu[1] <- exp(x[1])
    sd[1] <- sqrt(mu[1])
    
    be <- del[1]
    ph <- del[2]
    the <- del[3]
    
    
    for(i in 2:n) {
      x[i] <- be + ph*x[i-1] + the*(y[i-1] - mu[i-1])/sd[i-1]
      mu[i] <- exp(x[i])
      sd[i] <- sqrt(mu[i])
      y[i] <- rpois(1, mu[i])
    }
    
    return(list(x=x, mu=mu, y=y))
  }
  
  ##########################################################################
  # log-likelihood and gradient functions
  #########################################################################
  
  
  # function for loglikelihood and (outer)gradient 
  # of Poisson log-linear autoregression
  
  logl.p <- function(delta){   
    
    g <- 0
    
    beta <- delta[1]   # parameters
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    eta <- vector()
    
    eta[1] <- 0      # starting value
    
    #log-likelihood
    
    for (i in 2:n){
      
      eta[i] <- beta + phi1*log(y[i-1]+1) + theta1*eta[i-1]  
      g <- g + y[i]*eta[i]-exp(eta[i])-lfactorial(y[i])       
      
    }
    
    g <- -g/n
    return(g)
    
  }
  
  gradd.p <- function(delta){  
    
    gr <- rep(0, 3)
    out <- matrix(0, 3, 3)
    
    beta <- delta[1]
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    eta <- vector()
    
    eta[1] <- 0
    
    for (i in 2:n){
      
      eta[i] <- beta + phi1*log(y[i-1]+1) + theta1*eta[i-1]
      
      # gradient
      
      gr[1] <- y[i]-exp(eta[i])
      gr[2] <- (y[i]-exp(eta[i]))*log(1+y[i-1])
      gr[3] <- (y[i]-exp(eta[i]))*eta[i-1]
      
      # outer product
      
      out <- out + gr%*%t(gr)
      
      
    }
    
    out <- out/n   # print outer product
    
    return(out)
    
  }
  
  # function for loglikelihood and (outer)gradient
  # of Poisson GARMA model
  
  garmal.p <- function(delta){  
    
    g <- 0
    
    beta <- delta[1]
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    y_ <- vector()
    gy <- vector()
    eta <- vector()
    c <- 0.1
    
    y_[1] <- max(y[1], c)
    gy[1] <- log(y_[1])
    eta[1] <- 0
    
    for (i in 2:n){
      
      y_[i] <- max(y[i], c)
      gy[i] <- log(y_[i])
      eta[i] <- beta + phi1*gy[i-1] + theta1*(gy[i-1] - eta[i-1])# + phi2*gy[i-2]
      g <- g + y[i]*eta[i]-exp(eta[i])-lfactorial(y[i])
      
    }
    
    g <- -g/n
    return(g)
    
  }
  
  garmagr.p <- function(delta){  
    
    gr <- rep(0, 3)
    out <- matrix(0, 3, 3)
    
    beta <- delta[1]
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    y_ <- vector()
    gy <- vector()
    eta <- vector()
    c <- 0.1
    
    y_[1] <- max(y[1], c)
    gy[1] <- log(y_[1])
    eta[1] <- 0
    
    for (i in 2:n){
      
      y_[i] <- max(y[i], c)
      gy[i] <- log(y_[i])
      
      eta[i] <- beta + phi1*gy[i-1] + theta1*(gy[i-1] - eta[i-1])
      
      gr[1] <- y[i]-exp(eta[i])
      gr[2] <- (y[i]-exp(eta[i]))*gy[i-1]
      gr[3] <- (y[i]-exp(eta[i]))*(gy[i-1] - eta[i-1])
      
      out <- out + gr%*%t(gr)
      
    }
    
    out <- out/n
    
    return(out)
    
  }
  
  # function for loglikelihood and (outer)gradient
  # of Poisson GLARMA model
  
  glarmal.p <- function(delta){ 
    
    g <- 0
    
    beta <- delta[1]
    gamma1 <- delta[2]
    theta1 <-  delta[3]
    
    lambda <- vector()
    eta <- vector()
    
    lambda[1] <- 1
    eta[1] <- 0
    
    for (i in 2:n){
      
      eta[i] <- beta + gamma1*eta[i-1] + theta1*(y[i-1] - lambda[i-1])/sqrt(lambda[i-1])# + gamma2*eta[i-2]
      lambda[i] <- exp(eta[i])
      
      g <- g + y[i]*eta[i]-exp(eta[i])-lfactorial(y[i])
      
    }
    
    g <- -g/n
    return(g)
    
  }
  
  glarmagr.p <- function(delta){
    
    gr <- rep(0, 3)
    out <- matrix(0, 3, 3)
    
    beta <- delta[1]
    gamma1 <- delta[2]
    theta1 <-  delta[3]
    
    lambda <- vector()
    eta <- vector()
    
    lambda[1] <- 1
    eta[1] <- 0
    
    for (i in 2:n){
      
      eta[i] <- beta + gamma1*eta[i-1] + theta1*(y[i-1] - lambda[i-1])/sqrt(lambda[i-1])# + gamma2*eta[i-2]
      lambda[i] <- exp(eta[i])
      
      gr[1] <- y[i]-exp(eta[i])
      gr[2] <- (y[i]-exp(eta[i]))*eta[i-1]
      gr[3] <- (y[i]-exp(eta[i]))*(y[i-1] - lambda[i-1])/sqrt(lambda[i-1])
      
      out <- out + gr%*%t(gr)
      
    }
    
    out <- out/n
    
    return(out)
    
  }
  
  # function for loglikelihood and (outer)gradient 
  # of Negative Binomial log-linear autoregression
  
  logl.nb <- function(delta){
    
    g <- 0
    
    beta <- delta[1]
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    eta <- vector()
    
    eta[1] <- 0
    
    for (i in 2:n){
      
      eta[i] <- beta + phi1*log(y[i-1]+1) + theta1*eta[i-1]
      g <- g + y[i]*eta[i]-(y[i]+v)*log(v+exp(eta[i]))+v*log(v)+
        lgamma(v+y[i])-lgamma(v)-lgamma(y[i]+1)
      
    }
    
    g <- -g/n
    return(g)
    
  }
  
  gradd.nb <- function(delta){
    
    gr <- rep(0, 3)
    out <- matrix(0, 3, 3)
    
    beta <- delta[1]
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    eta <- vector()
    
    eta[1] <- 0
    
    for (i in 2:n){
      
      eta[i] <- beta + phi1*log(y[i-1]+1) + theta1*eta[i-1]
      
      gr[1] <- y[i]-((y[i]+v)*exp(eta[i]))/(v+exp(eta[i]))
      gr[2] <- (y[i]-((y[i]+v)*exp(eta[i]))/(v+exp(eta[i])))*log(1+y[i-1])
      gr[3] <- (y[i]-((y[i]+v)*exp(eta[i]))/(v+exp(eta[i])))*eta[i-1]
      
      out <- out + gr%*%t(gr)
      
    }
    
    out <- out/n
    return(out)
    
  }
  
  # function for loglikelihood and (outer)gradient 
  # of Negative Binomial GARMA model
  
  garmal.nb <- function(delta){ 
    
    g <- 0
    
    beta <- delta[1]
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    y_ <- vector()
    gy <- vector()
    eta <- vector()
    c <- 0.1
    
    y_[1] <- max(y[1], c)
    gy[1] <- log(y_[1])
    eta[1] <- 0
    
    for (i in 2:n){
      
      y_[i] <- max(y[i], c)
      gy[i] <- log(y_[i])
      eta[i] <- beta + phi1*gy[i-1] + theta1*(gy[i-1] - eta[i-1])
      g <- g + y[i]*eta[i]-(y[i]+v)*log(v+exp(eta[i]))+v*log(v)+
        lgamma(v+y[i])-lgamma(v)-lgamma(y[i]+1)
      
    }
    
    g <- -g/n
    return(g)
    
  }
  
  garmagr.nb <- function(delta){
    
    gr <- rep(0, 3)
    out <- matrix(0, 3, 3)
    
    beta <- delta[1]
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    y_ <- vector()
    gy <- vector()
    eta <- vector()
    c <- 0.1
    
    y_[1] <- max(y[1], c)
    gy[1] <- log(y_[1])
    eta[1] <- 0
    
    for (i in 2:n){
      
      y_[i] <- max(y[i], c)
      gy[i] <- log(y_[i])
      
      eta[i] <- beta + phi1*gy[i-1] + theta1*(gy[i-1] - eta[i-1])
      
      gr[1] <- y[i]-((y[i]+v)*exp(eta[i]))/(v+exp(eta[i]))
      gr[2] <- (y[i]-((y[i]+v)*exp(eta[i]))/(v+exp(eta[i])))*gy[i-1]
      gr[3] <- (y[i]-((y[i]+v)*exp(eta[i]))/(v+exp(eta[i])))*(gy[i-1] - eta[i-1])
      
      out <- out + gr%*%t(gr)
      
    }
    
    out <- out/n
    return(out)
    
  }
  
  
  # function for loglikelihood and (outer)gradient 
  # of Negative Binomial GLARMA model
  
  glarmal.nb <- function(delta){
    
    g <- 0
    
    beta <- delta[1]
    gamma1 <- delta[2]
    theta1 <-  delta[3]
    
    lambda <- vector()
    eta <- vector()
    var <- vector()
    
    lambda[1] <- 1
    eta[1] <- 0
    var[1] <- 1+1/v
    
    for (i in 2:n){
      
      eta[i] <- beta + gamma1*eta[i-1] + theta1*(y[i-1] - lambda[i-1])/sqrt(var[i-1])
      lambda[i] <- exp(eta[i])
      var[i] <- lambda[i]*(1+lambda[i]/v)
      
      g <- g + y[i]*eta[i]-(y[i]+v)*log(v+exp(eta[i]))+v*log(v)+
        lgamma(v+y[i])-lgamma(v)-lgamma(y[i]+1)
      
    }
    
    g <- -g/n
    return(g)
    
  }
  
  glarmagr.nb <- function(delta){ 
    
    gr <- rep(0, 3)
    out <- matrix(0, 3, 3)
    
    beta <- delta[1]
    gamma1 <- delta[2]
    theta1 <-  delta[3]
    
    lambda <- vector()
    eta <- vector()
    var <- vector()
    
    lambda[1] <- 1
    eta[1] <- 0
    var[1] <- 1+1/v
    
    for (i in 2:n){
      
      eta[i] <- beta + gamma1*eta[i-1] + theta1*(y[i-1] - lambda[i-1])/sqrt(var[i-1])
      lambda[i] <- exp(eta[i])
      var[i] <- lambda[i]*(1+lambda[i]/v)
      
      gr[1] <- y[i]-((y[i]+v)*exp(eta[i]))/(v+exp(eta[i]))
      gr[2] <- (y[i]-((y[i]+v)*exp(eta[i]))/(v+exp(eta[i])))*eta[i-1]
      gr[3] <- (y[i]-((y[i]+v)*exp(eta[i]))/(v+exp(eta[i])))*(y[i-1] - lambda[i-1])/sqrt(var[i-1])
      
      out <- out + gr%*%t(gr)
      
    }
    
    out <- out/n
    return(out)
    
  }
  
  
  
  ##########################################################################
  # generate data
  ##########################################################################
  
  for(j in 1:s){
    
    # choose generetaing model
    
    model1p <- g_logar.p(delt)  
    #model1p <- g_garma.p(delt)
    #model1p <- g_glarma.p(delt)
    
    
    s_1p[,j] <- model1p$y
    la_1p[,j] <- model1p$mu
    e_1p[,j] <- model1p$x
    
    y <- s_1p[,j]
    
    
    # estimate models
    
    r <- 2
    
    d0 <- rep(0.01, 3)   # set initial values of the parameters
    
    # estimate log-linear AR
    
    d0[1] <- mean(y)
    ui <- rbind(c(0,-1,-1),c(0,1,0),c(0,0,1))
    ci <- c(-1,0,0)
    num_mle <- constrOptim(theta=d0, f=logl.p, ui=ui, ci=ci, method="Nelder-Mead",
                           control=list(trace=1, maxit=10000, reltol=sqrt(.Machine$double.eps)))
    num_mle$hessian <- optimHess(num_mle$par, logl.p)
    
    num_par <- num_mle$par        # store result of estimation
    num_l <- -num_mle$value       # value of maximum (quasi) log-likelihood
    
    num_J <- solve(num_mle$hessian)  #compute Hessian                              
    
    num_I <- gradd.p(num_par)          #compute conditional information matrix
    
    num_var <- num_J%*%num_I%*%num_J/n # sandiwhich estimator
    num_V <- diag(num_var)             # estimated variance
    num_V <- t(num_V)
    num_SE <- sqrt(num_V)              # estimated SE
    num_t <- num_par/num_SE
    num_p <- pt(abs(num_t), n-length(num_mle$par), lower.tail=FALSE)  # t test
    
    num_aic <- 2*length(num_mle$par) - 2*num_l      # compute IC     #AIC
    num_bic <- log(n)*length(num_mle$par) - 2*num_l                  #BIC
    
    x <- vector()
    x[1] <- 0
    lam <- vector()
    lam[1] <- exp(x[1])
    
    for (i in 2:n){
      
      x[i] <- num_par[1] + num_par[2]*log(y[i-1]+1) + num_par[3]*x[i-1]
      lam[i] <- exp(x[i])
      
    }
    
    # store results
    
    x_1pp[,j] <- x
    la_1pp[,j] <- lam
    d_1pp[j,] <- num_par
    t_1pp[j,] <- num_t
    l_1p[j, 1] <- num_l
    aic_1p[j, 1] <- num_aic
    bic_1p[j, 1] <- num_bic
    
    
    # estimate GARMA 
    
    d0[1] <- mean(y)
    num_mle <- optim(d0, garmal.p, method="BFGS", hessian=T)
    
    num_par <- num_mle$par        # store result of estimation
    num_l <- -num_mle$value       # value of maximum (quasi) log-likelihood
    
    num_J <- solve(num_mle$hessian)  #compute Hessian                              
    
    num_I <- garmagr.p(num_par)          #compute conditional information matrix
    
    num_var <- num_J%*%num_I%*%num_J/n # sandiwhich estimator
    num_V <- diag(num_var)             # estimated variance
    num_V <- t(num_V)
    num_SE <- sqrt(num_V)              # estimated SE
    num_t <- num_par/num_SE
    num_p <- pt(abs(num_t), n-length(num_mle$par), lower.tail=FALSE)  # t test
    
    num_aic <- 2*length(num_mle$par) - 2*num_l      # compute IC    #AIC
    num_bic <- log(n)*length(num_mle$par) - 2*num_l                 #BIC
    
    x <- vector()
    x[1] <- 0
    lam <- vector()
    lam[1] <- 1
    y_ <- vector()
    gy <- vector()
    c <- 0.1
    
    y_[1] <- max(y[1], c)
    gy[1] <- log(y_[1])
    
    
    for (i in 2:n){
      
      y_[i] <- max(y[i], c)
      gy[i] <- log(y_[i])
      x[i] <- num_par[1] + num_par[2]*gy[i-1] + num_par[3]*(gy[i-1] - x[i-1])
      lam[i] <- exp(x[i])
      
    }
    
    
    # store results
    
    x_1pg[,j] <- x
    la_1pg[,j] <- lam
    d_1pg[j,] <- num_par
    t_1pg[j,] <- num_t
    l_1p[j, 2] <- num_l
    aic_1p[j, 2] <- num_aic
    bic_1p[j, 2] <- num_bic
    
    # estimate GLARMA
    
    
    d0[1] <- mean(y)
    
    ui <- rbind(c(0,-1,-1))
    ci <- -1
    num_mle <- constrOptim(theta=d0, f=glarmal.p, ui=ui, ci=ci, method="Nelder-Mead",
                           control=list(trace=1, maxit=10000, reltol=sqrt(.Machine$double.eps)))
    num_mle$hessian <- optimHess(num_mle$par, glarmal.p)
    
    num_par <- num_mle$par        # store result of estimation
    num_l <- -num_mle$value       # value of maximum (quasi) log-likelihood
    
    num_aic <- 2*length(num_mle$par) - 2*num_l      # compute IC    #AIC
    num_bic <- log(n)*length(num_mle$par) - 2*num_l                 #BIC
    
    x <- vector()
    x[1] <- 0
    lam <- vector()
    lam[1] <- 1
    sv <- vector()
    sv[1] <- sqrt(lam[1])
    
    for (i in 2:n){
      
      x[i] <- num_par[1] + num_par[2]*x[i-1] + num_par[3]*(y[i-1] - lam[i-1])/sv[i-1]
      lam[i] <- exp(x[i])
      sv[i] <- sqrt(lam[i])
    }
    
    # store results
    
    x_1pgl[,j] <- x
    la_1pgl[,j] <- lam
    d_1pgl[j,] <- num_par
    t_1pgl[j,] <- 1
    l_1p[j, 3] <- num_l
    aic_1p[j, 3] <- num_aic
    bic_1p[j, 3] <- num_bic
    
    
  }
  
  
  #summary statistics
  
  mean_1pp[1,] <- colMeans(d_1pp)          #mean of parameters estimates
  mean_1pg[1,] <- colMeans(d_1pg)
  mean_1pgl[1,] <- colMeans(d_1pgl)
  
  mean_1pp[2,] <- apply(d_1pp, 2, sd)      #standard deviations
  mean_1pg[2,] <- apply(d_1pg, 2, sd)
  mean_1pgl[2,] <- apply(d_1pgl, 2, sd)
  
  freq<- function(x){
    sum(x>2.58)/s*100
  }
  
  mean_1pp[3,] <- apply(t_1pp, 2, FUN=freq)   # freq of rejecting tests
  mean_1pg[3,] <- apply(t_1pg, 2, FUN=freq)
  mean_1pgl[3,] <- apply(t_1pgl, 2, FUN=freq)
  
  mean <- cbind(mean_1pp, mean_1pg, mean_1pgl)
  
  aiccp <- aic_1p==apply(aic_1p, 1, min)   # freq of correct selection
  icp[1,] <- 100*apply(aiccp, 2, mean)
  biccp <- bic_1p==apply(bic_1p, 1, min)
  icp[2,] <- 100*apply(biccp, 2, mean)
  
  
  results <- list(eta=e_1p, lambda=la_1p, simulation=s_1p,
                  d_1pp=d_1pp, d_1pg=d_1pg, d_1pgl=d_1pgl,
                  mean=mean, icp=icp, logl=l_1p, aic=aic_1p, bic=bic_1p)
  
  return(results)
  
}

# Results

set.seed(12346)

#################log-ar

simmp_logar1 <- simm_p(n=250, s=1000, delt=c(0.2, 0.4, 0.3))
simmp_logar2 <- simm_p(n=500, s=1000, delt=c(0.2, 0.4, 0.3))
simmp_logar3 <- simm_p(n=1000, s=1000, delt=c(0.2, 0.4, 0.3))

round(simmp_logar1$mean, 4)
round(simmp_logar2$mean, 4)
round(simmp_logar3$mean, 4)
round(simmp_logar1$icp, 4)
round(simmp_logar2$icp, 4)
round(simmp_logar3$icp, 4)


##################garma 

simmp_garma1 <- simm_p(n=250, s=1000, delt=c(0.2, 0.4, 0.2))
simmp_garma2 <- simm_p(n=500, s=1000, delt=c(0.2, 0.4, 0.2))
simmp_garma3 <- simm_p(n=1000, s=1000, delt=c(0.2, 0.4, 0.2))

round(simmp_garma1$mean, 4)
round(simmp_garma2$mean, 4)
round(simmp_garma3$mean, 4)
round(simmp_garma1$icp, 4)
round(simmp_garma2$icp, 4)
round(simmp_garma3$icp, 4)


#################glarma

simmp_glarma1 <- simm_p(n=250, s=1000, delt=c(0.2, 0.3, 0.4))
simmp_glarma2 <- simm_p(n=500, s=1000, delt=c(0.2, 0.3, 0.4))
simmp_glarma3 <- simm_p(n=1000, s=1000, delt=c(0.2, 0.3, 0.4))

round(simmp_glarma1$mean, 4)
round(simmp_glarma2$mean, 4)
round(simmp_glarma3$mean, 4)
round(simmp_glarma1$icp, 4)
round(simmp_glarma2$icp, 4)
round(simmp_glarma3$icp, 4)

