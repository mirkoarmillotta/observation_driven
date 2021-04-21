##################################################################
# R code for Armillotta, Luati and Lupparelli 
# Observation driven models for discrete-valued time series
##################################################################

####################################################################
## Simulation study for model selection, Binomial distribution
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

simm<-function(a, n, s, delt){
  
  s_1b <- matrix(0, n, s)
  p_1b <- matrix(0, n, s)
  e_1b <- matrix(0, n, s)
  
  x_1bb <- matrix(0, n, s)
  pro_1bb <- matrix(0, n, s)
  x_1bg <- matrix(0, n, s)
  pro_1bg <- matrix(0, n, s)
  x_1bgl <- matrix(0, n, s)
  pro_1bgl <- matrix(0, n, s)
  x_1bl <- matrix(0, n, s)
  pro_1bl <- matrix(0, n, s)
  
  d_1bb <- matrix(nrow=s, ncol=3)
  d_1bg <- matrix(nrow=s, ncol=3)
  d_1bgl <- matrix(nrow=s, ncol=3)
  
  t_1bb <- matrix(nrow=s, ncol=3)
  t_1bg <- matrix(nrow=s, ncol=3)
  t_1bgl <- matrix(nrow=s, ncol=3)
  
  l_1b <- matrix(nrow=s, ncol=3)
  aic_1b <- matrix(nrow=s, ncol=3)
  bic_1b <- matrix(nrow=s, ncol=3)
  
  mean_1bb <- matrix(nrow=3, ncol=3)
  mean_1bg <- matrix(nrow=3, ncol=3)
  mean_1bgl <- matrix(nrow=3, ncol=3)
  
  ic <- matrix(nrow=2, ncol=3)
  
  #######################################################################
  # functions for generating the data
  #######################################################################
  
  g_barma.b <- function(a, del){
    
    y <- vector()
    x <-vector()
    pr <- vector()
    mu <- vector()
    y[1] <- 0
    x[1] <- 0
    pr[1] <- exp(x[1])/(1 + exp(x[1]))
    mu[1] <- a*pr[1]
    
    be <- del[1]
    ph <- del[2]
    the <- del[3]
    
    
    for(i in 2:n) {
      x[i] <- be + ph*y[i-1] + the*(y[i-1] - mu[i-1])
      pr[i] <- exp(x[i])/(1 + exp(x[i]))
      mu[i] <- a*pr[i]
      y[i] <- rbinom(1, a, pr[i])
    }
    
    return(list(x=x, pr=pr, mu=mu, y=y))
  }
  
  
  g_garma.b <- function(a, del){
    
    y <- vector()
    x <-vector()
    pr <- vector()
    y[1] <- 0
    x[1] <- 0
    pr[1] <- exp(x[1])/(1 + exp(x[1]))
    
    be <- del[1]
    ph <- del[2]
    the <- del[3]
    
    y_ <- vector()
    gy <- vector()
    c <- 0.1
    
    y_[1] <- min(max(y[1], c), a-c)
    gy[1] <- log(y_[1]/(a-y_[1]))
    
    
    for(i in 2:n) {
      x[i] <- be + ph*gy[i-1] + the*(gy[i-1] - x[i-1])
      pr[i] <- exp(x[i])/(1 + exp(x[i]))
      y[i] <- rbinom(1, a, pr[i])
      y_[i] <- min(max(y[i], c), a-c)
      gy[i] <- log(y_[i]/(a-y_[i]))
    }
    
    return(list(x=x, pr=pr, y=y))
  }
  
  g_glarma.b <- function(a, del){
    
    y <- vector()
    x <-vector()
    pr <- vector()
    mu <- vector()
    sd <- vector()
    y[1] <- 0
    x[1] <- 0
    pr[1] <- exp(x[1])/(1 + exp(x[1]))
    mu[1] <- a*pr[1]
    sd[1] <- sqrt(a*pr[1]*(1-pr[1]))
    
    be <- del[1]
    ph <- del[2]
    the <- del[3]
    
    
    for(i in 2:n) {
      x[i] <- be + ph*x[i-1] + the*(y[i-1] - mu[i-1])/sd[i-1]
      pr[i] <- exp(x[i])/(1 + exp(x[i]))
      mu[i] <- a*pr[i]
      sd[i] <- sqrt(a*pr[i]*(1-pr[i]))
      y[i] <- rbinom(1, a, pr[i])
    }
    
    return(list(x=x, pr=pr, mu=mu, y=y))
  }
  
  
  ##########################################################################
  # log-likelihood and gradient functions
  #########################################################################
  
  
  # loglikelihood and gradient for BARMA
  
  barl.b <- function(delta){
    
    g <- 0
    
    eta <- vector()
    prob <- vector()
    m <- vector()
    
    eta[1] <- 0
    prob[1] <- exp(eta[1])/(1 + exp(eta[1]))
    m[1] <- a*prob[1]
    
    beta <- delta[1]
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    for (i in 2:n){
      
      eta[i] <- beta + phi1*y[i-1] + theta1*(y[i-1] - m[i-1])
      prob[i] <- exp(eta[i])/(1 + exp(eta[i]))
      m[i] <- a*prob[i]
      g <- g + log(choose(a,y[i])) + y[i]*log(prob[i]) + (a-y[i])*log(1-prob[i])
    }
    g <- -1 * g
    return(g)
  }
  
  
  bargr.b <- function(delta){  # gradient
    
    gr <- rep(0, 3)
    out <- matrix(0, 3, 3)
    
    
    eta <- vector()
    m <- vector()
    prob <- vector()
    
    eta[1] <- 0
    prob[1] <- exp(eta[1])/(1 + exp(eta[1]))
    m[1] <- a*prob[1]
    
    beta <- delta[1]
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    for (i in 2:n){
      
      eta[i] <- beta + phi1*y[i-1] + theta1*(y[i-1] - m[i-1])
      prob[i] <- exp(eta[i])/(1 + exp(eta[i]))
      m[i] <- a*prob[i]
      
      gr[1] <- y[i]-m[i]
      gr[2] <- (y[i]-m[i])*y[i-1]
      gr[3] <- (y[i]-m[i])*(y[i-1] - m[i-1])
      
      out <- out + gr%*%t(gr)
      
    }
    
    out <- out/n
    
    return(out)
    
  }
  
  
  # loglikelihood and gradient for Binomial GARMA
  
  garmal.b <- function(delta){
    
    g <- 0
    
    beta <- delta[1]
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    eta <- vector()
    m <- vector()
    prob <- vector()
    
    eta[1] <- 0
    prob[1] <- exp(eta[1])/(1 + exp(eta[1]))
    m[1] <- a*prob[1]
    
    y_ <- vector()
    gy <- vector()
    c <- 0.1
    
    y_[1] <- min(max(y[1], c), a-c)
    gy[1] <- log(y_[1]/(a-y_[1]))
    
    
    for(i in 2:n) {
      y_[i] <- min(max(y[i], c), a-c)
      gy[i] <- log(y_[i]/(a-y_[i]))
      eta[i] <- beta + phi1*gy[i-1] + theta1*(gy[i-1] - eta[i-1])
      prob[i] <- exp(eta[i])/(1 + exp(eta[i]))
      g <- g + log(choose(a,y[i])) + y[i]*log(prob[i]) + (a-y[i])*log(1-prob[i])
    }
    g <- -1 * g
    return(g)
  }
  
  
  garmagr.b <- function(delta){  # gradient
    
    gr <- rep(0, 3)
    out <- matrix(0, 3, 3)
    
    
    eta <- vector()
    m <- vector()
    prob <- vector()
    
    eta[1] <- 0
    prob[1] <- exp(eta[1])/(1 + exp(eta[1]))
    m[1] <- a*prob[1]
    
    y_ <- vector()
    gy <- vector()
    c <- 0.1
    
    y_[1] <- min(max(y[1], c), a-c)
    gy[1] <- log(y_[1]/(a-y_[1]))
    
    beta <- delta[1]
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    for (i in 2:n){
      
      y_[i] <- min(max(y[i], c), a-c)
      gy[i] <- log(y_[i]/(a-y_[i]))
      eta[i] <- beta + phi1*gy[i-1] + theta1*(gy[i-1] - eta[i-1])
      prob[i] <- exp(eta[i])/(1 + exp(eta[i]))
      m[i] <- a*prob[i]
      
      gr[1] <- y[i]-m[i]
      gr[2] <- (y[i]-m[i])*gy[i-1]
      gr[3] <- (y[i]-m[i])*(gy[i-1] - eta[i-1])
      
      out <- out + gr%*%t(gr)
      
    }
    
    out <- out/n
    
    return(out)
    
  }
  
  
  # loglikelihood and gradient for Binomial GLARMA
  
  glarmal.b <- function(delta){
    
    g <- 0
    
    eta <- vector()
    prob <- vector()
    m <- vector()
    sdv <- vector()
    
    eta[1] <- 0
    prob[1] <- exp(eta[1])/(1 + exp(eta[1]))
    m[1] <- a*prob[1]
    sdv[1] <- sqrt(a*prob[1]*(1-prob[1]))
    
    beta <- delta[1]
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    for (i in 2:n){
      
      eta[i] <- beta + phi1*eta[i-1] + theta1*(y[i-1] - m[i-1])/sdv[i-1]
      prob[i] <- exp(eta[i])/(1 + exp(eta[i]))
      m[i] <- a*prob[i]
      sdv[i] <- sqrt(a*prob[i]*(1-prob[i]))
      g <- g + log(choose(a,y[i])) + y[i]*log(prob[i]) + (a-y[i])*log(1-prob[i])
    }
    g <- -1 * g
    return(g)
  }
  
  
  glarmagr.b <- function(delta){  # gradient
    
    gr <- rep(0, 3)
    out <- matrix(0, 3, 3)
    
    
    eta <- vector()
    m <- vector()
    prob <- vector()
    sdv <- vector()
    
    eta[1] <- 0
    prob[1] <- exp(eta[1])/(1 + exp(eta[1]))
    m[1] <- a*prob[1]
    sdv[1] <- sqrt(a*prob[1]*(1-prob[1]))
    
    beta <- delta[1]
    phi1 <- delta[2]
    theta1 <-  delta[3]
    
    for (i in 2:n){
      
      eta[i] <- beta + phi1*eta[i-1] + theta1*(y[i-1] - m[i-1])/sdv[i-1]
      prob[i] <- exp(eta[i])/(1 + exp(eta[i]))
      m[i] <- a*prob[i]
      sdv[i] <- sqrt(a*prob[i]*(1-prob[i]))
      
      gr[1] <- y[i]-m[i]
      gr[2] <- (y[i]-m[i])*eta[i-1]
      gr[3] <- (y[i]-m[i])*(y[i-1] - m[i-1])/sdv[i-1]
      
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
    
    model1b <- g_barma.b(a, delt)                
    #model1b <- g_garma.b(a, delt)
    #model1b <- g_glarma.b(a, delt)
    
    s_1b[,j] <- model1b$y
    p_1b[,j] <- model1b$pr
    e_1b[,j] <- model1b$x
    
    y <- s_1b[,j]
    
    
    # estimate models
    
    r <- 2
    
    d0 <- rep(0.01, 3)   # set initial values of the parameters
    
    # estimate BARMA
    
    
    d0[1] <- mean(y)
    
    ui <- rbind(c(0,0,-1),c(0,1,0),c(0,0,1))
    ci <- c(-1,0,0)

    num_mle <- constrOptim(theta=d0, f=barl.b, ui=ui, ci=ci, method="Nelder-Mead",
                           control=list(trace=1, maxit=10000, reltol=sqrt(.Machine$double.eps)))
    num_mle$hessian <- optimHess(num_mle$par, barl.b)
    
    num_par <- num_mle$par        # store result of estimation
    num_l <- -num_mle$value       # value of maximum (quasi) log-likelihood
    
    num_J <- solve(num_mle$hessian)  #compute Hessian                              
    
    num_I <- bargr.b(num_par)          #compute conditional information matrix
    
    num_var <- num_J%*%num_I%*%num_J/n # sandiwhich estimator
    num_V <- diag(num_var)             # estimated variance
    num_V <- t(num_V)
    num_SE <- sqrt(num_V)              # estimated SE
    num_t <- num_par/num_SE
    num_p <- pt(abs(num_t), n-length(num_mle$par), lower.tail=FALSE)  # t test
    
    num_aic <- 2*length(num_mle$par) - 2*num_l      # compute IC                      #AIC
    num_bic <- log(n)*length(num_mle$par) - 2*num_l                                   #BIC
    
    x <- vector()
    x[1] <- 0
    pro <- vector()
    mm <- vector()
    pro[1] <- exp(x[1])/(1+exp(x[1]))
    mm[1] <- a*pro[1]
    
    for (i in 2:n){
      
      x[i] <- num_par[1] + num_par[2]*y[i-1] + num_par[3]*(y[i-1]-mm[i-1])
      pro[i] <- exp(x[i])/(1+exp(x[i]))
      mm[i] <- a*pro[i]
      
    }
    
    
    # store results
    
    x_1bb[,j] <- x
    pro_1bb[,j] <- pro
    d_1bb[j,] <- num_par
    t_1bb[j,] <- num_t
    l_1b[j, 1] <- num_l
    aic_1b[j, 1] <- num_aic
    bic_1b[j, 1] <- num_bic
    
    
    # estimate GARMA 
    
    d0[1] <- mean(y)
    num_mle <- optim(d0, garmal.b, method="BFGS", hessian=T)
    
    num_par <- num_mle$par        # store result of estimation
    num_l <- -num_mle$value       # value of maximum (quasi) log-likelihood
    
    num_J <- solve(num_mle$hessian)  #compute Hessian                              
    
    num_I <- garmagr.b(num_par)          #compute conditional information matrix
    
    num_var <- num_J%*%num_I%*%num_J/n # sandiwhich estimator
    num_V <- diag(num_var)             # estimated variance
    num_V <- t(num_V)
    num_SE <- sqrt(num_V)              # estimated SE
    num_t <- num_par/num_SE
    num_p <- pt(abs(num_t), n-length(num_mle$par), lower.tail=FALSE)  # t test
    
    num_aic <- 2*length(num_mle$par) - 2*num_l      # compute IC                      #AIC
    num_bic <- log(n)*length(num_mle$par) - 2*num_l                                   #BIC
    
    x <- vector()
    x[1] <- 0
    y_ <- vector()
    gy <- vector()
    c <- 0.1
    
    y_[1] <- min(max(y[1], c), a-c)
    gy[1] <- log(y_[1]/(a-y_[1]))
    
    
    for (i in (r+1):n){
      
      y_[i] <- min(max(y[i], c), a-c)
      gy[i] <- log(y_[i]/(a-y_[i]))
      x[i] <- num_par[1] + num_par[2]*gy[i-1] + num_par[3]*(gy[i-1] - x[i-1])
      
    }
    
    pr <- exp(x)/(1+exp(x))
    
    # store results
    
    x_1bg[,j] <- x
    pro_1bg[,j] <- pr
    d_1bg[j,] <- num_par
    t_1bg[j,] <- num_t
    l_1b[j, 2] <- num_l
    aic_1b[j, 2] <- num_aic
    bic_1b[j, 2] <- num_bic
    
    # estimate GLARMA
    
    
    d0[1] <- mean(y)
    
    ui <- c(0,-1,-1)
    ci <- -1
    num_mle <- constrOptim(theta=d0, f=glarmal.b, ui=ui, ci=ci, method="Nelder-Mead",
                           control=list(trace=1, maxit=10000, reltol=sqrt(.Machine$double.eps)))
    num_mle$hessian <- optimHess(num_mle$par, glarmal.b)
    
    num_par <- num_mle$par        # store result of estimation
    num_l <- -num_mle$value       # value of maximum (quasi) log-likelihood
    
    num_aic <- 2*length(num_mle$par) - 2*num_l      # compute IC                      #AIC
    num_bic <- log(n)*length(num_mle$par) - 2*num_l                                   #BIC

    
    x <- vector()
    x[1] <- 0
    pro <- vector()
    mm <- vector()
    sv <- vector()
    
    pro[1] <- exp(x[1])/(1 + exp(x[1]))
    mm[1] <- a*pro[1]
    sv[1] <- sqrt(a*pro[1]*(1-pro[1]))
    
    for (i in 2:n){
      
      x[i] <- num_par[1] + num_par[2]*x[i-1] + num_par[3]*(y[i-1] - mm[i-1])/sv[i-1]
      pro[i] <- exp(x[i])/(1 + exp(x[i]))
      mm[i] <- a*pro[i]
      sv[i] <- sqrt(a*pro[i]*(1-pro[i]))
    }
    
    # store results
    
    x_1bgl[,j] <- x
    pro_1bgl[,j] <- pro
    d_1bgl[j,] <- num_par
    t_1bgl[j,] <- 1
    l_1b[j, 3] <- num_l
    aic_1b[j, 3] <- num_aic
    bic_1b[j, 3] <- num_bic
    
  }
  
  
  #summary statistics
  
  mean_1bb[1,] <- colMeans(d_1bb)          #mean of parametr estimate
  mean_1bg[1,] <- colMeans(d_1bg)
  mean_1bgl[1,] <- colMeans(d_1bgl)
  
  mean_1bb[2,] <- apply(d_1bb, 2, sd)      #stadard deviations
  mean_1bg[2,] <- apply(d_1bg, 2, sd)
  mean_1bgl[2,] <- apply(d_1bgl, 2, sd)
  
  freq<- function(x){
    sum(x>2.58)/s*100
  }
  
  mean_1bb[3,] <- apply(t_1bb, 2, FUN=freq)   # freq of rejecting tests
  mean_1bg[3,] <- apply(t_1bg, 2, FUN=freq)
  mean_1bgl[3,] <- apply(t_1bgl, 2, FUN=freq)
  
  mean <- cbind(mean_1bb, mean_1bg, mean_1bgl)
  
  aicc <- aic_1b==apply(aic_1b, 1, min)   # freq of correct selection
  ic[1,] <- 100*apply(aicc, 2, mean)
  bicc <- bic_1b==apply(bic_1b, 1, min)
  ic[2,] <- 100*apply(bicc, 2, mean)
  
  
  results <- list(eta=e_1b, prob=p_1b, simulation=s_1b,
                  d_1bb=d_1bb, d_1bg=d_1bg, d_1bgl=d_1bgl,
                  mean=mean, ic=ic, logl=l_1b, aic=aic_1b, bic=bic_1b)

  return(results)
  
}

# Results

set.seed(12346)

##################barma


simm_barma1 <- simm(a=5, n=250, s=1000, delt=c(0.2, 0.4, 0.2))
simm_barma2 <- simm(a=5, n=500, s=1000, delt=c(0.2, 0.4, 0.2))
simm_barma3 <- simm(a=5, n=1000, s=1000, delt=c(0.2, 0.4, 0.2))

round(simm_barma1$mean, 4)
round(simm_barma2$mean, 4)
round(simm_barma3$mean, 4)
round(simm_barma1$ic, 4)
round(simm_barma2$ic, 4)
round(simm_barma3$ic, 4)


####################garma


simm_garma1 <- simm(a=5, n=250, s=1000, delt=c(0.2, 0.4, 0.2))
simm_garma2 <- simm(a=5, n=500, s=1000, delt=c(0.2, 0.4, 0.2))
simm_garma3 <- simm(a=5, n=1000, s=1000, delt=c(0.2, 0.4, 0.2))

round(simm_garma1$mean, 4)
round(simm_garma2$mean, 4)
round(simm_garma3$mean, 4)
round(simm_garma1$ic, 4)
round(simm_garma2$ic, 4)
round(simm_garma3$ic, 4)


##################glarma


simm_glarma1 <- simm(a=5, n=250, s=1000, delt=c(0.2, 0.4, 0.3))
simm_glarma2 <- simm(a=5, n=500, s=1000, delt=c(0.2, 0.4, 0.3))
simm_glarma3 <- simm(a=5, n=1000, s=1000, delt=c(0.2, 0.4, 0.3))

round(simm_glarma1$mean, 4)
round(simm_glarma2$mean, 4)
round(simm_glarma3$mean, 4)
round(simm_glarma1$ic, 4)
round(simm_glarma2$ic, 4)
round(simm_glarma3$ic, 4)