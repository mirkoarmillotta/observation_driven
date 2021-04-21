##################################################################
# R code for Armillotta, Luati and Lupparelli 
# Observation driven models for discrete-valued time series
##################################################################

###################################################################
# Functions
###################################################################

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
  
  gr <- rep(0, k)
  out <- matrix(0, k, k)
  
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
  
  gr <- rep(0, k)
  out <- matrix(0, k, k)
  
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
  
  gr <- rep(0, k)
  out <- matrix(0, k, k)
  
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
  
  gr <- rep(0, k)
  out <- matrix(0, k, k)
  
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
  
  gr <- rep(0, k)
  out <- matrix(0, k, k)
  
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
  
  gr <- rep(0, k)
  out <- matrix(0, k, k)
  
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


#-----------------------------------------------------------------

# function for QMLE estimation

ml <- function(k, y, n, d0, logl, gradd, mod, met){
  
  r <- max(p,q)
  
  # select different optimization algorithm
  
  if(met=="nm")
    
    num_mle <- optim(par=d0, fn=logl, method="Nelder-Mead",
                     hessian=TRUE,
                     control=list(trace=1, maxit=10000))
  
  if(met=="cg")
    
    num_mle <- optim(par=d0, fn=logl, method="CG",
                     hessian=TRUE,
                     control=list(trace=1, maxit=10000))
  
  if(met=="l-bfgs-b")
    
    num_mle <- optim(par=d0, fn=logl, method="L-BFGS-B",
                     lower=c(0, 0, 0, 0),
                     upper=c(Inf, Inf, Inf, Inf),
                     hessian=TRUE,
                     control=list(trace=1, maxit=100000))
  
  if(met=="bfgs")
    
    num_mle <- optim(par=d0, fn=logl, method="BFGS",
                     hessian=TRUE,
                     control=list(trace=1, maxit=10000))
  
  if(met=="constr"){
    
    ui <- rbind(c(0,-1,-1,0),c(0,-1,0,-1),c(0,1,-1,0))
    ci <- c(-1,-1,-1)
    num_mle <- constrOptim(theta=d0, f=logl,
                           ui=ui, ci=ci,
                           method="Nelder-Mead",
                           control=list(trace=1, maxit=10000))
    
    num_mle$hessian <- optimHess(num_mle$par, logl)
    
  }
  
  if(met=="constr2"){
    
    ui <- c(0, -1, -1)
    ci <- -0.999
    num_mle <- constrOptim(theta=d0, f=logl, 
                           ui=ui, ci=ci, method="Nelder-Mead",
                           control=list(trace=1, maxit=10000))
    
    num_mle$hessian <- hessian(logl, x=num_mle$par)
    
  }
  
  # store result of estimation
  # value of maximum (quasi) log-likelihood
  # compute IC
  
  num_par <- num_mle$par        
  num_l <- -num_mle$value       
  num_aic <- 2*k - 2*num_l      
  num_bic <- log(n)*k - 2*num_l
  
  #compute inverse Hessian (matrix J^(-1)) 
  #compute outer product (Matrix I)
  # sandiwhich estimator i.e. matrix J^(-1)IJ^(-1)
  # estimated variance
  # estimated SE
  # t test
  
  num_J <- solve(num_mle$hessian)                               
  
  num_I <- gradd(num_par)          
  
  num_var <- num_J%*%num_I%*%num_J/n 
  num_V <- diag(num_var)             
  num_V <- t(num_V)
  num_SE <- sqrt(num_V)              
  num_t <- num_par/num_SE
  num_p <- pt(abs(num_t), n-k, lower.tail=FALSE) 
  
  # compute fitted values
  
  x <- vector()
  x[1] <- 0
  
  if(mod=="log.p" || mod=="log.nb")
    
    for (i in 2:n){
      
      x[i] <- num_par[1] + num_par[2]*log(y[i-1]+1) +
        num_par[3]*x[i-1]
      
    }
  
  if(mod=="garma.p" || mod=="garma.nb"){
    
    y_ <- vector()
    gy <- vector()
    c <- 0.1
    
    y_[1] <- max(y[1], c)
    gy[1] <- log(y_[1])
    
    for (i in 2:n){
      
      y_[i] <- max(y[i], c)
      gy[i] <- log(y_[i])
      x[i] <- num_par[1] + num_par[2]*gy[i-1] +
        num_par[3]*(gy[i-1] - x[i-1])
      
    }
    
  }
  
  if(mod=="glarma.p"){
    
    lambda <- vector()
    lambda[1] <- 1
    
    for (i in 2:n){
      
      x[i] <- num_par[1] + num_par[2]*x[i-1] +
        num_par[3]*
        (y[i-1] - lambda[i-1])/sqrt(lambda[i-1])
      
      lambda[i] <- exp(x[i])
      
    }
    
  }
  
  if(mod=="glarma.nb"){
    
    lambda <- vector()
    var <- vector()
    
    lambda[1] <- 1
    var[1] <- 1+1/v
    
    for (i in 2:n){
      
      x[i] <- num_par[1] + num_par[2]*x[i-1] +
        num_par[3]*
        (y[i-1] - lambda[i-1])/sqrt(var[i-1])
      
      lambda[i] <- exp(x[i])
      var[i] <- lambda[i]*(1+lambda[i]/v)
      
    }
    
  }
  
  
  # compute residuals and print results
  
  
  mu <- exp(x)
  
  if(mod=="log.p" || mod=="garma.p" ||
     mod=="glarma.p" ||
     mod=="gen.p" ||
     mod=="gen1.p") 
    
    e <- (y-mu)/sqrt(mu)
  
  else e <- (y-mu)/sqrt(mu*(1+mu/v))
  
  
  return(list(est=num_par, std=num_SE, t.test=num_t,
              p.value=num_p, mean=mu, pred=x,
              res=e, max=num_l,
              aic=num_aic, bic=num_bic))
  
}

#--------------------------------------------------------------

### function for nonrandomized PIT histogram 
###
### input: 
###   x    observed data 
###   Px   CDF at x 
###   Px1  CDF at x-1 

### This function is taken from the R code
### by Czado, Gneiting and Held, Biometrics (2009)

pit <- function(x, Px, Px1, n.bins=10, y.max=2.75, my.title="PIT Histogram")
{
  a.mat <- matrix(0,n.bins,length(x))
  k.vec <- pmax(ceiling(n.bins*Px1),1)
  m.vec <- ceiling(n.bins*Px)
  d.vec <- Px-Px1
  for (i in 1:length(x))
  {
    if (k.vec[i]==m.vec[i]) {a.mat[k.vec[i],i]=1}
    else 
    { 
      a.mat[k.vec[i],i]=((k.vec[i]/n.bins)-Px1[i])/d.vec[i]
      if ((k.vec[i]+1)<=(m.vec[i]-1))
      {for (j in ((k.vec[i]+1):(m.vec[i]-1))) {a.mat[j,i]=(1/(n.bins*d.vec[i]))}}
      a.mat[m.vec[i],i]=(Px[i]-((m.vec[i]-1)/n.bins))/d.vec[i]     
    }
  }
  a <- apply(a.mat,1,sum)
  a <- (n.bins*a)/(length(x))
  p <- (0:n.bins)/n.bins
  PIT <- "PIT"
  RF <- "Freq."
  plot(p, p, ylim=c(0,y.max), type="n", xlab=PIT, ylab=RF, main=my.title) 
  temp1 <- ((1:n.bins)-1)/n.bins
  temp2 <- ((1:n.bins)/n.bins)
  o.vec <- rep(0,n.bins)
  segments(temp1,o.vec,temp1,a)
  segments(temp1,a,temp2,a)
  segments(temp2,o.vec,temp2,a)
  segments(0,0,1,0)
}

#-------------------------------------------------------------------