##################################################################
# R code for Armillotta, Luati and Lupparelli 
# Observation driven models for discrete-valued time series
##################################################################

##################################################################
# Estimation, results and plots for applications
##################################################################

# load required libraries

library(numDeriv)
#library(pracma)
library(tscount)

#Choose one of the two datasets:

#################################################################

# Number of storms in the North Atlantic Basin

hurricane <- read.delim2("hurricane.txt", header=FALSE)
hurricane <- hurricane[1:168,]

y<-hurricane$V3
n<-length(y)

time<-seq(1:length(y))
tt <- seq(from=1851, to=2018, by=1)

# plot series

plot(y ~ tt, type = "l", xlab="Time", ylab="counts",
     main="Named storms counts")

z <- seq(1, max(y), by=10)

######################################################################

# Disease cases of Escherichia coli in North Rhine-Westphalia

y<-ecoli[,3]
n<-length(y)

time<-seq(1:length(y))
tta <- seq(as.Date("2001/01/01"), as.Date("2013/05/15"), by = "week")

# plot series

plot(y ~ tta, type = "l", xlab="Time", ylab="counts",     
     main="E.coli counts", xaxt = "n")
axis.Date(1, at = seq(as.Date("2001/1/1"), max(tta)+30, "year"))

z <- seq(1, max(y), by=20)

######################################################################

#select lag one models

p<-1
q<-1
r <- max(p,q)
k <- p+q+1 

# set initial values of the parameters

d0 <- rep(0, k)

# estimate Poisson models

log.pois11 <- ml(k, y, n, d0, logl.p, gradd.p, mod="log.p", me="bfgs")
garma.pois11 <- ml(k, y, n, d0, garmal.p, garmagr.p, mod="garma.p", me="bfgs")
glarma.pois11 <- ml(k, y, n, d0, glarmal.p, glarmagr.p, mod="glarma.p", me="bfgs")


# estimate Negative Binomial models
# and calibrate dispersion parameter

log.mean <- vector()
garma.mean <- vector()
glarma.mean <- vector()

log.v1 <- 0
garma.v1 <- 0
glarma.v1 <- 0

for(j in 1:n){
  
  log.mean[j] <- log.pois11$mean[j]
  garma.mean[j] <- garma.pois11$mean[j]
  glarma.mean[j] <- glarma.pois11$mean[j]
  
  log.v1 <- log.v1 + ((y[j]- log.mean[j])^2 - log.mean[j])/( log.mean[j]^2)
  garma.v1 <- garma.v1 + ((y[j]- garma.mean[j])^2 - garma.mean[j])/( garma.mean[j]^2)
  glarma.v1 <- glarma.v1 + ((y[j]- glarma.mean[j])^2- glarma.mean[j])/( glarma.mean[j]^2)
  
}

log.v1 <- log.v1/n
log.v1 <- 1/log.v1
garma.v1 <- garma.v1/n
garma.v1 <- 1/garma.v1
glarma.v1 <- glarma.v1/n
glarma.v1 <- 1/glarma.v1

log.v <- log.v1
garma.v <- garma.v1
glarma.v <- glarma.v1

old.v <- c(log.v1, garma.v1, glarma.v1)

tolerance <- 1e-3
tol <- matrix(tolerance, 3, 1)
con <- FALSE
while(con == FALSE){
  
  v <- log.v
  log.nb11 <- ml(k, y, n, d0, logl.nb, gradd.nb, mod="log.nb", me="bfgs")
  
  v <- garma.v
  garma.nb11 <- ml(k, y, n, d0, garmal.nb, garmagr.nb, mod="garma.nb", me="bfgs")
  
  v <- glarma.v
  glarma.nb11 <- ml(k, y, n, d0, glarmal.nb, glarmagr.nb, mod="glarma.nb", me="bfgs")

  log.v1 <- 0
  garma.v1 <- 0
  glarma.v1 <- 0
  
  for(j in 1:n){
    
    log.mean[j] <- log.nb11$mean[j]
    garma.mean[j] <- garma.nb11$mean[j]
    glarma.mean[j] <- glarma.nb11$mean[j]
    
    log.v1 <- log.v1 + ((y[j]- log.mean[j])^2 - log.mean[j])/( log.mean[j]^2)
    garma.v1 <- garma.v1 + ((y[j]- garma.mean[j])^2 - garma.mean[j])/( garma.mean[j]^2)
    glarma.v1 <- glarma.v1 + ((y[j]- glarma.mean[j])^2- glarma.mean[j])/( glarma.mean[j]^2)
    
  }
  
  log.v1 <- log.v1/n
  log.v1 <- 1/log.v1
  garma.v1 <- garma.v1/n
  garma.v1 <- 1/garma.v1
  glarma.v1 <- glarma.v1/n
  glarma.v1 <- 1/glarma.v1
  
  log.v <- log.v1
  garma.v <- garma.v1
  glarma.v <- glarma.v1
  
  v.vec <- c(log.v, garma.v, glarma.v)
  
  ifelse(abs(v.vec - old.v) < tol, con <- TRUE, con <- FALSE)
  
  old.v <- v.vec
}

v <- log.v
log.nb11 <- ml(k, y, n, d0, logl.nb, gradd.nb, mod="log.nb", me="bfgs")

v <- garma.v
garma.nb11 <- ml(k, y, n, d0, garmal.nb, garmagr.nb, mod="garma.nb", me="bfgs")

v <- glarma.v
glarma.nb11 <- ml(k, y, n, d0, glarmal.nb, glarmagr.nb, mod="glarma.nb", me="bfgs")

# Probabilistic calibration

log.pois.Px <- vector()
log.pois.Px1 <- vector()
log.nb.Px <- vector()
log.nb.Px1 <- vector()
log.pois.px <- vector()
log.nb.px <- vector()

garma.pois.Px <- vector()
garma.pois.Px1 <- vector()
garma.nb.Px <- vector()
garma.nb.Px1 <- vector()
garma.pois.px <- vector()
garma.nb.px <- vector()

glarma.pois.Px <- vector()
glarma.pois.Px1 <- vector()
glarma.nb.Px <- vector()
glarma.nb.Px1 <- vector()
glarma.pois.px <- vector()
glarma.nb.px <- vector()

for(i in 1:n){
  
  log.pois.Px[i] <- ppois(y[i], log.pois11$mean[i])
  log.pois.Px1[i] <- ppois(y[i]-1, log.pois11$mean[i])
  log.nb.Px[i] <- pnbinom(y[i], size=log.v, prob=log.v/(log.v+log.nb11$mean[i]))
  log.nb.Px1[i] <- pnbinom(y[i]-1, size=log.v, prob=log.v/(log.v+log.nb11$mean[i]))
  log.pois.px[i] <- dpois(y[i], log.pois11$mean[i])
  log.nb.px[i] <- dnbinom(y[i], size=log.v, prob=log.v/(log.v+log.nb11$mean[i]))
  
  garma.pois.Px[i] <- ppois(y[i], garma.pois11$mean[i])
  garma.pois.Px1[i] <- ppois(y[i]-1, garma.pois11$mean[i])
  garma.nb.Px[i] <- pnbinom(y[i], size=garma.v, prob=garma.v/(garma.v+garma.nb11$mean[i]))
  garma.nb.Px1[i] <- pnbinom(y[i]-1, size=garma.v, prob=garma.v/(garma.v+garma.nb11$mean[i]))
  garma.pois.px[i] <- dpois(y[i], garma.pois11$mean[i])
  garma.nb.px[i] <- dnbinom(y[i], size=garma.v, prob=garma.v/(garma.v+garma.nb11$mean[i]))
  
  glarma.pois.Px[i] <- ppois(y[i], glarma.pois11$mean[i])
  glarma.pois.Px1[i] <- ppois(y[i]-1, glarma.pois11$mean[i])
  glarma.nb.Px[i] <- pnbinom(y[i], size=glarma.v, prob=glarma.v/(glarma.v+glarma.nb11$mean[i]))
  glarma.nb.Px1[i] <- pnbinom(y[i]-1, size=glarma.v, prob=glarma.v/(glarma.v+glarma.nb11$mean[i]))
  glarma.pois.px[i] <- dpois(y[i], glarma.pois11$mean[i])
  glarma.nb.px[i] <- dnbinom(y[i], size=glarma.v, prob=glarma.v/(glarma.v+glarma.nb11$mean[i]))
  
}


# Marginal calibration


log.Px.p <- rep(0, length(z))
log.Px.nb <- rep(0, length(z))

Gx <- rep(0, length(z))

garma.Px.p <- rep(0, length(z))
garma.Px.nb <- rep(0, length(z))

glarma.Px.p <- rep(0, length(z))
glarma.Px.nb <- rep(0, length(z))


for (j in 1:length(z)){
  
  for(i in 1:n){
    
    log.Px.p[j] <- log.Px.p[j] + ppois(z[j], log.pois11$mean[i])
    log.Px.nb[j] <- log.Px.nb[j] + pnbinom(z[j], size=log.v, prob=log.v/(log.v+log.nb11$mean[i]))
    
    garma.Px.p[j] <- garma.Px.p[j] + ppois(z[j], garma.pois11$mean[i])
    garma.Px.nb[j] <- garma.Px.nb[j] + pnbinom(z[j], size=garma.v, prob=garma.v/(garma.v+garma.nb11$mean[i]))
    
    glarma.Px.p[j] <- glarma.Px.p[j] + ppois(z[j], glarma.pois11$mean[i])
    glarma.Px.nb[j] <- glarma.Px.nb[j] + pnbinom(z[j], size=glarma.v, prob=glarma.v/(glarma.v+glarma.nb11$mean[i]))
    
    Gx[j] <- Gx[j] + ifelse(y[i]<=z[j],1,0)
    
  }
  
  log.Px.p[j] <- log.Px.p[j]/n
  log.Px.nb[j] <- log.Px.nb[j]/n
  
  garma.Px.p[j] <- garma.Px.p[j]/n
  garma.Px.nb[j] <- garma.Px.nb[j]/n
  
  glarma.Px.p[j] <- glarma.Px.p[j]/n
  glarma.Px.nb[j] <- glarma.Px.nb[j]/n
  
  Gx[j] <- Gx[j]/n
  
}

log.mc.p <- log.Px.p - Gx
log.mc.nb <- log.Px.nb - Gx

garma.mc.p <- garma.Px.p - Gx
garma.mc.nb <- garma.Px.nb - Gx

glarma.mc.p <- glarma.Px.p - Gx
glarma.mc.nb <- glarma.Px.nb - Gx


# Sharpness via scoring rules

kk <- 1000                           

log.pois.logs <- - log(log.pois.px)
garma.pois.logs <- - log(garma.pois.px)
glarma.pois.logs <- - log(glarma.pois.px)

log.nb.logs <- - log(log.nb.px)
garma.nb.logs <- - log(garma.nb.px)
glarma.nb.logs <- - log(glarma.nb.px)

log.pois.norm <- rep(0, n)
garma.pois.norm <- rep(0, n)
glarma.pois.norm <- rep(0, n)

log.nb.norm <- rep(0, n)
garma.nb.norm <- rep(0, n)
glarma.nb.norm <- rep(0, n)

for(i in 1:n){
  
  for(j in 1:kk){
    
    log.pois.norm[i] <- log.pois.norm[i] + dpois(j, log.pois11$mean[i])^2
    log.nb.norm[i] <- log.nb.norm[i] + dnbinom(j, size=log.v, prob=log.v/(log.v+log.nb11$mean[i]))^2
    
    garma.pois.norm[i] <- garma.pois.norm[i] + dpois(j, garma.pois11$mean[i])^2
    garma.nb.norm[i] <- garma.nb.norm[i] + dnbinom(j, size=garma.v, prob=garma.v/(garma.v+garma.nb11$mean[i]))^2
    
    glarma.pois.norm[i] <- glarma.pois.norm[i] + dpois(j, glarma.pois11$mean[i])^2
    glarma.nb.norm[i] <- glarma.nb.norm[i] + dnbinom(j, size=glarma.v, prob=glarma.v/(glarma.v+glarma.nb11$mean[i]))^2
    
    
  }
  
}

log.pois.qs <- - 2*log.pois.px + log.pois.norm
garma.pois.qs <- - 2*garma.pois.px + garma.pois.norm
glarma.pois.qs <- - 2*glarma.pois.px + glarma.pois.norm

log.nb.qs <- - 2*log.nb.px + log.nb.norm
garma.nb.qs <- - 2*garma.nb.px + garma.nb.norm
glarma.nb.qs <- - 2*glarma.nb.px + glarma.nb.norm

log.pois.sphs <- - log.pois.px / sqrt(log.pois.norm)
garma.pois.sphs <- - garma.pois.px / sqrt(garma.pois.norm)
glarma.pois.sphs <- - glarma.pois.px / sqrt(glarma.pois.norm)

log.nb.sphs <- - log.nb.px / sqrt(log.nb.norm)
garma.nb.sphs <- - garma.nb.px / sqrt(garma.nb.norm)
glarma.nb.sphs <- - glarma.nb.px / sqrt(glarma.nb.norm)

log.diff.p <- rep(0, n)
log.diff.nb <- rep(0, n)
diffx <- rep(0, n)
log.rps.p <- rep(0, n)
log.rps.nb <- rep(0, n)

garma.diff.p <- rep(0, n)
garma.diff.nb <- rep(0, n)
garma.rps.p <- rep(0, n)
garma.rps.nb <- rep(0, n)

glarma.diff.p <- rep(0, n)
glarma.diff.nb <- rep(0, n)
glarma.rps.p <- rep(0, n)
glarma.rps.nb <- rep(0, n)

for(i in 1:n){
  
  for(j in 1:kk){
    
    log.diff.p[i] <- ppois(j, log.pois11$mean[i])
    log.diff.nb[i] <- pnbinom(j, size=log.v, prob=log.v/(log.v+log.nb11$mean[i]))
    diffx[i] <- ifelse(y[i]<=j,1,0)
    log.rps.p[i]  <- log.rps.p[i] + (log.diff.p[i] - diffx[i])^2
    log.rps.nb[i]  <- log.rps.nb[i] + (log.diff.nb[i] - diffx[i])^2
    
    garma.diff.p[i] <- ppois(j, garma.pois11$mean[i])
    garma.diff.nb[i] <- pnbinom(j, size=garma.v, prob=garma.v/(garma.v+garma.nb11$mean[i]))
    garma.rps.p[i]  <- garma.rps.p[i] + (garma.diff.p[i] - diffx[i])^2
    garma.rps.nb[i]  <- garma.rps.nb[i] + (garma.diff.nb[i] - diffx[i])^2
    
    glarma.diff.p[i] <- ppois(j, glarma.pois11$mean[i])
    glarma.diff.nb[i] <- pnbinom(j, size=glarma.v, prob=glarma.v/(glarma.v+glarma.nb11$mean[i]))
    glarma.rps.p[i]  <- glarma.rps.p[i] + (glarma.diff.p[i] - diffx[i])^2
    glarma.rps.nb[i]  <- glarma.rps.nb[i] + (glarma.diff.nb[i] - diffx[i])^2
    
  }
  
}

log.pois.rps <- log.rps.p
garma.pois.rps <- garma.rps.p
glarma.pois.rps <- glarma.rps.p

log.nb.rps <- log.rps.nb
garma.nb.rps <- garma.rps.nb
glarma.nb.rps <- glarma.rps.nb


######################### Results & plots #####################

# QMLE results

round(c(log.pois11$est, 0 , log.pois11$aic, log.pois11$bic), 3)
round(log.pois11$std,3)
round(c(garma.pois11$est, 0, garma.pois11$aic, garma.pois11$bic), 3)
round(garma.pois11$std,3)
round(c(glarma.pois11$est, 0, glarma.pois11$aic, glarma.pois11$bic), 3)
round(glarma.pois11$std,3)
round(c(log.nb11$est, log.v1, log.nb11$aic, log.nb11$bic), 3)
round(log.nb11$std,3)
round(c(garma.nb11$est, garma.v1, garma.nb11$aic, garma.nb11$bic), 3)
round(garma.nb11$std,3)
round(c(glarma.nb11$est, glarma.v1, glarma.nb11$aic, glarma.nb11$bic), 3)
round(glarma.nb11$std,3)

#compute QIC 

log.pois11$hessian <- hessian(logl.p, x=log.pois11$est)
log.pois11$var <- solve(log.pois11$hessian)%*%gradd.p(log.pois11$est)%*%solve(log.pois11$hessian)/n
log.pois11$qic <- 2*sum(diag(log.pois11$hessian%*%log.pois11$var*n)) - 2*log.pois11$max   #QIC

garma.pois11$hessian <- hessian(garmal.p, x=garma.pois11$est)
garma.pois11$var <- solve(garma.pois11$hessian)%*%garmagr.p(garma.pois11$est)%*%solve(garma.pois11$hessian)/n
garma.pois11$qic <- 2*sum(diag(garma.pois11$hessian%*%garma.pois11$var*n)) - 2*garma.pois11$max    #QIC

glarma.pois11$hessian <- hessian(glarmal.p, x=glarma.pois11$est)
glarma.pois11$var <- solve(glarma.pois11$hessian)%*%glarmagr.p(glarma.pois11$est)%*%solve(glarma.pois11$hessian)/n
glarma.pois11$qic <- 2*sum(diag(glarma.pois11$hessian%*%glarma.pois11$var*n)) - 2*glarma.pois11$max    #QIC

log.nb11$hessian <- hessian(logl.nb, x=log.nb11$est)
log.nb11$var <- solve(log.nb11$hessian)%*%gradd.nb(log.nb11$est)%*%solve(log.nb11$hessian)/n
log.nb11$qic <- 2*sum(diag(log.nb11$hessian%*%log.nb11$var*n)) - 2*log.nb11$max   #QIC

garma.nb11$hessian <- hessian(garmal.nb, x=garma.nb11$est)
garma.nb11$var <- solve(garma.nb11$hessian)%*%garmagr.nb(garma.nb11$est)%*%solve(garma.nb11$hessian)/n
garma.nb11$qic <- 2*sum(diag(garma.nb11$hessian%*%garma.nb11$var*n)) - 2*garma.nb11$max    #QIC

glarma.nb11$hessian <- hessian(glarmal.nb, x=glarma.nb11$est)
glarma.nb11$var <- solve(glarma.nb11$hessian)%*%glarmagr.nb(glarma.nb11$est)%*%solve(glarma.nb11$hessian)/n
glarma.nb11$qic <- 2*sum(diag(glarma.nb11$hessian%*%glarma.nb11$var*n)) - 2*glarma.nb11$max    #QIC

round(c(log.pois11$qic, garma.pois11$qic , glarma.pois11$qic), 3)
round(c(log.nb11$qic, garma.nb11$qic , glarma.nb11$qic), 3)

#result of scoring rules

# logarithmic  # quadratic  # spherical  # ranked probability 

round(c(mean(log.pois.logs), mean(log.pois.qs), mean(log.pois.sphs), mean(log.pois.rps)), 4)
round(c(mean(log.nb.logs), mean(log.nb.qs), mean(log.nb.sphs), mean(log.nb.rps)), 4)

round(c(mean(garma.pois.logs), mean(garma.pois.qs), mean(garma.pois.sphs), mean(garma.pois.rps)), 4)
round(c(mean(garma.nb.logs), mean(garma.nb.qs), mean(garma.nb.sphs), mean(garma.nb.rps)), 4)

round(c(mean(glarma.pois.logs), mean(glarma.pois.qs), mean(glarma.pois.sphs), mean(glarma.pois.rps)), 4)
round(c(mean(glarma.nb.logs), mean(glarma.nb.qs), mean(glarma.nb.sphs), mean(glarma.nb.rps)), 4)


# plot of autocorrelation function (acf) and marginal calibration (mc)

par(mfrow=c(1,1))
acf<- acf(y, xlab="Lag", ylab="ACF", main="ACF")
plot(log.mc.p ~ z, type = "l", lty=2, xlab="x", ylab="", main="log-AR mc plot", ylim=c(-0.04,0.03))
lines(log.mc.nb ~ z, type = "l", lty=1)
abline(h=0, lty = 3)
plot(glarma.mc.p ~ z, type = "l", lty=2, xlab="x", ylab="", main="GLARMA mc plot", ylim=c(-0.04,0.03))
lines(glarma.mc.nb ~ z, type = "l", lty=1)
abline(h=0, lty = 3)


# plot PIT

par(mfrow=c(2,3))
pit(y, log.pois.Px, log.pois.Px1, n.bins=10, y.max=3, my.title="PIT Poisson log-AR")
pit(y, garma.pois.Px, garma.pois.Px1, n.bins=10, y.max=3, my.title="PIT Poisson GARMA")
pit(y, glarma.pois.Px, glarma.pois.Px1, n.bins=10, y.max=3, my.title="PIT Poisson GLARMA")
pit(y, log.nb.Px, log.nb.Px1, n.bins=10, y.max=3, my.title="PIT NB log-AR")
pit(y, garma.nb.Px, garma.nb.Px1, n.bins=10, y.max=3, my.title="PIT NB GARMA")
pit(y, glarma.nb.Px, glarma.nb.Px1, n.bins=10, y.max=3, my.title="PIT NB GLARMA")

#-----------------------------------------------------------


