rm(list = ls())

source("DataGen1.R")
data1 = simulateiv(size=1050, rhoxz = 0.1, rhoxe = c(.1,.2,0.3,.4,.5))

maxlik <- function(datmat, exo=1, instrument=1){
  
  n =1000
  r8 <- matrix(1, nrow=nrow(datmat), ncol= ncol(datmat))
  results <- matrix(0, nrow=ncol(datmat), ncol= 1)
  
  c8 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
  coverage <- matrix(0, nrow=ncol(datmat), ncol= 1)
  
  code <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
  
  e8 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat)) 
  endo <- matrix(0, nrow=ncol(datmat), ncol = 1)
  
  cover <- function(estimate, se){
    upper <- estimate + 1.96*se
    lower <- estimate - 1.96*se
    if (.5 > lower & .5 < upper){
      return(1)}
    else{
      return(0)}
  }
  
  lik=function(theta){
    
    alpha1<-theta[1]
    beta<-theta[2]
    alpha2<-theta[3]
    logsig<-theta[4]
    iht<-theta[5]
    gamma<-as.matrix(theta[6:(5+exo)], ncol=1, nrow= exo)
    pi1<- as.matrix(theta[(6+exo):(5+2*exo)], ncol = 1, nrow=exo)
    pi2<- as.matrix(theta[(6+2*exo):(5+2*exo+instrument)], ncol =1, nrow=instrument)
    
    y = y_values
    
    rho = (exp(2*iht)-1)/(exp(2*iht)+1)
    #rho = tanh(iht)
    m = (alpha1+x*beta+xo%*%gamma+(rho*(x-alpha2-xo%*%pi1-z%*%pi2)/exp(logsig))/(sqrt(1-(rho^2))))
    
    logl = y*log(pnorm(m,mean=0,sd=1))+(1-y)*log(1-pnorm(m,mean=0,sd=1))+log(dnorm(((x-alpha2-xo%*%pi1-z%*%pi2)/exp(logsig)),mean=0,sd=1))-logsig
    
    return(-sum(logl))
  }
  
  maxboots <- function(bootsize=100, eevs= 1, exo=1, instrument=1){
    bootse <- c()
    lol <- c()
    bootrho <- c()
    for (p in 1:bootsize){
      rows = sample(1:n, 1000, replace = TRUE)
      bootdat <- as.data.frame(dat[rows, ])
      booty <- as.matrix(bootdat[, 1], ncol=1, nrow=n)
      bootx <- as.matrix(bootdat[, 2], ncol=1, nrow=n)
      bootxo <- as.matrix(bootdat[, (3+instrument):(2+exo+instrument)], ncol=exo, nrow=n)
      bootz <- as.matrix(bootdat[, 3:(2+instrument)], ncol=instrument, nrow=n)
      
      bootlik=function(theta){
        alpha1<-theta[1]
        beta<-theta[2]
        alpha2<-theta[3]
        logsig<-theta[4]
        iht<-theta[5]
        gamma<-as.matrix(theta[6:(5+exo)], ncol=1, nrow= exo)
        pi1<- as.matrix(theta[(6+exo):(5+2*exo)], ncol = 1, nrow=exo)
        pi2<- as.matrix(theta[(6+2*exo):(5+2*exo+instrument)], ncol =1, nrow=instrument)
        y = booty[1]
        rho = (exp(2*iht)-1)/(exp(2*iht)+1)
        #rho = tanh(iht)
        m = (alpha1+bootx*beta+bootxo%*%gamma+(rho*(bootx-alpha2-bootxo%*%pi1-bootz%*%pi2)/exp(logsig))/(sqrt(1-(rho^2))))
        logl = y*log(pnorm(m,mean=0,sd=1))+(1-y)*log(1-pnorm(m,mean=0,sd=1))+log(dnorm(((bootx-alpha2-bootxo%*%pi1-bootz%*%pi2)/exp(logsig)),mean=0,sd=1))-logsig
        return(-sum(logl))
      }
      
      tryCatch({
        out<-nlm(bootlik,start,iterlim=1000)
        eststore<-out$estimate
        estgrad<-out$gradient
        convcode<-out$code}, error=function(e){cat("Error:",conditionMessage(e), "\n")})
      if (out$code == 1){
        lol[p] <- out$code
        bootse[p] <- out$estimate[2] 
        bootrho[p] <- out$estimate[5]
      } else {
        bootse[p] <- NA
        lol[p] <- NA
        bootrho[p] <- NA
        bootsize = bootsize+1  
      }
    }
    bootse <- bootse[which(complete.cases(bootse))]
    lol <-  lol[which(complete.cases(lol))]
    bootrho <- bootrho[which(complete.cases(bootrho))]
    return(list(est = sd(bootse), code=max(lol), rhovar= var(bootrho)))
  }
  
  for (j in 1:ncol(datmat)){  
    for (i in 1:nrow(datmat)){
      dat= unlist(datmat[i,j])
      dat = matrix(dat, ncol=5 ,nrow = 1000)
      n = 1000
      y_values = dat[,1]
      ypre = dat[,2]
      x <- dat[, 3]
      ##Obtain IV z (excluded exogenous regressor)
      z <- dat[, (4):(3+instrument)]
      ##Obtain included exogenous regressor
      xo <- dat[, (4+instrument):(3+ instrument+exo)]
      
      dat =as.data.frame(cbind(y_values, x,z,xo))
      
      # Vector of starting values
      start<-c(0,0,0,1,.5)
      start[6:(5+exo)] <- 0
      start[(6+exo):(5+2*exo)] <- 0
      start[(6+2*exo):(5+2*exo+instrument)] <- 0
      
      # tryCatch silences error messages to ensure that loop doesn't end prematurely 
      tryCatch({
        out<-nlm(lik,start,iterlim=1000)
        eststore<-out$estimate
        estgrad<-out$gradient
        convcode<-out$code}, error=function(e){cat("Error:",conditionMessage(e), "\n")})
      
      cfse = maxboots()
      
      if (convcode == 1){
        r8[i, j] <- out$estimate[2]
        
        c8[i,j] <- cover(estimate = r8[i,j], se=cfse$est)     
        
        rho <- out$estimate[5]
        critval = (rho)^2/cfse$rhovar
        e8[i,j] <- ifelse(critval > 3.841,1, 0)
        
        code[i,j] <- convcode} else{
        c8[i,j] <- NA
        e8[i,j] <- NA
        r8[i, j] <- NA
      }
    }  
    r8 <- r8[which(complete.cases(r8)),]
    c8 <- c8[which(complete.cases(c8)),]
    e8 <- e8[which(complete.cases(e8)),]
    r8 <- r8[1:1000,]
    c8 <- c8[1:1000,]
    e8 <- e8[1:1000,]
    results[j, 1] <- mean(abs(r8[, j]-0.5))
    
    coverage[j, 1] <- sum(c8[,j])
    endo[j,1] = sum(e8[,j])
    
  }
  return(list(results =results, coverage=coverage, code =max(code), endo=endo, outersize = nrow(e8) ))  
}

sink("NULL")
mad1 <- maxlik(datmat=data1)
sink()
mad1$results
mad1$coverage
mad1$code
mad1$endo
mad1$outersize