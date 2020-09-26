rm(list = ls())

source("DataGen4.R")

library(lmtest)
library(ivpack)
probitboots <- function(trial){
  bootse <- c()
  rows = sample(1:1000, 1000, replace = TRUE)
  bootdat <- dat[rows, ]
  booty <- as.matrix(bootdat[, 1], ncol=1, nrow=n)
  bootx <- as.matrix(bootdat[, 2], ncol=1, nrow=n)
  bootxo <- as.matrix(bootdat[, (3+instrument):(2+exo+instrument)], ncol=exo, nrow=n)
  bootz <- as.matrix(bootdat[, 3:(2+instrument)], ncol=instrument, nrow=n)
  bootcf1 <- lm(bootx~bootz+bootxo)
  bootv = bootcf1$residuals
  bootcf2 <- glm(booty ~ bootx + bootxo + bootv , family = binomial(link = "probit"), control=list(epsilon = 1e-8, maxit = 100))
  bootse <- bootcf2$coefficients[2]
  return(bootse)
}

probit <- function(datmat, exo=1, instrument=1){
  r6 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
  results <- matrix(0, nrow=ncol(datmat), ncol= 1)
  
  c6 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
  coverage <- matrix(0, nrow=ncol(datmat), ncol= 1)
  
  cover <- function(estimate, se){
    upper <- estimate + 1.96*se
    lower <- estimate - 1.96*se
    if (.5 > lower & .5 < upper){
      return(1)}
    else{
      return(0)}
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
      
      dat2 =as.data.frame(cbind(y_values, x,z,xo))
      
      probitcf1 <- lm(x~z+xo)
      ## Collect residuals
      v = probitcf1$residuals
      cfdata = as.data.frame(cbind(y_values, x, xo, v))
      probitcf2 <- glm(y_values ~. , family = binomial(link = "probit"), data = cfdata, control=list(epsilon = 1e-8, maxit = 100))
      r6[i, j] <- probitcf2$coefficients[2]
      
      trials <- seq(1:500)
      results <- as.data.frame(lapply(trials,'probitboots'))
      cprobit <- sd(results[1,])
      c6[i,j] <- cover(estimate = r6[i,j], se=cprobit)      
    }  
    results[j, 1] <- mean(abs(r6[, j]-0.5))
    
    coverage[j, 1] <- sum(c6[,j])
    
  }
  return(list(results =results, coverage=coverage ))  
}

mad1 <- probit(datmat=data1)
mad1$results
mad1$coverage