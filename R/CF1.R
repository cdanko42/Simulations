rm(list = ls())

source("DataGen1.R")
library(ks)
library(lmtest)
library(ivpack)

specreg <- function(datmat, exo=1, instrument=1){
  cfboots <- function(bootsize=599, eevs= 1, exo=1, instrument=1){
    bootse <- c()
    bootse2 <- c()
    for (p in 1:bootsize){
      rows = sample(1:n, 1000, replace = TRUE)
      bootdat <- as.data.frame(dat[rows, ])
      booty <- as.matrix(bootdat[, 1], ncol=1, nrow=n)
      bootx <- as.matrix(bootdat[, 2], ncol=1, nrow=n)
      bootxo <- as.matrix(bootdat[, (3+instrument):(2+exo+instrument)], ncol=exo, nrow=n)
      bootz <- as.matrix(bootdat[, 3:(2+instrument)], ncol=instrument, nrow=n)
      bootdat$demeanbootxo <- (bootxo[,1]-mean(bootxo[,1]))
      demeanbootxo <- bootdat$demeanbootxo
      
      ##Obtain residuals, which are our U_i
      specialregdataboot = as.data.frame(cbind(demeanbootxo, bootx, bootz))
      if (exo > 1){
        specialregdataboot = cbind(specialregdata, bootxo[, 2:exo])
      }
      
      fssr <- lm(demeanbootxo ~., data=specialregdata)
      bootdat$u <- fssr$residuals
      
      ##Step 2
      fhat <- kde(bootdat$u,eval.points=bootdat$u)$estimate
      
      ##Step 3
      ## Create vector that will collect our I(Vi > 0)
      bootdat$idv <- rep(0,nrow(bootdat))
      ## Collect identity function of v
      bootdat$idv <- replace(bootdat$idv,which(bootdat$demeanbootxo>=0),1)
      
      ##Obtain Ti; y_values=Di, fi and idv is defined above. 
      bootdat$t <- (bootdat[,1]- bootdat$idv)/fhat
      
      ## Trim extreme values
      lower = .025*n + 1
      upper = .975*n
      
      sortdatboot <- bootdat[order(bootdat$t),]
      trimdatboot <- sortdatboot[lower:upper,]
      firststagespecboot <- trimdatboot[, (3):(2 + instrument + exo)]
      firststagespecboot <- firststagespecboot[, -(instrument +1)]
      firststagespecboot <- as.data.frame(cbind(firststagespecboot, trimdat$x))
      colnames(firststagespecboot)[instrument + exo ] = "x"
      
      ## Step 4
      ##Conduct 2SLS
      specialregboot <- lm(x~., data=firststagespecboot)
      specHatboot <- specialregboot$residuals
      ## Estimate beta
      
      secondstagespecboot <- trimdatboot[, (3):(2 + instrument + exo)]
      firststagespecboot$t <- trimdatboot$t
      firststagespecboot$specHat  <- specHat
      specialreg2boot <- lm(t~x+specHat, data=firststagespecboot)
      bootse[p] <- specialreg2boot$coefficients[2]
      bootse2[p] <- specialreg2boot$coefficients[3]
    }
    return(list(se= sd(bootse), se2  = sd(bootse2)))
  }
  n =1000
  r7 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
  results <- matrix(0, nrow=ncol(datmat), ncol= 1)
  e7 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
  c7 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
  endogeneity <- matrix(0, nrow=ncol(datmat), ncol= 1)
  coverage <- matrix(0, nrow=ncol(datmat), ncol= 1)
  se <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
  se2 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
  t <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
  est7 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
  cover <- function(estimate, se){
    upper <- estimate + 1.96*se
    lower <- estimate - 1.96*se
    if (.5 > lower & .5 < upper){
      return(1)}
    else{
      return(0)}
  }

  test <- function(s){
    if (abs(s) > 1.96){
      return(1)
    }
    else {
      return(0)  
    }
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
      
      xo <- as.data.frame(xo)
      ##demean our special regressor, the exogenous variable "xo". This is our V
      dat$demeanxo <- xo[,1]-mean(xo[,1])
      demeanxo <- dat$demeanxo
      
      ##Obtain residuals, which are our U_i
      specialregdata = as.data.frame(cbind(demeanxo, x, z))
      if (exo > 1){
        specialregdata = cbind(specialregdata, xo[, 2:exo])
      }
      
      fssr <- lm(demeanxo ~., data=specialregdata)
      dat$u <- fssr$residuals
      
      ##Step 2
      fhat <- kde(dat$u,eval.points=dat$u)$estimate
      
      
      ##Step 3
      ## Create vector that will collect our I(Vi > 0)
      dat$idv <- rep(0,nrow(dat))
      ## Collect identity function of v
      dat$idv <- replace(dat$idv,which(dat$demeanxo>=0),1)
      
      ##Obtain Ti; y_values=Di, fi and idv is defined above. 
      dat$t <- (dat$y_values- dat$idv)/fhat
      
      ## Trim extreme values
      lower = .025*n + 1
      upper = .975*n
      
      sortdat <- dat[order(dat$t),]
      trimdat <- sortdat[lower:upper,]
      firststagespec <- trimdat[, (3):(2 + instrument + exo)]
      firststagespec <- firststagespec[, -(instrument +1)]
      firststagespec <- as.data.frame(cbind(firststagespec, trimdat$x))
      colnames(firststagespec)[instrument + exo ] = "x"
      
      ## Step 4
      ##Conduct 2SLS
      specialreg <- lm(x~., data=firststagespec)
      specHat <- specialreg$residuals
      ## Estimate beta
      
      secondstagespec <- trimdat[, (3):(2 + instrument + exo)]
      firststagespec$t <- trimdat$t
      firststagespec$specHat  <- specHat
      specialreg2 <- lm(t~x+specHat, data=firststagespec)
      r7[i, j] <- specialreg2$coefficients[2]
      
      cfse = cfboots()
      c7[i,j] <- cover(estimate = r7[i,j], se=cfse$se)
      se[i,j] <- cfse$se
      se2[i,j] <- cfse$se2 
      est7[i,j] <- specialreg2$coefficients[3]
      
      tstat = est7[i,j]/cfse$se2  
      t[i,j] <- tstat
      e7[i,j] = test(tstat)
    }  
    results[j, 1] <- mean(abs(r7[, j]-0.5))
    endogeneity[j,1] <- sum(e7[,j])
    coverage[j, 1] <- sum(c7[,j])
    results[j, 1] <- mean(abs(r7[, j]-0.5))
    
  }
  return(list(results =results, coverage=coverage, endogeneity=endogeneity ))   
}

mad1 <- specreg(datmat=data1)
setwd("..")
bias <- mad1$results
coverage <- mad1$coverage
endogeneity <- mad1$endogeneity
write.csv(bias, "Data/bias16.csv")
write.csv(coverage, "Data/coverage16.csv")
write.csv(endogeneity, "Data/endogeneity16.csv")