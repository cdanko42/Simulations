rm(list = ls())

source("DataGen3.R")
library(ks)
library(lmtest)
library(ivpack)
library(ggplot2)
library(gridExtra)


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
      specialregdata = as.data.frame(cbind(demeanbootxo, bootx, bootz))
      if (exo > 1){
        specialregdata = cbind(specialregdata, bootxo[, 2:exo])
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
      bootse[p] <- specialreg2$coefficients[2]
      bootse2[p] <- specialreg2$coefficients[3]
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

  yuh = list()
  plot=list()
  plot2= list()
  plots=list()
  plots2=list()
  plot3 = list()
  plot4=list()
  plots3=list()
  plots4=list()
  rhoxe= c(.1,.3,.5,.7,.9)
  
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
      
      tstat = specialreg2$coefficients[3]/cfse$se2  
      t[i,j] <- tstat
      e7[i,j] = test(tstat)
    }  
    results[j, 1] <- mean(abs(r7[, j]-0.5))
    endogeneity[j,1] <- sum(e7[,j])
    coverage[j, 1] <- sum(c7[,j])
    
    yuh[[j]] <- as.data.frame(cbind(abs(r7[,j]-.5), se[,j], e7[,j], abs(t[,j]), r7[,j], c7[,j], (r7[,j]-.5)/se[,j], est7[,j]))
    colnames(yuh[[j]])[1:8] <- cbind("AbsoluteDeviation", "StndErr", "Endogeneity",  "TStat", "Beta", "Coverage", "SDistance", "Rho")
    
    rho = "Cor(x, z) ="
    
    p3 <- ggplot(yuh[[j]], aes(y=Rho, x=Beta))+geom_point(aes(color= factor(Endogeneity)), show.legend = ifelse(j==1, TRUE, FALSE))+
      xlim(min(yuh[[j]]$Beta)-.05,max(yuh[[j]]$Beta)+.05)+
      ylim(min(yuh[[j]]$Rho-.005),max(yuh[[j]]$Rho+.005))+
      xlab("Beta")+
      ylab("Rho")+ labs(color = "Endogeneity\nTests") +scale_fill_manual(values = c("#4bfffd", "#030057"))+
      geom_vline(xintercept = .5)+
      ggtitle(paste(rho, rhoxe[j]))+
      scale_colour_manual(values = c("0" = "#74BDCB", "1" = "#FFA384"))+theme_bw()
    
    plots3[[j]] <- p3
    
    p4 <- ggplot(yuh[[j]], aes(x=Beta, y=TStat))+geom_point(aes(color= factor(Endogeneity)),show.legend = FALSE)+
      xlim(min(yuh[[j]]$Beta)-.05,max(yuh[[j]]$Beta)+.05)+
      ylim(0,max(yuh[[j]]$TStat + .05))+
      ylab("T-Stat on Rho")+
      xlab("Beta")+
      geom_hline(yintercept = 1.96)+
      ggtitle(paste(rho, rhoxe[j]))+
      scale_colour_manual(values = c("0" = "#74BDCB", "1" = "#FFA384"))+theme_bw()
    plots4[[j]] <- p4
    
    
  }
  plot3[[1]] <- plots3[[1]]
  plot3[[2]] <- plots4[[1]]
  plot3[[3]] <- plots3[[2]]
  plot3[[4]] <- plots4[[2]]
  plot3[[5]] <- plots3[[3]]
  plot3[[6]] <- plots4[[3]]
  plot4[[1]] <- plots3[[4]]
  plot4[[2]] <- plots4[[4]]
  plot4[[3]] <- plots3[[5]]
  plot4[[4]] <- plots4[[5]]
  
  options(bitmapType='cairo')
  
  png("Graphics/Graph15.png")
  grid.arrange(grobs= plot3,
               widths = c(4,1.65,4),
               heights= unit(c(1.8,1.8, 1.8), c("in", "in")),
               top = "Relationship Between Beta and Standard Error",
               layout_matrix=rbind(c(1,1, 2), c(3,NA, 4), c(5,NA, 6)))
  dev.off()
  
  png("Graphics/Graph6.png")
  grid.arrange(grobs= plot4,
               widths = c(4,1.65,4),
               heights= unit(c(1.8,1.8), c("in", "in")),
               layout_matrix=rbind(c(1,NA, 2), c(3,NA,4)))
  dev.off()
  
  return(list(results =results, coverage=coverage, endogeneity=endogeneity ))  

  }

mad1 <- specreg(datmat=data1)
mad1$results
mad1$coverage
mad1$endogeneity