rm(list = ls())
library(ivpack)
library(ggplot2)
library(gridExtra)
setwd("~/Simulations/R")
data1 <- source("DataGen1.R")
data2 <- source("DataGen2.R")
data3 <- source("DataGen3.R")
data4 <- source("DataGen4.R")
data5 <- source("DataGen5.R")
setwd("..")

data2 <- data2$value
data3 <- data3$value
data4 <- data4$value
data5 <- data5$value


data = as.matrix(data1, ncol=5, nrow=1000)
data[,2] <- data2[, 4]
data[,3] <- data3[, 4]
data[,4] <- data4[,4]
data[,5] <- data5[,4]

regressions <- function(datmat, exo=1, instrument=1){
  n = 1000  
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
  
  test <- function(s){
    if (abs(s) > 1.96){
      return(0)
    }
    else {
      return(1)  
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
  r5 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))  
  c5 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat)) 
  se <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
  e5 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat)) 
  est5 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
  t <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
  for (j in 1:ncol(datmat)){  
    for (i in 1:nrow(datmat)){
      dat= unlist(datmat[i,j])
      dat = matrix(dat, ncol=5 ,nrow = 1000)
      
      y_values = dat[,1]
      ypre = dat[,2]
      x <- dat[, 3]
      ##Obtain IV z (excluded exogenous regressor)
      z <- dat[, (4):(3+instrument)]
      ##Obtain included exogenous regressor
      xo <- dat[, (4+instrument):(3+ instrument+exo)]
      
      cover <- function(estimate, se){
        upper <- estimate + 1.96*se
        lower <- estimate - 1.96*se
        if (.5 > lower & .5 < upper){
          return(1)}
        else{
          return(0)}
      }
      
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

        r5[i, j] <- out$estimate[2]
        
        c5[i,j] <- cover(estimate = r5[i,j], se=cfse$est)     
        
        rho <- out$estimate[5]
        se[,j] <- cfse$est
        critval = (rho)^2/cfse$rhovar
        t[i,j] <- critval
        e5[i,j] <- ifelse(critval > 3.841,0, 1)
      
    }      
    
    yuh[[j]] <- as.data.frame(cbind(abs(r5[,j]-.5), se[,j], e5[,j], abs(t[,j]), r5[,j]))
    colnames(yuh[[j]])[1:5] <- cbind("AbsoluteDeviation", "StndErr", "Endogeneity",  "TStat", "Beta")
    
    rho = "Cor(x, z) ="
    
    p3 <- ggplot(yuh[[j]], aes(y=StndErr, x=Beta))+geom_point(aes(color= factor(Endogeneity)), show.legend = ifelse(j==1, TRUE, FALSE))+
      xlim(min(yuh[[j]]$Beta)-.05,max(yuh[[j]]$Beta)+.05)+
      ylim(min(yuh[[j]]$StndErr-.005),max(yuh[[j]]$StndErr+.005))+
      xlab("Beta")+
      ylab("Standard Error")+ labs(color = "Endogeneity\nTests") +scale_fill_manual(values = c("#4bfffd", "#030057"))+
      geom_vline(xintercept = .5)+
      ggtitle(paste(rho, rhoxe[j]))+
      scale_colour_manual(values = c("0" = "#74BDCB", "1" = "#FFA384"))+theme_bw()
    
    plots3[[j]] <- p3
    
    p4 <- ggplot(yuh[[j]], aes(x=Beta, y=TStat))+geom_point(aes(color= factor(Endogeneity)),show.legend = FALSE)+
      xlim(min(yuh[[j]]$Beta)-.05,max(yuh[[j]]$Beta)+.05)+
      ylim(0,max(yuh[[j]]$TStat + .05))+
      ylab("Critical Value")+
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
  
  png("Graphics/Graph5.png")
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
}

sink("NULL")
mad1 <- regressions(datmat=data)
sink()