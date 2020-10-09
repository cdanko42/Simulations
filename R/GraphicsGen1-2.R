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
  se2 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
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
      
      ivyval <- ivreg(y_values~x+ xo, ~z + xo)
      r5[i, j] <- ivyval$coefficients[2]
      invisible(iv2se <- robust.se(ivyval)[2,2])
      c5[i,j] <- cover(estimate = r5[i,j], se=iv2se)
      
      se[i,j] <- iv2se
      
      firststage <- (lm(x~z+xo))$residuals
      secondstep <- lm(y_values~x+xo +firststage)
      est5[i,j] <-summary(secondstep)$coefficients[4,1]
      se2[i,j] <-summary(secondstep)$coefficients[4,2]
      s <- summary(secondstep)$coefficients[4,3]
      t[i,j] <- abs(s)
      e5[i,j] <- test(s=s)
      
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
      ylab("T-Stat")+
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

png("Graphics/Graph1.png")
  grid.arrange(grobs= plot3,
               widths = c(4,1.65,4),
               heights= unit(c(1.8,1.8, 1.8), c("in", "in")),
               top = "Relationship Between Beta and Standard Error",
               layout_matrix=rbind(c(1,1, 2), c(3,NA, 4), c(5,NA, 6)))
  dev.off()
  
  png("Graphics/Graph2.png")
  grid.arrange(grobs= plot4,
               widths = c(4,1.65,4),
               heights= unit(c(1.8,1.8), c("in", "in")),
               layout_matrix=rbind(c(1,NA, 2), c(3,NA,4)))
  dev.off()
}

sink("NULL")
mad1 <- regressions(datmat=data)
sink()