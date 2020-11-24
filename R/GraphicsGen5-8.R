rm(list = ls())
library(ivpack)
library(ggplot2)
library(gridExtra)
setwd("~/Simulations/R")
data1 <- source("DataGen1.R")
data1 <- data1$value
data = as.matrix(data1, ncol=5, nrow=1000)
data2 <- source("DataGen2.R")
data3 <- source("DataGen3.R")
data4 <- source("DataGen4.R")
data5 <- source("DataGen5.R")
setwd("..")

data2 <- data2$value
data3 <- data3$value
data4 <- data4$value
data5 <- data5$value

data[,2] <- data2[, 4]
data[,3] <- data3[, 4]
data[,4] <- data4[,4]
data[,5] <- data5[,4]

regressions <- function(datmat, exo=1, instrument=1){
  probitboots <- function(bootsize=599){
    bootse <- c()
    vse <- c()
    for (p in 1:bootsize){
      rows = sample(1:1000, 1000, replace = TRUE)
      bootdat <- dat2[rows, ]
      booty <- as.matrix(bootdat[, 1], ncol=1, nrow=n)
      bootx <- as.matrix(bootdat[, 2], ncol=1, nrow=n)
      bootxo <- as.matrix(bootdat[, (3+instrument):(2+exo+instrument)], ncol=exo, nrow=n)
      bootz <- as.matrix(bootdat[, 3:(2+instrument)], ncol=instrument, nrow=n)
      bootcf1 <- lm(bootx~bootz+bootxo)
      bootv = bootcf1$residuals
      bootcf2 <- glm(booty ~ bootx + bootxo + bootv , family = binomial(link = "probit"), control=list(epsilon = 1e-8, maxit = 100))
      bootse[p] <- bootcf2$coefficients[2]
      vse[p] <- bootcf2$coefficients[4]
    }
    return(list(bootse =sd(bootse), vse = sd(vse)))
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
      
      dat2 =as.data.frame(cbind(y_values, x,z,xo))
      
      probitcf1 <- lm(x~z+xo)
      ## Collect residuals
      v = probitcf1$residuals
      cfdata = as.data.frame(cbind(y_values, x, xo, v))
      probitcf2 <- glm(y_values ~. , family = binomial(link = "probit"), data = cfdata, control=list(epsilon = 1e-8, maxit = 100))
      r5[i, j] <- probitcf2$coefficients[2]
      
      cprobit <- probitboots()
      c5[i,j] <- cover(estimate = r5[i,j], se=cprobit$bootse)      
      se[i,j] <- cprobit$bootse
      
      tstat = probitcf2$coefficients[4]/cprobit$vse
      est5[i,j] = probitcf2$coefficients[4]
      t[i,j] <- tstat
      e5[i,j] = test(tstat)
      
    }  
    
    
    yuh[[j]] <- as.data.frame(cbind(abs(r5[,j]-.5), se[,j], e5[,j], abs(t[,j]), r5[,j], c5[,j], (r5[,j]-.5)/se[,j], est5[,j]))
    colnames(yuh[[j]])[1:8] <- cbind("AbsoluteDeviation", "StndErr", "Endogeneity",  "TStat", "Beta", "Coverage", "SDistance", "Rho")
    
    yuh[[j]]$Endogeneity[which(yuh[[j]]$Endogeneity == 1)] <- "Fail to Reject"
    yuh[[j]]$Endogeneity[which(yuh[[j]]$Endogeneity == 0)] <- "Reject"

    rho = "Cor(x, z) ="
    
    p <- ggplot(yuh[[j]], aes(y=StndErr, x=Beta))+geom_point(aes(color= factor(Coverage)), show.legend = ifelse(j==1, TRUE, FALSE))+
      xlim(min(yuh[[j]]$Beta-.01),max(yuh[[j]]$Beta+.01))+
      ylim(min(yuh[[j]]$StndErr-.005),max(yuh[[j]]$StndErr+.005))+
      xlab("Beta")+
      ylab("Standard Error")+ labs(color = "Endogeneity\nTests") +scale_fill_manual(values = c("#4bfffd", "#030057"))+
      geom_vline(xintercept = mean(yuh[[j]]$AbsoluteDeviation))+
      ggtitle(paste(rho, rhoxe[j]))+
      scale_colour_manual(values = c("0" = "#74BDCB", "1" = "#FFA384"))+theme_bw()
    
    plots[[j]] <- p
    
    p2 <- ggplot(yuh[[j]], aes(x=Beta, y=SDistance))+geom_point(aes(color= factor(Coverage)),show.legend = FALSE)+
      xlim(min(yuh[[j]]$Beta-.05),max(yuh[[j]]$Beta+.05))+
      ylim(min(yuh[[j]]$SDistance - .05),max(yuh[[j]]$SDistance + .05))+
      ylab("T-Stat")+
      xlab("Beta")+
      geom_hline(yintercept = c(1.96, -1.96))+
      ggtitle(paste(rho, rhoxe[j]))+
      scale_colour_manual(values = c("0" = "#74BDCB", "1" = "#FFA384"))+theme_bw()
    plots2[[j]] <- p2
    
    p3 <- ggplot(yuh[[j]], aes(y=Rho, x=Beta))+geom_point(aes(color= factor(Endogeneity)), show.legend = ifelse(j==1, TRUE, FALSE))+
      xlim(-3,3)+
      ylim(-3,3.5)+
      xlab("Beta")+
      ylab("Rho")+ labs(color = "Endogeneity\nTests") +scale_fill_manual(values = c("#4bfffd", "#030057"))+
      geom_vline(xintercept = .5)+
      ggtitle(paste(rho, rhoxe[j]))+
      scale_colour_manual(values = c("Reject" = "#74BDCB", "Fail to Reject" = "#FFA384"))+theme_bw()
    
    plots3[[j]] <- p3
    
    p4 <- ggplot(yuh[[j]], aes(x=Beta, y=TStat))+geom_point(aes(color= factor(Endogeneity)),show.legend = FALSE)+
      xlim(-3,3)+
      ylim(0,12)+
      ylab("T-Stat on Rho")+
      xlab("Beta")+
      geom_hline(yintercept = 1.96)+
      ggtitle(paste(rho, rhoxe[j]))+
      scale_colour_manual(values = c("Reject" = "#74BDCB", "Fail to Reject" = "#FFA384"))+theme_bw()
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
  
  plot[[1]] <- plots[[1]]
  plot[[2]] <- plots2[[1]]
  plot[[3]] <- plots[[2]]
  plot[[4]] <- plots2[[2]]
  plot[[5]] <- plots[[3]]
  plot[[6]] <- plots2[[3]]
  plot2[[1]] <- plots[[4]]
  plot2[[2]] <- plots2[[4]]
  plot2[[3]] <- plots[[5]]
  plot2[[4]] <- plots2[[5]]
  
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
  grid.arrange(grobs= plot,
               widths = c(4,1.65,4),
               heights= unit(c(1.8,1.8, 1.8), c("in", "in")),
               top = "Relationship Between Bias and Standard Error",
               layout_matrix=rbind(c(1,1, 2), c(3,NA, 4), c(5,NA, 6)))
  dev.off()
  
  png("Graphics/Graph6.png")
  grid.arrange(grobs= plot2,
               widths = c(4,1.65,4),
               heights= unit(c(1.8,1.8), c("in", "in")),
               layout_matrix=rbind(c(1,NA, 2), c(3,NA,4)))
  dev.off()
  
  png("Graphics/Graph7.png")
  grid.arrange(grobs= plot3,
               widths = c(4,1.65,4),
               heights= unit(c(1.8,1.8, 1.8), c("in", "in")),
               top = "Relationship Between Beta and Rho",
               layout_matrix=rbind(c(1,1, 2), c(3,NA, 4), c(5,NA, 6)))
  dev.off()
  
  png("Graphics/Graph8.png")
  grid.arrange(grobs= plot4,
               widths = c(4,1.65,4),
               heights= unit(c(1.8,1.8), c("in", "in")),
               layout_matrix=rbind(c(1,NA, 2), c(3,NA,4)))
  dev.off()
}

sink("NULL")
mad1 <- regressions(datmat=data)
sink()