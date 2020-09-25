rm(list = ls())

source("R/DataGen3.R")

library(lmtest)
library(ivpack)

regressions <- function(datmat, exo=1, instrument=1){
r1 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
r2 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
r3 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
r4 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
r5 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))  
results <- matrix(0, nrow=ncol(datmat), ncol= 5)
  
c1 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
c2 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
c3 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
c4 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat))
c5 <- matrix(0, nrow=nrow(datmat), ncol= ncol(datmat)) 
coverage <- matrix(0, nrow=ncol(datmat), ncol= 5)
  
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

olspyre <- lm(ypre ~ x + xo)
r1[i, j] <- olspyre$coefficients[2]
cols <- coeftest(olspyre)[2, 2]
cover <- function(estimate, se){
  upper <- estimate + 1.96*se
  lower <- estimate - 1.96*se
  if (.5 > lower & .5 < upper){
    return(1)}
  else{
    return(0)}
}

c1[i, j] <- cover(estimate= r1[i,j], se = cols)

ivpre <- ivreg(ypre~x+ xo, ~z + xo)
r2[i,j] <- ivpre$coefficients[2]
invisible(ivse <- robust.se(ivpre)[2,2])
c2[i, j] <- cover(estimate = r2[i,j], se=ivse)

yvaldata = as.data.frame(cbind(y_values, x, xo))
olsyval <- lm(y_values ~., data=yvaldata)
r3[i, j] <- olsyval$coefficients[2]
invisible(cols3 <- coeftest(olsyval)[2, 2])
c3[i, j] <- cover(estimate = r3[i,j], se=cols3)

dat = as.data.frame(cbind(y_values, x,z,xo))
probyval <- glm(y_values ~., family = binomial(link = "probit"), data = yvaldata)
r4[i, j] <- probyval$coefficients[2]
invisible(seprobit <- coeftest(probyval)[2,2])
c4[i, j] <- cover(estimate = r4[i,j], se=seprobit)

ivyval <- ivreg(y_values~x+ xo, ~z + xo)
r5[i, j] <- ivyval$coefficients[2]
invisible(iv2se <- robust.se(ivyval)[2,2])
c5[i,j] <- cover(estimate = r5[i,j], se=iv2se)


  
}  
results[j, 1] <- mean(abs(r1[, j]-0.5))
results[j, 2] <- mean(abs(r2[, j]-0.5))  
results[j, 3] <- mean(abs(r3[, j]-0.5))  
results[j, 4] <- mean(abs(r4[, j]-0.5))  
results[j, 5] <- mean(abs(r5[, j]-0.5))  

coverage[j, 1] <- sum(c1[,j])
coverage[j, 2] <- sum(c2[,j])  
coverage[j, 3] <- sum(c3[,j])  
coverage[j, 4] <- sum(c4[,j])
coverage[j, 5] <- sum(c5[,j])  
  
}
return(list(results =results, coverage=coverage ))  
}

mad1 <- regressions(datmat=data1)
mad1$results
mad1$coverage
