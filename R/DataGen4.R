library(MASS)

set.seed(8675309)

simulateiv <- function(n=1000, size=1000, rhoxz, rhoxe, eevs= 1, exo =1, instrument =1 ){
  ##rho = correlation of instrumental variable with x
  ##  Initialize matrices
  rhoxz <- matrix(rhoxz, nrow=length(rhoxz), ncol=eevs)
  rhoxe <- matrix(rhoxe, nrow= length(rhoxe), ncol =eevs)
  
  initmatrix <- matrix(0, nrow = n, ncol = eevs+exo+instrument+2)
  bigdat <- array(list(initmatrix) , c(size, length(rhoxe)))
  
  innergen <- function(){
    ##Create model vcov matrix
    # Let's add in one continuous exogenous regressor - we'll most surely want this later on
    # For now, there's no reason to assume multicollinearity between endog var and exog var
    sig <- matrix(0, nrow=eevs+exo+instrument+1, ncol = eevs+exo+instrument+1)
    for (k in 1:ncol(rhoxe)){
      sig[1, k+1] = rhoxe[j, k]
      sig[k+1, 1] = rhoxe[j, k]
      sig[k+1, (eevs+exo+2):(eevs+exo+instrument+1)] = rhoxz[, k]
      sig[(eevs+exo+2):(eevs+exo+instrument+1), k+1] = rhoxz[, k]
    }
    ## Add 1 diagonal, initialize vector
    variable = c()
    for (h in 1:(eevs+exo+instrument+1)){
      sig[h,h]=1
      variable[h] = 0
    }
    xez <- mvrnorm(n, variable,Matrix::nearPD(sig, TRUE, TRUE)$mat)
    ##Obtain error, e
    e <- xez[, 1]
    ##Obtain EEV, x
    x <- xez[, 2:(1+eevs)]
    ##Obtain IV z (excluded exogenous regressor)
    z <- xez[, (2+eevs+exo):(1+ eevs+ instrument+ exo)]
    ##Obtain included exogenous regressor
    xo <- xez[, (2+eevs):(1+ eevs+exo)]
    xo <- as.data.frame(xo)
    
    ##Specify and sample from true model
    ypre = 1 + 0.5*x + e
    
    for (g in 1:exo){
      ypre = ypre - (.3/exo)*xo[,g]
    }
    
    xo <- xez[, (2+eevs):(1+ eevs+exo)]
    
    y_values = rep(0,n)
    y_values=replace(y_values,which(ypre>0),1)
    
    dat = as.matrix(cbind(y_values, ypre, x,z,xo))
    
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
    
    start<-c(0,0,0,1,.5)
    start[6:(5+exo)] <- 0
    start[(6+exo):(5+2*exo)] <- 0
    start[(6+2*exo):(5+2*exo+instrument)] <- 0
    
    tryCatch({
      out<-nlm(lik,start,iterlim=1000)
      eststore<-out$estimate
      estgrad<-out$gradient
      convcode<-out$code}, error=function(e){cat("Error:",conditionMessage(e), "\n")})
    return(list(data=dat,code = out$code))
  }
  
  for (j in 1:(nrow(rhoxe))){
    
    for (i in 1:size){
      
      code = 2
      
      tryCatch({
      inner <- innergen()
      code = inner$code
      dat = inner$data}, error=function(e){cat("Error:",conditionMessage(e), "\n")})
      
      while (code != 1){
        inner <- innergen()
        code = inner$code
        dat = inner$data  
      }
      bigdat[i, j] <- list(dat)  
    }
  }
  return(bigdat)
}

data1 = simulateiv(size = 1000, rhoxz = 0.7, rhoxe = c(0,0.1, .2,.3,.4,.5))