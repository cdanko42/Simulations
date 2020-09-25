library(MASS)

set.seed(8675309)
simulateiv <- function(n=1000, size=1000, rhoxz, rhoxe, eevs= 1, exo =1, instrument =1 ){
  ##rho = correlation of instrumental variable with x
  ##  Initialize matrices
  rhoxz <- matrix(rhoxz, nrow=length(rhoxz), ncol=eevs)
  rhoxe <- matrix(rhoxe, nrow= length(rhoxe), ncol =eevs)
  
  initmatrix <- matrix(0, nrow = n, ncol = eevs+exo+instrument+2)
  bigdat <- array(list(initmatrix) , c(size, length(rhoxe)))
  
  for (j in 1:(nrow(rhoxe))){
    
    for (i in 1:size){
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
      
      y_values = rep(0,n)
      y_values=replace(y_values,which(ypre>0),1)
      
      dat = as.matrix(cbind(y_values, ypre, x,z,xo))
      bigdat[i, j] <- list(dat)
      
    }
  }
  return(bigdat)
}

data1 = simulateiv(rhoxz = 0.5, rhoxe = c(.1,.2,0.3,.4,.5))

