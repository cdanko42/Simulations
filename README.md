# Simulations

## Data Generation

The R folder contains all of the code required to produce both the graphics and the results. DataGen1-5 simulates data sets with corresponding correlations. For example, DataGen1.R produces five sets of data, with cor(x,err) = (.1, .3, .5, .7, .9) respectively. For all five data sets, cor(x,z) = .1. DataGen2.R produces five sets of data with cor(x, err) identitical to DataGen1.R and cor(x,z) =.3. For DataGen3-5, cor(x,z) = (.5, .7, .9), respectively. 

## Methods

Regressions1-5 produces the estimates for the OLS, 2SLS, Probit, and 2SLPM models on the simulated data. Regressions1 uses the data sets produced by DataGen1. Regressions2 uses the data sets produced by DataGen2, and so on. CF1-5 produces the estimates for the Special Regressor model on the simulated data. Note that the simulated data is exactly the same across methods. Probit1-5 produces the estimates for the Control Function model on the simulated data. Max1-5 produces the Maximum Likelihood estimates for the simulated data.

## Results
regout1-5 contains the result from Regressions1-5. cfout1-5 contains the result from CF1-5. probout1-5 contains the results from Probit1-5. maxout1-5 contains the results from Max1-5. 

## Graphics
GraphicsGen1-2 produces Graph1 and Graph2, found in the graphics folder. GraphicsGen3-4 produces Graph3 and Graph4. 
