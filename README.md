# mv_spectral_confounder
This repository contains code and plots for "A spectral confounder adjustment for spatial regression with multiple exposures and outcomes." The paper can be found on [arXiv](https://arxiv.org/abs/2506.09325).

## Overview 

The scientific data management plan is described in the Word document `Prim_SDMP.docx`. All relevant code, except for the health outcome data that is proprietary, is included in this repository.

## Simulation code

R functions for simulation are included in `functions_sim.R`. The code to run simulations is included in `simulation.R`. The two matrices used to generate data, `B1` and `B2`, are placed in `imulatin/inputs`.

To test the function, below is an example:

```{r}
source("functions_sim.R")
Qmat <- makeQ(20)
EI <- eigen(Qmat)
Gamma <- EI$vectors
W <- EI$values
B <- bSpline(W,df=10,intercept=T)
load("simulation/inputs/B1.RData")
load("simulation/inputs/B2.RData")
beta <- matrix(c(3,-3,2,3,5,4,-2,4,5,7,2,-3,3,4,6,1,-4,2,3,4,5,1,3,4,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), nrow = 10, ncol = 5, byrow = T)/10
set <- create_set7(n=400, nz=10, E=10, R=5, Q=Qmat, sig2X=c(3,3,3,3,3,3,3,3,3), sig2Z=rep(16,10), sig2theta=c(2,2,2,2,2), tau2r=c(4,4,4,4,4), lambda0=.9999, lambda1=.9, lambda2=.9, betas=beta, B1=B1, B2=B2, rho1=.9, rho2=.9, n1=20, bXZ=1, bw=1)
Ystar <- t(Gamma)%*%set$Y
Xstar <- t(Gamma)%*%set$X
fit <- MCMC_model4_8_v7(Ystar, Xstar, Qmat = Qmat, E = 10, n = 400, R = 5, Q = 5, L=10, K=5, B = B, iters = 5000, burn = 1000, full_result = T, W=W)
```


## Data analysis code

Code preparing data is in `data_preparation.R`, code for running data analysis is in `final_alayses.R`, and code for making plots is in `final_results_plots.R`. SVI data can be downloaded using the R package `findSVI`, and the health outcome data from Medicare is not provided here. 

## Citation 

The arXiv paper can be cited as:

Prim, S.-N., Guan, Y., Yang, S., Rappold, A. G., Hill, K. L., Tsai, W.-L., Keeler, C. and Reich, B. J. (2025) A Spectral Confounder Adjustment for Spatial Regression with Multiple Exposures and Outcomes. *arXiv*:2506.09325.


