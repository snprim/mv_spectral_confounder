library(tidyverse)

# these scripts contain necessary functions 
source("functions_analysis.R")
source("cross_validation.R")

# read in data 

Y_star <- readRDS("Y_star.rds")
X_star <- readRDS("X_star.rds")
Z_star <- readRDS("Z_star.rds")

Qmat_south <- readRDS("Qmat_south.rds")
Gamma_south <- readRDS("Gamma_south.rds")
W_south <- readRDS("W_south.rds")

exposure_names <- c("Theme 1", "Theme 2", "Theme 3", "Theme 4")
outcome_names <- c("Hypertension", "CKD", "Hyperlipidemia", "CHF", "Diabetes")

set.seed(919)

B5_2 <-  bSpline(rank(W_south)/length(W_south),df=5, intercept = T)
B10_2 <-  bSpline(rank(W_south)/length(W_south),df=10, intercept = T)

# use standardized basis function
# exposures are 4 svi
# covariates are log(population) and urbanicity
# final model
fit_svi_K10_L10_stB <- MCMC_model4_8_3_XZ(Y_star, X_star, Z_star, Qmat_south, E=4, n=nrow(Y_star), R=5, Q=5, L=10, K=10, B=B10_2, iters = 100000, burn = 10000, W_south, thin = 10)


### naive 

# final naive model
fit_svi_naive <- MCMC_model3_2_XZ(Y_star, X_star, Z_star, Qmat_south, p=4, n=nrow(Y_star), R=5, Q=5, iters=100000, burn=10000, 1, nrow(Y_star), W_south, thin = 10)

