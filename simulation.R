
#############################################
##### Author: Shih-Ni Prim ##################
#############################################

# This is the final simulations 
# 100 datasets for main simulation
# 20 datasets for sensitivity test

# MCMC_model4_8 full model with HS prior
# MCMC_model4_8_2 full model with IG prior
# MCMC_model4_9 univariate with HS prior
# MCMC_model4_9_2 univariate with IG prior
# MCMC_model3_2 does not allow coefficients to vary by frequency

# all settings with lambda = .9 and rho = .9
# vary bXZ and bw between 1 and 2 (4 combinations)

# model 1: full model with HS prior, L = 10
# model 2: full model with HS prior, L = 20
# model 3: full model with IG prior, L = 10
# model 4: full model with IG prior, L = 20
# model 5: coefficients do not vary by frequency
# model 6: one step spatial plus 1:320
# model 7: one step spatial plus 1:160
# model 8: univariate model with HS prior, L = 10
# use old basis functions that are uncorrelated 

args = commandArgs(TRUE)
ind = as.numeric(args[1])

source("functions_sim.R")

iters <- 5000
burn <- 1000
sig2X <- c(3,3,3,3,3,3,3,3,3)
sig2Z <- rep(16,10)
sig2theta <- c(2,2,2,2,2)
tau2r <- c(4,4,4,4,4)
lambdaZ <- .9999
lambda1 <- .9
lambda2 <- .9
rho1 <- .9
rho2 <- .9
n <- 400
n1 <- 20
E <- 10
R <- 5
nz <- 10

L <- 10
K <- 5

Qmat <- makeQ(20)
W <- eigen(Qmat)$values

beta <- matrix(c(3,-3,2,3,5,4,-2,4,5,7,2,-3,3,4,6,1,-4,2,3,4,5,1,3,4,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), nrow = E, ncol = R, byrow = T)/10

reps <- 20

load("simulation/inputs/B1.RData")
load("simulation/inputs/B2.RData")

###########################################################
################ main simulation ##########################
###########################################################

# extra smoothing, bw = 1, bXZ = 1
# B1 and B2
if (ind == 1){
  sim1_par <- foreach (i=1:reps, .combine = rbind) %dopar% {
    list(sim_func4_1_tuneL_6(n, nz, E, R, L=10, K=5, Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters, burn, num1 = 1, num2 = 320, n1, bXZ = 1, bw = 1, i))
  }
  sim1 <- flatten(sim1_par, reps, 8, n, 7, E, R)
  save(sim1, file = "sim1.RData")
}

# extra smoothing, bw = 1, bXZ = 2
# B1 and B2
if (ind == 2){
  sim2_par <- foreach (i=1:reps, .combine = rbind) %dopar% {
    list(sim_func4_1_tuneL_6(n, nz, E, R, L=10, K=5, Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters, burn, num1 = 1, num2 = 320, n1, bXZ = 2, bw = 1, i))
  }
  sim2 <- flatten(sim2_par, reps, 8, n, 7, E, R)
  save(sim2, file = "sim2.RData")
}

# extra smoothing, bw = 2, bXZ = 1
# B1 and B2
if (ind == 3){
  sim3_par <- foreach (i=1:reps, .combine = rbind) %dopar% {
    list(sim_func4_1_tuneL_6(n, nz, E, R, L=10, K=5, Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters, burn, num1 = 1, num2 = 320, n1, bXZ = 1, bw = 2, i))
  }
  sim3 <- flatten(sim3_par, reps, 8, n, 7, E, R)
  save(sim3, file = "sim3_betaL.RData")
}

# extra smoothing, bw = 2, bXZ = 2
# B1 and B2
if (ind == 4){
  sim4_par <- foreach (i=1:reps, .combine = rbind) %dopar% {
    list(sim_func4_1_tuneL_6(n, nz, E, R, L=10, K=5, Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters, burn, num1 = 1, num2 = 320, n1, bXZ = 2, bw = 2, i))
  }
  sim4 <- flatten(sim4_par, reps, 8, n, 7, E, R)
  save(sim4, file = "sim4.RData")
}

###############################################################
############### sensitivity tests #############################
###############################################################

reps <- 20

Lvec <- c(5, 10, 20)
Kvec <- c(2, 5, 10)

# test combinations of L and K
# extra smoothing, bw = 1, bXZ = 2
# B1 and B2

if (ind == 1){
  sim_par <- foreach (i=1:reps, .combine = rbind) %dopar% {
    list(sim_func4_1_tuneL_6(n, nz, E, R, L=Lvec[1], K=Kvec[1], Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters, burn, num1 = 1, num2 = 320, n1, bXZ = 2, bw = 1, i))
  }
  sim <- flatten(sim_par, reps, 8, n, 7, E, R)
  name <- paste0("sim_L",Lvec[1],"_K",Kvec[1],".RData")
  save(sim, file = name)
}

if (ind == 2){
  sim_par <- foreach (i=1:reps, .combine = rbind) %dopar% {
    list(sim_func4_1_tuneL_6(n, nz, E, R, L=Lvec[1], K=Kvec[2], Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters, burn, num1 = 1, num2 = 320, n1, bXZ = 2, bw = 1, i))
  }
  sim <- flatten(sim_par, reps, 8, n, 7, E, R)
  name <- paste0("sim_L",Lvec[1],"_K",Kvec[2],".RData")
  save(sim, file = name)
}

if (ind == 3){
  sim_par <- foreach (i=1:reps, .combine = rbind) %dopar% {
    list(sim_func4_1_tuneL_6(n, nz, E, R, L=Lvec[1], K=Kvec[3], Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters, burn, num1 = 1, num2 = 320, n1, bXZ = 2, bw = 1, i))
  }
  sim <- flatten(sim_par, reps, 8, n, 7, E, R)
  name <- paste0("sim_L",Lvec[1],"_K",Kvec[3],".RData")
  save(sim, file = name)
}

if (ind == 4){
  sim_par <- foreach (i=1:reps, .combine = rbind) %dopar% {
    list(sim_func4_1_tuneL_6(n, nz, E, R, L=Lvec[2], K=Kvec[1], Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters, burn, num1 = 1, num2 = 320, n1, bXZ = 2, bw = 1, i))
  }
  sim <- flatten(sim_par, reps, 8, n, 7, E, R)
  name <- paste0("sim_L",Lvec[2],"_K",Kvec[1],".RData")
  save(sim, file = name)
}

if (ind == 5){
  sim_par <- foreach (i=1:reps, .combine = rbind) %dopar% {
    list(sim_func4_1_tuneL_6(n, nz, E, R, L=Lvec[2], K=Kvec[2], Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters, burn, num1 = 1, num2 = 320, n1, bXZ = 2, bw = 1, i))
  }
  sim <- flatten(sim_par, reps, 8, n, 7, E, R)
  name <- paste0("sim_L",Lvec[2],"_K",Kvec[2],".RData")
  save(sim, file = name)
}

if (ind == 6){
  sim_par <- foreach (i=1:reps, .combine = rbind) %dopar% {
    list(sim_func4_1_tuneL_6(n, nz, E, R, L=Lvec[2], K=Kvec[3], Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters, burn, num1 = 1, num2 = 320, n1, bXZ = 2, bw = 1, i))
  }
  sim <- flatten(sim_par, reps, 8, n, 7, E, R)
  name <- paste0("sim_L",Lvec[2],"_K",Kvec[3],".RData")
  save(sim, file = name)
}

if (ind == 7){
  sim_par <- foreach (i=1:reps, .combine = rbind) %dopar% {
    list(sim_func4_1_tuneL_6(n, nz, E, R, L=Lvec[3], K=Kvec[1], Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters, burn, num1 = 1, num2 = 320, n1, bXZ = 2, bw = 1, i))
  }
  sim <- flatten(sim_par, reps, 8, n, 7, E, R)
  name <- paste0("sim_L",Lvec[3],"_K",Kvec[1],".RData")
  save(sim, file = name)
}

if (ind == 8){
  sim_par <- foreach (i=1:reps, .combine = rbind) %dopar% {
    list(sim_func4_1_tuneL_6(n, nz, E, R, L=Lvec[3], K=Kvec[2], Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters, burn, num1 = 1, num2 = 320, n1, bXZ = 2, bw = 1, i))
  }
  sim <- flatten(sim_par, reps, 8, n, 7, E, R)
  name <- paste0("sim_L",Lvec[3],"_K",Kvec[2],".RData")
  save(sim, file = name)
}

if (ind == 9){
  sim_par <- foreach (i=1:reps, .combine = rbind) %dopar% {
    list(sim_func4_1_tuneL_6(n, nz, E, R, L=Lvec[3], K=Kvec[3], Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, beta, B1, B2, rho1, rho2, iters, burn, num1 = 1, num2 = 320, n1, bXZ = 2, bw = 1, i))
  }
  sim <- flatten(sim_par, reps, 8, n, 7, E, R)
  name <- paste0("sim_L",Lvec[3],"_K",Kvec[3],".RData")
  save(sim, file = name)
}

