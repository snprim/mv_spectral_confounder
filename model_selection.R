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

####################################
######## cross validation ##########
####################################

# set training set 
cv_num <- 5
set.seed(919)
ind.test <- sample(1:cv_num, nrow(Y_star), replace = TRUE)
Y_train_lst <- list()
for (i in 1:cv_num){
  Y_train_lst[[i]] <- Y_star
  for (r in 1:ncol(Y_star)){
    Y_train_lst[[i]][,r] <- ifelse(ind.test==i, NA, Y_star[,r])  
  }
}

# only do the first fold

# K=5, L=5
fit_svi_K5_L5_cv <- list()
log_score_svi_K5_L5_cv <- rep(NA, cv_num)
for (i in 1:1){
  fit_svi_K5_L5_cv[[i]] <- MCMC_model4_8_3_XZ_cv(Y_train_lst[[i]], X_star, Z_star, Qmat_south, E=4, n=nrow(Y_star), R=5, Q=5, L=5, K=5, B=B5_2, iters = 20000, burn = 2000, W_south, thin = 10)
  miss <- which(is.na(Y_train_lst[[i]][,1]))
  log_score <- rep(NA, ncol(Y_star)) 
  for (r in 1:ncol(Y_star)){
    log_score[r] <- mean(dnorm(Y_star[miss,r], fit_svi_K5_L5_cv[[i]]$y_mean[miss,r], sqrt(fit_svi_K5_L5_cv[[i]]$y_var[miss,r]), log = TRUE))
  }
  log_score_svi_K5_L5_cv[i] <- mean(log_score)
}

# K=5, L=10
fit_svi_K5_L10_cv <- list()
log_score_svi_K5_L10_cv <- rep(NA, cv_num)
for (i in 1:1){
  fit_svi_K5_L10_cv[[i]] <- MCMC_model4_8_3_XZ_cv(Y_train_lst[[i]], X_star, Z_star, Qmat_south, E=4, n=nrow(Y_star), R=5, Q=5, L=10, K=5, B=B10_2, iters = 20000, burn = 2000, W_south, thin = 2)
  miss <- which(is.na(Y_train_lst[[i]][,1]))
  log_score <- rep(NA, ncol(Y_star)) 
  for (r in 1:ncol(Y_star)){
    log_score[r] <- mean(dnorm(Y_star[miss,r], fit_svi_K5_L10_cv[[i]]$y_mean[miss,r], sqrt(fit_svi_K5_L10_cv[[i]]$y_var[miss,r]), log = TRUE))
  }
  log_score_svi_K5_L10_cv[i] <- mean(log_score)
}

# K=10, L=5
fit_svi_K10_L5_cv <- list()
log_score_svi_K10_L5_cv <- rep(NA, cv_num)
for (i in 1:1){
  fit_svi_K10_L5_cv[[i]] <- MCMC_model4_8_3_XZ_cv(Y_train_lst[[i]], X_star, Z_star, Qmat_south, E=4, n=nrow(Y_star), R=5, Q=5, L=5, K=10, B=B5_2, iters = 20000, burn = 2000, W_south, thin = 10)
  miss <- which(is.na(Y_train_lst[[i]][,1]))
  log_score <- rep(NA, ncol(Y_star)) 
  for (r in 1:ncol(Y_star)){
    log_score[r] <- mean(dnorm(Y_star[miss,r], fit_svi_K10_L5_cv[[i]]$y_mean[miss,r], sqrt(fit_svi_K10_L5_cv[[i]]$y_var[miss,r]), log = TRUE))
  }
  log_score_svi_K10_L5_cv[i] <- mean(log_score)
}

# K=10, L=10
fit_svi_K10_L10_cv <- list()
log_score_svi_K10_L10_cv <- rep(NA, cv_num)
for (i in 1:1){
  fit_svi_K10_L10_cv[[i]] <- MCMC_model4_8_3_XZ_cv(Y_train_lst[[i]], X_star, Z_star, Qmat_south, E=4, n=nrow(Y_star), R=5, Q=5, L=10, K=10, B=B10_2, iters = 20000, burn = 2000, W_south, thin = 2)
  miss <- which(is.na(Y_train_lst[[i]][,1]))
  log_score <- rep(NA, ncol(Y_star)) 
  for (r in 1:ncol(Y_star)){
    log_score[r] <- mean(dnorm(Y_star[miss,r], fit_svi_K10_L10_cv[[i]]$y_mean[miss,r], sqrt(fit_svi_K10_L10_cv[[i]]$y_var[miss,r]), log = TRUE))
  }
  log_score_svi_K10_L10_cv[i] <- mean(log_score)
}

# check convergence of cross validation 

# calculate beta from Tl, Te, Tr
n <- nrow(Y_star)
Test_beta_cv <- array(0, c(n,E,R,1800))
model_cv <- fit_svi_K10_L10_cv[[1]]
for (i in 1:1800){
  Test_beta_cv[,,,i] <- bbygamma(B10_2, find_tol_gamma(model_cv$Tl[,,i], model_cv$Te[,,i], model_cv$Tr[,,i]), n) 
}

### check convergence only for cross validation ###
mean_cv <-  rowMeans(Test_beta_cv[1,,,], dims = 2)
upper_cv <- lower_cv <- array(NA, c(E,R))
prob_cv <- array(NA, c(E,R))

for (e in 1:E){
  for (r in 1:R){
    lower_cv[e,r] <- quantile(Test_beta_cv[1,e,r,], probs = .025)
    upper_cv[e,r] <- quantile(Test_beta_cv[1,e,r,], probs = .975)
    prob_cv[e,r] <- mean(Test_beta_cv[1,e,r,]>0)
  }
}

rownames(prob_cv) <- exposure_names
colnames(prob_cv) <- outcome_names
round(prob_cv,2)

par(mfrow = c(E,R), mar = c(2,2,2,2))
for (e in 1:E){
  for (r in 1:R){
    plot(Test_beta_cv[1,e,r,], type = "l", main = paste0(exposure_names[e],", ", outcome_names[r]))
    abline(h = mean_cv[e,r], col = "yellow")
    abline(h = upper_cv[e,r], col = "red")
    abline(h = lower_cv[e,r], col = "red")
    abline(h = 0, col = "green")
  }
}

