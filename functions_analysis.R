library(MCMCpack)
library(MASS)
library(splines2)
library(plyr)
library(LaplacesDemon)
library(doParallel)
library(fields)
library(spam)
registerDoParallel(cores=20)

########################################################
###### Functions for MCMC ##############################
########################################################

xbybeta3 <- function(X, beta, R){
  test <- (X*beta)
  estimate <- matrix(apply(test, FUN = function(X){rowSums(X)}, MARGIN = 3), ncol = R)
}
#  a function to multiply B (nxL) and gamma (LxExR) to get an nxExR array
bbygamma <- function(B, gamma, n){
  B2 <- matrix(B, nrow = n)
  E <- dim(gamma)[2]
  R <- dim(gamma)[3]
  Blist <- lapply(seq_len(R), function(X) B2)
  gammalist <- alply(gamma, 3)
  Bgamma <- mapply("%*%", Blist, gammalist, SIMPLIFY = F)
  Bgamma2 <- array(unlist(Bgamma), c(n, E, R))
  return(Bgamma2)
}
# a function to turn factor matrices (Tl, Te, Tr) into the gamma tensor (LxExR)
find_tol_gamma <- function(Tl, Te, Tr){
  Tllist <- alply(Tl, 2)
  Telist <- alply(Te, 2)
  Trlist <- alply(Tr, 2)
  tol_gamma <- mapply(function(x,y,z) outer(outer(x,y),z), Tllist, Telist, Trlist, SIMPLIFY = F)
  tol_gamma2 <- Reduce("+", tol_gamma)
  return(tol_gamma2)
}

# a function to turn factor matrices (Tl, Te, Tr) into the gamma tensor (LxExR) for univariate
find_tol_gamma2 <- function(Tl, Te, Tr){
  Tllist <- alply(Tl, 2)
  Telist <- alply(Te, 2)
  Trlist <- alply(Tr, 1)
  tol_gamma <- mapply(function(x,y,z) outer(outer(x,y),z), Tllist, Telist, Trlist, SIMPLIFY = F)
  tol_gamma2 <- Reduce("+", tol_gamma)
  return(tol_gamma2)
}

# for univariate
find_tol_gamma3 <- function(Tl, Te){
  Tllist <- alply(Tl, 2)
  Telist <- alply(Te, 2)
  tol_gamma <- mapply(function(x,y) outer(x,y), Tllist, Telist, SIMPLIFY = F)
  tol_gamma2 <- Reduce("+", tol_gamma)
  return(tol_gamma2)
}

find_estimate <- function(X_tensor, B, Tl, Te, Tr, n, R){
  tol_gamma <- find_tol_gamma(Tl, Te, Tr)
  Bgamma <- bbygamma(B, tol_gamma, n)
  estimate <- xbybeta3(X_tensor, Bgamma, R)
}

# for univariate
find_estimate2 <- function(X_tensor, B, Tl, Te, Tr, n, R){
  tol_gamma <- find_tol_gamma2(Tl, Te, Tr)
  Bgamma <- bbygamma(B, tol_gamma, n)
  estimate <- xbybeta3(X_tensor, Bgamma, R)
}

# for univariate
find_estimate3 <- function(X, B, Tl, Te){
  tol_gamma <- find_tol_gamma3(Tl, Te)
  # Bgamma <- bbygamma(B, tol_gamma, n)
  Bgamma <- B %*% tol_gamma
  # estimate <- xbybeta3(X_tensor, Bgamma, R)
  estimate <- rowSums( X * Bgamma )
}

makeQ <- function(m){
  s   <- expand.grid(1:m,1:m)
  A   <- ifelse(as.matrix(dist(s))==1,1,0)
  M   <- diag(rowSums(A))
  Q   <- M - A
  return(Q)
}
# make covariance matrix for CAR(sigma2, lambda)
makecov <- function(sigma2, lambda, Q){
  I <- diag(ncol(Q))
  cov <- sigma2*solve((1-lambda)*I+lambda*Q)
  return(cov)
}
# make AR1 structure
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  rho^exponent
}

# calculate multivariate normal log-likelihood
log_like <- function(K, Kinv, sig2_q){
  -0.5*sum(log(K))-1/(2*sig2_q)*Kinv
}
# make the covariance matrix of Y
makeK <- function(sig2_q, lambda, W){
  sig2_q/(1-lambda+lambda*W)
}
# make the inverse covariance matrix of Y
makeKinv <- function(Uq, lambda, W){
  sum(Uq^2*(1-lambda+lambda*W))
}

#########################################################
####### Functions for data generation ###################
#########################################################

# generate data in spatial domain 
# E means number of covariate plus the intercept
# This version adds extra smoothing on Z
create_set7 <- function(n, nz, E, R, Q, sig2X, sig2Z, sig2theta, tau2r, lambda0, lambda1, lambda2, betas, B1, B2, rho1, rho2, n1, bXZ, bw){
  library(fields)
  # AR1 structures
  ar1_e <- ar1_cor(E-1, rho1)
  ar1_r <- ar1_cor(R, rho2)
  # Z
  Z <- matrix(0, nrow = n, ncol = nz)
  s     <- expand.grid(1:n1,1:n1)
  W     <- exp(-(rdist(s)/bw)^2)
  W     <- sweep(W,1,rowSums(W),"/")
  for (i in 1:nz){
    Z[,i] <- mvrnorm(1, mu = rep(0, n), Sigma = makecov(sig2Z[i], lambda0, Q))
  }
  smoothZ <- bXZ*W%*%Z
  # M
  M1 <- matrix(0, ncol = E-1, nrow = n)
  for (i in 1:(E-1)){
    M1[,i] <- mvrnorm(1, mu = rep(0, n), Sigma = makecov(sig2X[i], lambda1, Q))
  }
  X <- matrix(0, ncol = E, nrow = n)
  X <- cbind(rep(1, n), smoothZ%*%B1 + M1%*%t(chol(ar1_e)))
  # theta star nxR
  M2 <- matrix(0, ncol = R, nrow = n)
  for (i in 1:R){
    M2[,i] <- mvrnorm(1, mu = rep(0, n), Sigma = makecov(sig2theta[i], lambda2, Q))
  }
  theta <- smoothZ%*%B2 + M2%*%t(chol(ar1_r)) 
  # estar nxR
  e <- matrix(0, nrow = n, ncol = R)
  for (i in 1:R){
    e[,i] <- rnorm(n, 0, sqrt(tau2r[i]))
  }
  Y <- X%*%betas + theta + e
  out <- list(Y = Y, X = X, theta = theta, Z = smoothZ)
  return(out)
}

#########################################################
############### MCMC Functions ##########################
#########################################################

# full model with inverse gamma prior for beta
# lambda1 \sim invgamma(1/2, 1/gamma1) # this is lambdaKHC
# lambda2 \sim invgamma(1/2, 1/gamma2) # this is lambdaEHC
# lambda3 (lambdaRHC) is the same as tauHC
# allow for covariates (Z) and treatments (X)
MCMC_model4_8_3_XZ <- function(Y, X, Z, Qmat, E, n, R, Q, L, K, B, iters = 1000, burn = 50, W, thin = 2){
  tik <- proc.time()
  # Aindex keeps track of the number of entries in each row of A to be updated
  Aindex <- rep(0, R)
  for (r in 1:R){
    if (r <= Q) Aindex[r] <- r-1
    else Aindex[r] <- Q
  }
  iters_thin <- (iters-burn)%/%thin
  # keepers_beta <- array(0, c(n, E, R, iters_thin))
  # the draws group by r (beta, A, tau2) and by q (U and sig2)
  # draws are stored in 3d arrays of the dimension rows by cols by iters
  # keepers_beta <- array(0, c(n, E, R, iters))
  # keepers_gamma <- array(0, c(L, E, R, iters))
  p <- ncol(Z)
  # keepers_alpha <- array(0, c(p, R, iters))
  # keepers_AU <- array(0, c(n,R,iters))
  Amat <- matrix(0, nrow = R, ncol = Q)
  for (r in 1:R){
    if (r <= Q) Amat[r,r] <- 1    
  }
  att <- acc <- 0
  # priors and initial values
  MH <- 0.1
  # W<- eigen(Qmat)$values
  a_tau <- .5
  b_tau <- .005
  Umat <- matrix(0, nrow = n, ncol = Q)
  keepers_Tl <- array(0, c(L, K, iters_thin))
  keepers_Te <- array(0, c(E, K, iters_thin))
  keepers_Tr <- array(0, c(R, K, iters_thin))
  Tl <- matrix(0, nrow = L, ncol = K)
  Te <- matrix(0, nrow = E, ncol = K)
  Tr <- matrix(0, nrow = R, ncol = K)
  sig2 <- rep(0.5, Q)
  tau2 <- rep(0.5, R)
  lambda <- 0.5
  lambdaKHC <- 1
  lambdaEHC <- 1
  tauHC <- 1
  tauQHC <- 1
  nuKHC <- nuEHC <- xiHC <- 1
  X_tensor <- array(X, c(n, E, R))
  estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
  alpha <- matrix(0, nrow = p, ncol = R)
  for (iter in 2:iters){
    print(paste0("This is the ", iter, "th iterations."))
    # update tensor margins for beta
    AU <- Umat%*%t(Amat)
    Zalpha <- Z%*%alpha
    for (k in 1:K){
      for (l in 1:L){
        Cl <- B[,l] * (X %*% outer(Te[,k], Tr[,k]))
        Vl <- sum(colSums(Cl^2)/tau2) + 1/lambdaKHC
        ### Vl <- sum(colSums(Cl^2)/tau2) + 1/10
        Tl_prev <- Tl[l,k]
        Ml <- sum(colSums(Cl*(Y - AU - Zalpha - estimate + Tl_prev*Cl))/tau2)
        Tl[l,k] <- rnorm(1, Ml/Vl, 1/sqrt(Vl))
        # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
        estimate <- estimate + (Tl[l,k] - Tl_prev)*Cl 
      }
      for (e in 1:E){
        Ce <- X[,e] * (B %*% outer(Tl[,k], Tr[,k]))
        Ve <- sum(colSums(Ce^2)/tau2) + 1/lambdaEHC
        ### Ve <- sum(colSums(Ce^2)/tau2) + 1/10
        Te_prev <- Te[e,k]
        Me <- sum(colSums(Ce*(Y - AU - Zalpha - estimate + Te_prev*Ce))/tau2)
        Te[e,k] <- rnorm(1, Me/Ve, 1/sqrt(Ve))
        # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
        estimate <- estimate + (Te[e,k] - Te_prev)*Ce
      }
      for (r in 1:R){
        Cr <- rowSums( X * (B %*% outer(Tl[,k], Te[,k])))
        Vr <- sum(Cr^2)/tau2[r] + 1/(tauHC*tau2[r])
        ### Vr <- sum(Cr^2)/tau2[r] + 1/10
        Tr_prev <- Tr[r,k]
        Mr <- sum(Cr*(Y[,r] - AU[,r] - Zalpha[,r] - estimate[,r] + Tr_prev*Cr)/tau2[r])
        Tr[r,k] <- rnorm(1, Mr/Vr, 1/sqrt(Vr))
        # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
        estimate[,r] <- estimate[,r] + (Tr[r,k] - Tr_prev)*Cr
      }
    }
    # update shrinkage prior parameters
    sumbyl <- rowSums(Tl^2/2)
    sumbye <- rowSums(Te^2/2)
    sumbyr <- rowSums(Tr^2/2)
    
    lambdaKHC <- rinvgamma(1, (L*K+1)/2, sum(sumbyl) + 1/nuKHC)
    lambdaEHC <- rinvgamma(1, (E*K+1)/2, sum(sumbye) + 1/nuEHC)
    tauHC <- rinvgamma(1, (K*R+1)/2, sum(sumbyr/tau2) + 1/xiHC)
    
    nuKHC <- rinvgamma(1, 1, 1/lambdaKHC + 1)
    nuEHC <- rinvgamma(1, 1, 1/lambdaEHC + 1)
    xiHC <- rinvgamma(1, 1, 1/tauHC + 1)
    
    # tauKHC <- rinvgamma(1, (K+1)/2, sum(1/lambdaKHC) + 1)
    # tauEHC <- rinvgamma(1, (E+1)/2, sum(1/lambdaEHC) + 1)
    for (i in 1:R){
      # update lambdaRHC[r] and nuRHC[r]
      # lambdaRHC[i] <- rinvgamma(1, (K+1)/2, sumbyr[i]/(tauHC*tau2[i]) + 1/nuRHC[i])
      ### lambdaRHC[i] <- rinvgamma(1, K/2+ivg_prior_a, sumbyr[i]/(tauHC*tau2[i]) + ivg_prior_b)
      # nuRHC[i] <- rinvgamma(1, 1, 1/lambdaRHC[i] + 1)
      # update tau2_r
      tau2[i] <- rinvgamma(1, (n+K)/2 + a_tau, 0.5*sum((Y[,i]-Zalpha[,i]-estimate[,i]-Umat%*%Amat[i,])^2) + sumbyr[i]/tauHC + b_tau)
      ### tau2[i] <- rinvgamma(1, (n+K)/2 + a_tau, 0.5*sum((Y[,i]-Zalpha[,i]-estimate[,i]-Umat%*%Amat[i,])^2) + b_tau)
      # update A_r
      z <- Aindex[i]
      if (z == 0) {
        next
      }
      # tUK <- sweep(t(Umat[,1:z]), 2, 1/tau2[i], "*")
      # AA <- tUK %*% Umat[,1:z] + diag(1/.25, z)
      tUK <- sweep(t(Umat), 2, 1/tau2[i], "*")
      AA <- tUK %*% Umat + diag(1/.25, Q)
      # solAA <- solve(AA)
      # solAA <- chol2inv(chol(AA))
      BB <- tUK %*% (Y[,i] - Zalpha[,i] - estimate[,i])
      # Amat[i,1:z] <- mvrnorm(1, solAA%*%BB, Sigma = solAA)[1:z]
      Amat[i,1:z] <- spam::rmvnorm.canonical(1, BB, AA)[1:z]
    }
    for (i in 1:R){
      # update alpha_r
      tZK <- sweep(t(Z), 2, 1/tau2[i], "*")
      AAZ <- tZK %*% Z + diag(1/100, p)
      BBZ <- tZK %*% (Y[,i] - AU[,i] - estimate[,i])
      alpha[,i] <- spam::rmvnorm.canonical(1, BBZ, AAZ)
    }
    for (j in 1:Q){
      # update sig2Q
      Kinv_q <- makeKinv(Umat[,j], lambda, W)
      # sig2[j] <- rinvgamma(1, n/2+a_sig2, Kinv_q/2+b_sig2)
      sig2[j] <- rinvgamma(1, (n+1)/2, Kinv_q/2 + 1/tauQHC)
      # update Uq
      Qindex <- setdiff(1:Q, j)
      sumR <- rep(0, n)
      sumArqdivTaur <- 0
      for (w in 1:R){
        sumP <- rep(0, n)
        for (qprime in Qindex){
          sumP <- sumP + Umat[,qprime] * Amat[w,qprime]
        }
        sumArqdivTaur <- sumArqdivTaur + ((Amat[w,j])^2/tau2[w])
        sumR <- sumR + ((Y[,w]-Zalpha[,w]-estimate[,w]-sumP)*Amat[w,j]/tau2[w])
      }
      vark <- sumArqdivTaur+(1-lambda+lambda*W)/sig2[j]
      for (d in 1:n){
        Umat[d,j] <- rnorm(1, sumR[d]/vark[d], 1/sqrt(vark[d]))
      }
    }
    tauQHC <- rinvgamma(1, (Q+1)/2, sum(1/sig2) + 1)
    # update lambdaU
    att <- att + 1
    can <- pnorm(rnorm(1,qnorm(lambda), MH))
    # calculate current and candidate loglikelihood
    curlp <- 0
    canlp <- 0
    for (m in 1:Q){
      curK <- makeK(sig2[m], lambda, W)
      canK <- makeK(sig2[m], can, W)
      curKinv <- makeKinv(Umat[,m], lambda, W)
      canKinv <- makeKinv(Umat[,m], can, W)
      curlp <- curlp + log_like(curK, curKinv, sig2[m])
      canlp <- canlp + log_like(canK, canKinv, sig2[m])
    }
    Rval <- canlp - curlp + dnorm(qnorm(can), log = TRUE) - dnorm(qnorm(lambda), log = TRUE)
    if (!is.na(Rval) & log(runif(1)) < Rval & lambda < .9999){
      acc <- acc + 1
      lambda <- can
    }
    if(iter < burn){
      if(att > 50){
        if(acc/att < 0.3){MH <- MH*0.8}
        if(acc/att > 0.5){MH <- MH*1.2}
        acc <- att <- 0
      }
    }
    
    if (iter > burn){
      if ((iter-burn)%%thin == 0){
        # keepers_beta[,,,(iter-burn)%/%thin] <- bbygamma(B, find_tol_gamma(Tl, Te, Tr), n)  
        keepers_Tl[,,(iter-burn)%/%thin] <- Tl
        keepers_Te[,,(iter-burn)%/%thin] <- Te
        keepers_Tr[,,(iter-burn)%/%thin] <- Tr
      }
    }
    
    # keepers_gamma[,,,iter] <- find_tol_gamma(Tl, Te, Tr)
    # keepers_beta[,,,iter] <- bbygamma(B, find_tol_gamma(Tl, Te, Tr), n)
    # keepers_beta[,,,iter] <- bbygamma(B, find_tol_gamma(Tl, Te, Tr), n)
    # keepers_alpha[,,iter] <- alpha
    # keepers_AU[,,iter] <- AU
  }
  
  tok <- proc.time()
  out <- list(Tl = keepers_Tl, Te = keepers_Te, Tr = keepers_Tr, acc_rate = acc/att, time = tok - tik, MH = MH)
  return(out)
}

# p = total number of alphas and betas
# n = number of sites
# r = number of response variables 
# q = number of factors 
# constant beta across frequency
# added covariates Z
MCMC_model3_2_XZ <- function(Y, X, Z, Qmat, p = 5, n, R, Q, iters = 1000, burn = 50, num1, num2, W, thin=2){
  if (num1 != 1 | num2 != n){
    Y <- Y[num1:num2,]
    X <- X[num1:num2,]
    Z <- Z[num1:num2,]
    n <- num2-num1+1
  }
  tik <- proc.time()
  # Aindex keeps track of the number of entries in each row of A to be updated
  Aindex <- rep(0, R)
  for (r in 1:R){
    if (r <= Q) Aindex[r] <- r-1
    else Aindex[r] <- Q
  }
  # the draws group by r (beta, A, tau2) and by q (U and sig2)
  # draws are stored in 3d arrays of the dimension rows by cols by iters
  iters_thin <- (iters-burn)%/%thin
  keepers_beta <- array(0, c(p, R, iters_thin))
  p2 <- ncol(Z)
  # keepers_alpha <- array(0, c(p2, R, iters))
  Amat <- matrix(0, nrow = R, ncol = Q)
  for (r in 1:R){
    if (r <= Q) Amat[r,r] <- 1    
  }
  att <- acc <- 0
  # priors and initial values
  MH <- 0.1
  # W <- eigen(Qmat)$values[num1:num2]
  cov_B_inv <- diag(1/100, p)
  a_sig2 <- b_sig2 <- .1
  a_tau <- .5
  b_tau <- .005
  # logpdf1 <- 0
  # logpdf2 <- 0
  Umat <- matrix(0, nrow = n, ncol = Q)
  betaz <- keepers_beta[,,1]
  # alpha <- keepers_alpha[,,1]
  alpha <- matrix(0, nrow = p2, ncol = R)
  sig2 <- rep(0.5, Q)
  tau2 <- rep(0.5, R)
  lambda <- 0.5
  
  for (iter in 2:iters){
    print(paste0("This is the ", iter, "th iterations."))
    Zalpha <- Z%*%alpha
    for (i in 1:R){
      # update beta_r
      KtX <- sweep(t(X), 2, 1/tau2[i], "*")
      A <- KtX %*% X + cov_B_inv
      # solA <- solve(A)
      B <- KtX%*%(Y[,i]-Umat%*%Amat[i,]-Zalpha[,i])
      # betaz[,i] <- mvrnorm(1, solA%*%B, solA)
      betaz[,i] <- spam::rmvnorm.canonical(1, B, A)
      # update tau2_r
      tau2[i] <- rinvgamma(1, n/2 + a_tau, 0.5*sum((Y[,i]-X%*%betaz[,i]-Umat%*%Amat[i,]-Zalpha[,i])^2) + b_tau)
      # update A_r
      z <- Aindex[i]
      if (z == 0) {
        next
      }
      # tUK <- sweep(t(Umat[,1:z]), 2, 1/tau2[i], "*")
      # AA <- tUK %*% Umat[,1:z] + diag(1/.25, z)
      tUK <- sweep(t(Umat[,]), 2, 1/tau2[i], "*")
      AA <- tUK %*% Umat[,] + diag(1/.25, Q)
      # solAA <- solve(AA)
      BB <- tUK %*% (Y[,i]-X%*%betaz[,i]-Zalpha[,i])
      # Amat[i,1:z] <- mvrnorm(1, solAA%*%BB, solAA)[1:z]
      Amat[i,1:z] <- spam::rmvnorm.canonical(1, BB, AA)[1:z]
    }
    for (i in 1:R){
      # update alpha_r
      tZK <- sweep(t(Z), 2, 1/tau2[i], "*")
      AAZ <- tZK %*% Z + diag(1/100, p2)
      BBZ <- tZK %*% (Y[,i] - Umat%*%Amat[i,] - X%*%betaz[,i])
      alpha[,i] <- spam::rmvnorm.canonical(1, BBZ, AAZ)
    }
    for (j in 1:Q){
      # update sig2Q
      Kinv_q <- makeKinv(Umat[,j], lambda, W)
      sig2[j] <- rinvgamma(1, n/2+a_sig2, Kinv_q/2+b_sig2)
      # update Uq
      Qindex <- setdiff(1:Q, j)
      sumR <- rep(0, n)
      sumArqdivTaur <- 0
      for (w in 1:R){
        sumP <- rep(0, n)
        for (qprime in Qindex){
          sumP <- sumP + Umat[,qprime] * Amat[w,qprime]
        }
        sumArqdivTaur <- sumArqdivTaur + ((Amat[w,j])^2/tau2[w])
        sumR <- sumR + ((Y[,w]-X%*%betaz[,w]-Zalpha[,w]-sumP)*Amat[w,j]/tau2[w])
      }
      vark <- sumArqdivTaur+(1-lambda+lambda*W)/sig2[j]
      for (k in 1:n){
        Umat[k,j] <- rnorm(1, sumR[k]/vark[k], 1/sqrt(vark[k]))
      }
    }
    # update lambdaU
    att <- att + 1
    can <- pnorm(rnorm(1,qnorm(lambda), MH))
    # calculate current and candidate loglikelihood
    curlp <- 0
    canlp <- 0
    for (m in 1:Q){
      curK <- makeK(sig2[m], lambda, W)
      canK <- makeK(sig2[m], can, W)
      curKinv <- makeKinv(Umat[,m], lambda, W)
      canKinv <- makeKinv(Umat[,m], can, W)
      curlp <- curlp + log_like(curK, curKinv, sig2[m])
      canlp <- canlp + log_like(canK, canKinv, sig2[m])
    }
    Rval <- canlp - curlp + dnorm(qnorm(can), log = TRUE) - dnorm(qnorm(lambda), log = TRUE)
    if (!is.na(Rval) & log(runif(1)) < Rval & lambda < .9999){
      acc <- acc + 1
      lambda <- can
    }
    if(iter < burn){
      if(att > 50){
        if(acc/att < 0.3){MH <- MH*0.8}
        if(acc/att > 0.5){MH <- MH*1.2}
        acc <- att <- 0
      }
    }
    
    if (iter > burn){
      if ((iter-burn)%%thin == 0){
        keepers_beta[,,(iter-burn)%/%thin] <- betaz  
      }
    }
    
    # keepers_beta[,,iter] <- betaz
    # keepers_alpha[,,iter] <- alpha
    # calculate posterior of the likelihood 
    # yhat <- X %*% betaz + Umat%*%t(Amat) + Z%*%alpha
    
    # curll <- matrix(0, nrow = n, ncol = R)
    # for (i in 1:R){
    #   curll[,i] <- -.5*log(tau2[i]) - .5*((Y[,i]-yhat[,i])^2/tau2[i])
    # }
    # 
    # if(iter > burn){
    #   logpdf1 <- logpdf1 + (curll)/(iters-burn)
    #   logpdf2 <- logpdf2 + (curll^2)/(iters-burn)
    # }
  }
  # WAIC computations
  # mn_logpdf  <- logpdf1
  # var_logpdf <- logpdf2 - logpdf1^2
  # pW         <- sum(var_logpdf)
  # WAIC       <- list(WAIC=-2*sum(mn_logpdf)+2*pW,pW=pW)
  
  tok <- proc.time()
  out <- list(beta = keepers_beta, acc_rate = acc/att, time = tok - tik, MH = MH)
  return(out)
}

