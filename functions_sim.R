# five models, shrinkage priors for univariate
# testing out scenarios that make M3 not work well 
# M3 does not allow coefficients to vary by frequency

library(MCMCpack)
library(MASS)
library(splines2)
library(plyr)
library(LaplacesDemon)
library(doParallel)
library(fields)
library(spam)
# library(mvtnorm)
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
# makecov <- function(sigma2, lambda, Q){
#   I <- diag(ncol(Q))
#   cov <- sigma2*solve((1-lambda)*I+lambda*Q)
#   return(cov)
# }

# make covariance matrix for CAR(sigma2, lambda)
makecov <- function(sigma2, lambda, Q){
  I <- diag(ncol(Q))
  cov <- solve((1-lambda)*I+lambda*Q)
  cov <- cov/mean(diag(cov))
  cov <- sigma2*cov
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

# calculate multivariate normal log-likelihood
log_like_theta <- function(theta, Sigma, Sigma_inv, lambda, w){
  -0.5*(log(det(Sigma)))-0.5*(1-lambda+lambda*w)*t(theta)%*%Sigma_inv%*%theta
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

# full MV method with shrinkage prior
MCMC_model4_8_v7 <- function(Y, X, Qmat, E, n, R, Q, L, K, B, iters = 1000, burn = 50, full_result = FALSE, W){
  tik <- proc.time()
  # Aindex keeps track of the number of entries in each row of A to be updated
  Aindex <- rep(0, R)
  for (r in 1:R){
    if (r <= Q) Aindex[r] <- r-1
    else Aindex[r] <- Q
  }
  # the draws group by r (beta, A, tau2) and by q (U and sig2)
  # draws are stored in 3d arrays of the dimension rows by cols by iters
  keepers_yhat <- array(0, c(n, R, iters))
  keepers_beta <- array(0, c(n, E, R, iters))
  keepers_AU <- array(0, c(n, R, iters))
  keepers_Tl <- array(0, c(L, K, iters))
  keepers_Te <- array(0, c(E, K, iters))
  keepers_Tr <- array(0, c(R, K, iters))
  # lambda_j^2 for horseshoe prior
  keepers_lambdaKHC <- matrix(1, nrow = iters, ncol = K)
  keepers_lambdaEHC <- matrix(1, nrow = iters, ncol = E)
  keepers_lambdaRHC <- matrix(1, nrow = iters, ncol = R)
  # nu for horseshoe prior
  keepers_nuRHC <- matrix(1, nrow = iters, ncol = R)
  # tau^2 for horseshoe prior
  keepers_tauHC <- rep(1, iters)
  keepers_tauKHC <- rep(1, iters)
  keepers_tauEHC <- rep(1, iters)
  # keepers_tauQHC <- rep(1, iters)
  # xi for horseshoe prior
  keepers_xiHC <- rep(1, iters)
  keepers_A <- array(0, c(R, Q, iters))
  for (r in 1:R){
    if (r <= Q) keepers_A[r,r,] <- 1    
  }
  keepers_U <- array(0, c(n, Q, iters))
  # keepers_varz contains sigma2Q, tau2R, lambdaU
  keepers_varz <- matrix(0, nrow = iters, ncol = Q+R+1)
  colnames(keepers_varz) <- c(paste0("sig2_q", 1:Q), paste0("tau2_r", 1:R), "lambda")
  keepers_varz[1,] <- 0.5
  att <- acc <- 0
  # priors and initial values
  MH <- 0.1
  # W<- eigen(Qmat)$values
  # a_tau <- b_tau <- .1
  a_tau <- .5
  b_tau <- .005
  # logpdf1 <- 0
  # logpdf2 <- 0
  Amat <- keepers_A[,,1]
  Umat <- keepers_U[,,1]
  Tl <- keepers_Tl[,,1]
  Te <- keepers_Te[,,1]
  Tr <- keepers_Tr[,,1]
  sig2 <- keepers_varz[1, 1:Q]
  tau2 <- keepers_varz[1, (Q+1):(Q+R)]
  lambda <- keepers_varz[1, Q+R+1]
  lambdaKHC <- keepers_lambdaKHC[1,]
  lambdaEHC <- keepers_lambdaEHC[1,]
  lambdaRHC <- keepers_lambdaRHC[1,]
  nuRHC <- keepers_nuRHC[1,]
  tauKHC <- tauEHC <- tauHC <- 1
  xiHC <- 1
  X_tensor <- array(X, c(n, E, R))
  estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
  for (iter in 2:iters){
    print(paste0("This is the ", iter, "th iterations."))
    # update tensor margins for beta
    AU <- Umat%*%t(Amat)
    for (k in 1:K){
      for (l in 1:L){
        # Cl <- xbybeta3(X_tensor, outer(B[,l],outer(Te[,k],Tr[,k])), R)
        Cl <- B[,l] * (X %*% outer(Te[,k], Tr[,k]))
        Vl <- sum(colSums(Cl^2)/tau2) + 1/(lambdaKHC[k])
        # Vl <- sum(colSums(Cl^2)/tau2) + 1
        Tl_prev <- Tl[l,k]
        Ml <- sum(colSums(Cl*(Y - AU - estimate + Tl_prev*Cl))/tau2)
        Tl[l,k] <- rnorm(1, Ml/Vl, 1/sqrt(Vl))
        # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
        estimate <- estimate + (Tl[l,k] - Tl_prev)*Cl
      }
      for (e in 1:E){
        # Ce <- xbybeta3(X[,e], array((B %*% outer(Tl[,k], Tr[,k])), c(n,1,R)), R)
        Ce <- X[,e] * (B %*% outer(Tl[,k], Tr[,k]))
        Ve <- sum(colSums(Ce^2)/tau2) + 1/(lambdaEHC[e])
        # Ve <- sum(colSums(Ce^2)/tau2) + 1
        Te_prev <- Te[e,k]
        Me <- sum(colSums(Ce*(Y - AU - estimate + Te_prev*Ce))/tau2)
        Te[e,k] <- rnorm(1, Me/Ve, 1/sqrt(Ve))
        # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
        estimate <- estimate + (Te[e,k] - Te_prev)*Ce
      }
      for (r in 1:R){
        Cr <- rowSums( X * (B %*% outer(Tl[,k], Te[,k])))
        Vr <- sum(Cr^2)/tau2[r] + 1/(lambdaRHC[r]*tauHC*tau2[r])
        # Vr <- sum(Cr^2)/tau2[r] + 1
        Tr_prev <- Tr[r,k]
        Mr <- sum(Cr*(Y[,r] - AU[,r] - estimate[,r] + Tr_prev*Cr)/tau2[r])
        Tr[r,k] <- rnorm(1, Mr/Vr, 1/sqrt(Vr))
        # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
        estimate[,r] <- estimate[,r] + (Tr[r,k] - Tr_prev)*Cr
      }
      # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
      # update lambdaKHC[k]
      # sum_tauHC[k] <- sum(Tr[,k]^2/(2*lambdaRHC*tau2))
      lambdaKHC[k] <- rinvgamma(1, (L+1)/2, sum(Tl[,k]^2)/2 + 1/tauKHC)
    }
    # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
    for (e in 1:E){
      # update lambdaEHC[e] and nuEHC[e]
      lambdaEHC[e] <- rinvgamma(1, (K+1)/2, sum(Te[e,]^2)/2 + 1/tauEHC)
    }
    # update shrinkage prior parameters
    sumbyr <- rowSums(Tr^2/2)
    tauHC <- rinvgamma(1, (K*R+1)/2, sum(sumbyr/(lambdaRHC*tau2)) + 1/xiHC)
    xiHC <- rinvgamma(1, 1, 1/tauHC + 1)
    tauKHC <- rinvgamma(1, (K+1)/2, sum(1/lambdaKHC) + 1)
    tauEHC <- rinvgamma(1, (E+1)/2, sum(1/lambdaEHC) + 1)
    for (i in 1:R){
      # update lambdaRHC[r] and nuRHC[r]
      lambdaRHC[i] <- rinvgamma(1, (K+1)/2, sumbyr[i]/(tauHC*tau2[i]) + 1/nuRHC[i])
      nuRHC[i] <- rinvgamma(1, 1, 1/lambdaRHC[i] + 1)
      # update tau2_r
      tau2[i] <- rinvgamma(1, (n+K)/2 + a_tau, 0.5*sum((Y[,i]-estimate[,i]-Umat%*%Amat[i,])^2) + sumbyr[i]/(lambdaRHC[i]*tauHC) + b_tau)
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
      BB <- tUK %*% (Y[,i] - estimate[,i])
      # Amat[i,1:z] <- mvrnorm(1, solAA%*%BB, Sigma = solAA)[1:z]
      Amat[i,1:z] <- spam::rmvnorm.canonical(1, BB, AA)[1:z]
    }
    for (j in 1:Q){
      # update sig2Q
      Kinv_q <- makeKinv(Umat[,j], lambda, W)
      # sig2[j] <- rinvgamma(1, n/2 + a_sig2, Kinv_q/2 + b_sig2)
      sig2[j] <- rinvgamma(1, n/2 + .5, Kinv_q/2 + .005)
      # sig2[j] <- rinvgamma(1, (n+1)/2, Kinv_q/2 + 1/tauQHC)
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
        sumR <- sumR + ((Y[,w]-estimate[,w]-sumP)*Amat[w,j]/tau2[w])
      }
      vark <- sumArqdivTaur+(1-lambda+lambda*W)/sig2[j]
      for (d in 1:n){
        Umat[d,j] <- rnorm(1, sumR[d]/vark[d], 1/sqrt(vark[d]))
      }
      # Umat <- theta_mat
    }
    # tauQHC <- rinvgamma(1, (Q+1)/2, sum(1/sig2) + 1)
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
    if (!is.na(Rval) & log(runif(1)) < Rval & lambda <.9999){
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
    # lambda <- .95
    AU <- Umat%*%t(Amat)
    yhat <- estimate + AU
    
    # calculate posterior of the likelihood 
    # curll <- matrix(0, nrow = n, ncol = R)
    # for (r in 1:R){
    #   curll[,r] <- -.5*log(tau2[r]) - .5*((Y[,r]-yhat[,r])^2/tau2[r])
    # }
    # 
    # if(iter > burn){
    #   logpdf1 <- logpdf1 + (curll)/(iters-burn)
    #   logpdf2 <- logpdf2 + (curll^2)/(iters-burn)
    # }
    
    keepers_Tl[,,iter] <- Tl
    keepers_Te[,,iter] <- Te
    keepers_Tr[,,iter] <- Tr
    keepers_lambdaKHC[iter,] <- lambdaKHC
    keepers_lambdaEHC[iter,] <- lambdaEHC
    keepers_lambdaRHC[iter,] <- lambdaRHC
    keepers_nuRHC[iter,] <- nuRHC
    keepers_tauHC[iter] <- tauHC
    keepers_tauKHC[iter] <- tauKHC
    keepers_tauEHC[iter] <- tauEHC
    # keepers_tauQHC[iter] <- tauQHC
    keepers_xiHC[iter] <- xiHC
    keepers_A[,,iter] <- Amat
    keepers_U[,,iter] <- Umat
    keepers_varz[iter,] <- c(sig2, tau2, lambda)
    keepers_yhat[,,iter] <- yhat
    keepers_beta[,,,iter] <- bbygamma(B, find_tol_gamma(Tl, Te, Tr), n)
    keepers_AU[,,iter] <- AU
  }
  # WAIC computations
  # mn_logpdf  <- logpdf1
  # var_logpdf <- logpdf2 - logpdf1^2
  # pW         <- sum(var_logpdf)
  # WAIC       <- list(WAIC=-2*sum(mn_logpdf)+2*pW,pW=pW)
  
  tok <- proc.time()
  # out <- list(Tl = keepers_Tl, Te = keepers_Te, Tr = keepers_Tr, lambdaKHC = keepers_lambdaKHC, lambdaEHC = keepers_lambdaEHC, lambdaRHC = keepers_lambdaRHC, tauKHC = keepers_tauKHC, tauEHC = keepers_tauEHC, nuRHC = keepers_nuRHC, tauHC = keepers_tauHC, xiHC = keepers_xiHC, A = keepers_A, U = keepers_U, yhat = keepers_yhat, beta = keepers_beta, AU = keepers_AU, varz = keepers_varz, Aindex = Aindex, acc_rate = acc/att, time = tok - tik, MH = MH)
  if (full_result){
    out <- list(beta = keepers_beta[,,,(burn+1):iters], varz = keepers_varz, AU = keepers_AU, acc_rate = acc/att, time = tok - tik, MH = MH)
  }
  else{
    out <- list(beta = keepers_beta[1,,,(burn+1):iters], varz = keepers_varz, AU = keepers_AU, acc_rate = acc/att, time = tok - tik, MH = MH)
  }
  return(out)
}

# univariate with shrinkage prior
MCMC_model4_9_v2 <- function(Y, X, Qmat, E, n, L, K, B, iters = 1000, burn = 50){
  tik <- proc.time()
  # the draws group by r (beta, A, tau2) and by q (U and sig2)
  # draws are stored in 3d arrays of the dimension rows by cols by iters
  # keepers_yhat <- array(0, c(n, R, iters))
  R <- Q <- 1
  keepers_beta <- array(0, c(n, E, iters))
  keepers_Tl <- array(0, c(L, K, iters))
  keepers_Te <- array(0, c(E, K, iters))
  keepers_tauKHC <- rep(1, iters)
  # keepers_tauEHC <- rep(1, iters)
  # lambda_j^2 for horseshoe prior
  keepers_lambdaKHC <- matrix(1, nrow = iters, ncol = K)
  keepers_lambdaEHC <- matrix(1, nrow = iters, ncol = E)
  keepers_tauHC <- rep(1, iters)
  # nu for horseshoe prior
  keepers_nuEHC <- matrix(1, nrow = iters, ncol = E)
  keepers_xiHC <- rep(1, iters)
  keepers_U <- array(0, c(n, iters))
  # keepers_varz contains sigma2Q, tau2R, lambdaU
  keepers_varz <- matrix(0, nrow = iters, ncol = Q+R+1)
  keepers_varz[1,] <- 0.5
  att <- acc <- 0
  # priors and initial values
  MH <- 0.1
  W<- eigen(Qmat)$values
  a_tau <- .5
  b_tau <- .005
  Umat <- keepers_U[,1]
  Tl <- keepers_Tl[,,1]
  Te <- keepers_Te[,,1]
  sig2 <- keepers_varz[1, 1:Q]
  tau2 <- keepers_varz[1, (Q+1):(Q+R)]
  lambda <- keepers_varz[1, Q+R+1]
  lambdaKHC <- keepers_lambdaKHC[1,]
  lambdaEHC <- keepers_lambdaEHC[1,]
  nuEHC <- keepers_nuEHC[1,]
  tauKHC <- 1
  tauHC <- keepers_tauHC[1]
  xiHC <- 1
  # X_tensor <- array(X, c(n, E, R))
  estimate <- find_estimate3(X, B, Tl, Te)
  for (iter in 2:iters){
    # update tensor margins for beta
    for (k in 1:K){
      for (l in 1:L){
        Cl <- B[,l] * (X %*% Te[,k])
        Vl <- sum(colSums(Cl^2)/tau2) + 1/(lambdaKHC[k])
        Tl_prev <- Tl[l,k]
        Ml <- sum(colSums(Cl*(Y - Umat - estimate + Tl_prev*Cl))/tau2)
        Tl[l,k] <- rnorm(1, Ml/Vl, 1/sqrt(Vl))
        estimate <- estimate + (Tl[l,k] - Tl_prev)*Cl 
      }
      for (e in 1:E){
        Ce <- X[,e] * (B %*% Tl[,k])
        # Ve <- sum(colSums(Ce^2)/tau2) + 1/(lambdaEHC[e]*tauHC*tau2)
        Ve <- sum(colSums(Ce^2)/tau2) + 1/(lambdaEHC[e]*tauHC)
        Te_prev <- Te[e,k]
        Me <- sum(colSums(Ce*(Y - Umat - estimate + Te_prev*Ce))/tau2)
        Te[e,k] <- rnorm(1, Me/Ve, 1/sqrt(Ve))
        estimate <- estimate + (Te[e,k] - Te_prev)*Ce
      }
      # estimate <- find_estimate3(X, B, Tl, Te)
      # update lambdaKHC[k]
      lambdaKHC[k] <- rinvgamma(1, (L+1)/2, sum(Tl[,k]^2)/2 + 1/tauKHC)
    }
    # update shrinkage prior parameters
    sumbye <- rowSums(Te^2/2)
    # tauHC <- rinvgamma(1, (K*E+1)/2, sum(sumbye/(lambdaEHC*tau2)) + 1/xiHC)
    tauHC <- rinvgamma(1, (K*E+1)/2, sum(sumbye/(lambdaEHC)) + 1/xiHC)
    xiHC <- rinvgamma(1, 1, 1/tauHC + 1)
    tauKHC <- rinvgamma(1, (K+1)/2, sum(1/lambdaKHC) + 1)
    
    for (e in 1:E){
      # update lambdaEHC[e] and nuEHC[e]
      # lambdaEHC[e] <- rinvgamma(1, (E+1)/2, sumbye[e]/(tauHC*tau2) + 1/nuEHC[e])
      lambdaEHC[e] <- rinvgamma(1, (K+1)/2, sumbye[e]/(tauHC) + 1/nuEHC[e])
      nuEHC[e] <- rinvgamma(1, 1, 1/lambdaEHC[e] + 1)
    }
    # update tau2_r
    # tau2 <- rinvgamma(1, (n+K*E)/2 + a_tau, 0.5*sum((Y-estimate-Umat)^2) + sum(sumbye/(lambdaEHC*tauHC)) + b_tau)
    tau2 <- rinvgamma(1, n/2 + a_tau, 0.5*sum((Y-estimate-Umat)^2) + b_tau)
    # update sig2Q
    Kinv_q <- makeKinv(Umat, lambda, W)
    sig2 <- rinvgamma(1, n/2 + .5, Kinv_q/2 + .005)
    # sig2 <- rinvgamma(1, (n+1)/2, Kinv_q/2 + 1/tauQHC)
    # update Uq
    # Qindex <- setdiff(1:Q, j)
    # sumR <- rep(0, n)
    # sumArqdivTaur <- 0
    # sumR <- sumR + ((Y-estimate)/tau2)
    sumR <- (Y-estimate)/tau2
    # for (w in 1:R){
    # sumP <- rep(0, n)
    # for (qprime in Qindex){
    #   sumP <- sumP + Umat[,qprime] * Amat[w,qprime]
    # }
    # sumArqdivTaur <- sumArqdivTaur + ((Amat[w,j])^2/tau2[w])
    # sumR <- sumR + ((Y-estimate[,w]-sumP)*Amat[w]/tau2[w])
    # }
    # vark <- sumArqdivTaur+(1-lambda+lambda*W)/sig2[j]
    vark <- 1/tau2 +(1-lambda+lambda*W)/sig2
    for (d in 1:n){
      Umat[d] <- rnorm(1, sumR[d]/vark[d], 1/sqrt(vark[d]))
    }
    # tauQHC <- rinvgamma(1, (Q+1)/2, sum(1/sig2) + 1)
    # update lambdaU
    att <- att + 1
    can <- pnorm(rnorm(1,qnorm(lambda), MH))
    # calculate current and candidate loglikelihood
    curlp <- 0
    canlp <- 0
    curK <- makeK(sig2, lambda, W)
    canK <- makeK(sig2, can, W)
    curKinv <- makeKinv(Umat, lambda, W)
    canKinv <- makeKinv(Umat, can, W)
    curlp <- curlp + log_like(curK, curKinv, sig2)
    canlp <- canlp + log_like(canK, canKinv, sig2)
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
    
    yhat <- estimate + Umat
    
    keepers_tauKHC[iter] <- tauKHC
    keepers_lambdaKHC[iter,] <- lambdaKHC
    keepers_lambdaEHC[iter,] <- lambdaEHC
    keepers_nuEHC[iter,] <- nuEHC
    keepers_tauHC[iter] <- tauHC
    keepers_xiHC[iter] <- xiHC
    keepers_varz[iter,] <- c(sig2, tau2, lambda)
    keepers_beta[,,iter] <- B %*% find_tol_gamma3(Tl, Te)
    keepers_U[,iter] <- Umat
  }
  
  tok <- proc.time()
  out <- list(beta = keepers_beta, tauKHC = keepers_tauKHC, tauHC = keepers_tauHC, nuEHC = keepers_nuEHC, lambdaKHC = keepers_lambdaKHC, lambdaEHC = keepers_lambdaEHC, xiHC = keepers_xiHC, varz = keepers_varz, U = keepers_U, acc_rate = acc/att, time = tok - tik, MH = MH)
  return(out)
}

# full model with inverse gamma prior for beta
MCMC_model4_8_IVG <- function(Y, X, Qmat, E, n, R, Q, L, K, B, iters = 1000, burn = 50){
  tik <- proc.time()
  # Aindex keeps track of the number of entries in each row of A to be updated
  ivg_prior_a <- .5
  ivg_prior_b <- .005
  Aindex <- rep(0, R)
  for (r in 1:R){
    if (r <= Q) Aindex[r] <- r-1
    else Aindex[r] <- Q
  }
  # the draws group by r (beta, A, tau2) and by q (U and sig2)
  # draws are stored in 3d arrays of the dimension rows by cols by iters
  # keepers_yhat <- array(0, c(n, R, iters))
  keepers_beta <- array(0, c(n, E, R, iters))
  keepers_AU <- array(0, c(n, R, iters))
  keepers_Tl <- array(0, c(L, K, iters))
  keepers_Te <- array(0, c(E, K, iters))
  keepers_Tr <- array(0, c(R, K, iters))
  # lambda_j^2 for horseshoe prior
  keepers_lambdaKHC <- matrix(1, nrow = iters, ncol = K)
  keepers_lambdaEHC <- matrix(1, nrow = iters, ncol = E)
  keepers_lambdaRHC <- matrix(1, nrow = iters, ncol = R)
  # nu for horseshoe prior
  # keepers_nuRHC <- matrix(1, nrow = iters, ncol = R)
  # tau^2 for horseshoe prior
  keepers_tauHC <- rep(1, iters)
  # keepers_tauKHC <- rep(1, iters)
  # keepers_tauEHC <- rep(1, iters)
  keepers_tauQHC <- rep(1, iters)
  # xi for horseshoe prior
  # keepers_xiHC <- rep(1, iters)
  keepers_A <- array(0, c(R, Q, iters))
  for (r in 1:R){
    if (r <= Q) keepers_A[r,r,] <- 1    
  }
  keepers_U <- array(0, c(n, Q, iters))
  # keepers_varz contains sigma2Q, tau2R, lambdaU
  keepers_varz <- matrix(0, nrow = iters, ncol = Q+R+1)
  colnames(keepers_varz) <- c(paste0("sig2_q", 1:Q), paste0("tau2_r", 1:R), "lambda")
  keepers_varz[1,] <- 0.5
  att <- acc <- 0
  # priors and initial values
  MH <- 0.1
  W<- eigen(Qmat)$values
  a_tau <- b_tau <- .1
  # ogpdf1 <- 0
  # logpdf2 <- 0
  Amat <- keepers_A[,,1]
  Umat <- keepers_U[,,1]
  Tl <- keepers_Tl[,,1]
  Te <- keepers_Te[,,1]
  Tr <- keepers_Tr[,,1]
  sig2 <- keepers_varz[1, 1:Q]
  tau2 <- keepers_varz[1, (Q+1):(Q+R)]
  lambda <- keepers_varz[1, Q+R+1]
  lambdaKHC <- keepers_lambdaKHC[1,]
  lambdaEHC <- keepers_lambdaEHC[1,]
  lambdaRHC <- keepers_lambdaRHC[1,]
  # nuRHC <- keepers_nuRHC[1,]
  # tauKHC <- tauEHC <- tauHC <- tauQHC <- 1
  # tauHC <- tauQHC <- 1
  # tauKHC <- keepers_tauKHC[1]
  # tauEHC <- keepers_tauEHC[1]
  tauHC <- keepers_tauHC[1]
  tauQHC <- keepers_tauQHC[1]
  # xiHC <- keepers_xiHC[1]
  # xiHC <- 1
  X_tensor <- array(X, c(n, E, R))
  estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
  for (iter in 2:iters){
    # update tensor margins for beta
    AU <- Umat%*%t(Amat)
    for (k in 1:K){
      for (l in 1:L){
        Cl <- B[,l] * (X %*% outer(Te[,k], Tr[,k]))
        Vl <- sum(colSums(Cl^2)/tau2) + 1/(lambdaKHC[k])
        # Vl <- sum(colSums(Cl^2)/tau2) + 1
        Tl_prev <- Tl[l,k]
        Ml <- sum(colSums(Cl*(Y - AU - estimate + Tl_prev*Cl))/tau2)
        Tl[l,k] <- rnorm(1, Ml/Vl, 1/sqrt(Vl))
        # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
        estimate <- estimate + (Tl[l,k] - Tl_prev)*Cl 
      }
      for (e in 1:E){
        Ce <- X[,e] * (B %*% outer(Tl[,k], Tr[,k]))
        Ve <- sum(colSums(Ce^2)/tau2) + 1/(lambdaEHC[e])
        # Ve <- sum(colSums(Ce^2)/tau2) + 1
        Te_prev <- Te[e,k]
        Me <- sum(colSums(Ce*(Y - AU - estimate + Te_prev*Ce))/tau2)
        Te[e,k] <- rnorm(1, Me/Ve, 1/sqrt(Ve))
        # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
        estimate <- estimate + (Te[e,k] - Te_prev)*Ce
      }
      for (r in 1:R){
        Cr <- rowSums( X * (B %*% outer(Tl[,k], Te[,k])))
        Vr <- sum(Cr^2)/tau2[r] + 1/(lambdaRHC[r]*tauHC*tau2[r])
        # Vr <- sum(Cr^2)/tau2[r] + 1
        Tr_prev <- Tr[r,k]
        Mr <- sum(Cr*(Y[,r] - AU[,r] - estimate[,r] + Tr_prev*Cr)/tau2[r])
        Tr[r,k] <- rnorm(1, Mr/Vr, 1/sqrt(Vr))
        # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
        estimate[,r] <- estimate[,r] + (Tr[r,k] - Tr_prev)*Cr
      }
      # update lambdaKHC[k]
      # sum_tauHC[k] <- sum(Tr[,k]^2/(2*lambdaRHC*tau2))
      # lambdaKHC[k] <- rinvgamma(1, (L+1)/2, sum(Tl[,k]^2)/2 + 1/tauKHC)
      lambdaKHC[k] <- rinvgamma(1, L/2+ivg_prior_a, sum(Tl[,k]^2)/2 + ivg_prior_b)
    }
    for (e in 1:E){
      # update lambdaEHC[e] and nuEHC[e]
      # lambdaEHC[e] <- rinvgamma(1, (K+1)/2, sum(Te[e,]^2)/2 + 1/tauEHC)
      lambdaEHC[e] <- rinvgamma(1, K/2 + ivg_prior_a, sum(Te[e,]^2)/2 + ivg_prior_b)
    }
    # update shrinkage prior parameters
    sumbyr <- rowSums(Tr^2/2)
    # tauHC <- rinvgamma(1, (K*R+1)/2, sum(sumbyr/(lambdaRHC*tau2)) + 1/xiHC)
    tauHC <- rinvgamma(1, K*R/2+ivg_prior_a, sum(sumbyr/(lambdaRHC*tau2)) + ivg_prior_b)
    # xiHC <- rinvgamma(1, 1, 1/tauHC + 1)
    # tauKHC <- rinvgamma(1, (K+1)/2, sum(1/lambdaKHC) + 1)
    # tauEHC <- rinvgamma(1, (E+1)/2, sum(1/lambdaEHC) + 1)
    for (i in 1:R){
      # update lambdaRHC[r] and nuRHC[r]
      # lambdaRHC[i] <- rinvgamma(1, (K+1)/2, sumbyr[i]/(tauHC*tau2[i]) + 1/nuRHC[i])
      lambdaRHC[i] <- rinvgamma(1, K/2+ivg_prior_a, sumbyr[i]/(tauHC*tau2[i]) + ivg_prior_b)
      # nuRHC[i] <- rinvgamma(1, 1, 1/lambdaRHC[i] + 1)
      # update tau2_r
      tau2[i] <- rinvgamma(1, (n+K)/2 + a_tau, 0.5*sum((Y[,i]-estimate[,i]-Umat%*%Amat[i,])^2) + sumbyr[i]/(lambdaRHC[i]*tauHC) + b_tau)
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
      BB <- tUK %*% (Y[,i] - estimate[,i])
      # Amat[i,1:z] <- mvrnorm(1, solAA%*%BB, Sigma = solAA)[1:z]
      Amat[i,1:z] <- spam::rmvnorm.canonical(1, BB, AA)[1:z]
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
        sumR <- sumR + ((Y[,w]-estimate[,w]-sumP)*Amat[w,j]/tau2[w])
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
    if (!is.na(Rval) & log(runif(1)) < Rval){
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
    
    AU <- Umat%*%t(Amat)
    yhat <- estimate + AU
    
    # calculate posterior of the likelihood 
    # curll <- matrix(0, nrow = n, ncol = R)
    # for (r in 1:R){
    #   curll[,r] <- -.5*log(tau2[r]) - .5*((Y[,r]-yhat[,r])^2/tau2[r])
    # }
    # 
    # if(iter > burn){
    #   logpdf1 <- logpdf1 + (curll)/(iters-burn)
    #   logpdf2 <- logpdf2 + (curll^2)/(iters-burn)
    # }
    
    keepers_Tl[,,iter] <- Tl
    keepers_Te[,,iter] <- Te
    keepers_Tr[,,iter] <- Tr
    keepers_lambdaKHC[iter,] <- lambdaKHC
    keepers_lambdaEHC[iter,] <- lambdaEHC
    keepers_lambdaRHC[iter,] <- lambdaRHC
    # keepers_nuRHC[iter,] <- nuRHC
    keepers_tauHC[iter] <- tauHC
    # keepers_tauKHC[iter] <- tauKHC
    # keepers_tauEHC[iter] <- tauEHC
    keepers_tauQHC[iter] <- tauQHC
    # keepers_xiHC[iter] <- xiHC
    keepers_A[,,iter] <- Amat
    keepers_U[,,iter] <- Umat
    keepers_varz[iter,] <- c(sig2, tau2, lambda)
    # keepers_yhat[,,iter] <- yhat
    keepers_beta[,,,iter] <- bbygamma(B, find_tol_gamma(Tl, Te, Tr), n)
    keepers_AU[,,iter] <- AU
  }
  # WAIC computations
  # mn_logpdf  <- logpdf1
  # var_logpdf <- logpdf2 - logpdf1^2
  # pW         <- sum(var_logpdf)
  # WAIC       <- list(WAIC=-2*sum(mn_logpdf)+2*pW,pW=pW)
  # 
  tok <- proc.time()
  # out <- list(beta = keepers_beta[,,,(burn+1):iters], acc_rate = acc/att, WAIC = WAIC, time = tok - tik, MH = MH)
  out <- list(Tl = keepers_Tl, Te = keepers_Te, Tr = keepers_Tr, lambdaKHC = keepers_lambdaKHC, lambdaEHC = keepers_lambdaEHC, lambdaRHC = keepers_lambdaRHC, tauHC = keepers_tauHC, tauQHC = keepers_tauQHC, A = keepers_A, U = keepers_U, beta = keepers_beta, AU = keepers_AU, varz = keepers_varz, Aindex = Aindex, acc_rate = acc/att, time = tok - tik, MH = MH)
  return(out)
}

# p = total number of alphas and betas
# n = number of sites
# r = number of response variables 
# q = number of factors 
# constant beta across frequency
MCMC_model3_2 <- function(Y, X, Q, p = 5, n, r, q, iters = 1000, burn = 50, num1, num2, W = NULL){
  if (num1 != 1 | num2 != n){
    Y <- Y[num1:num2,]
    X <- X[num1:num2,]
    n <- num2-num1+1
  }
  tik <- proc.time()
  # Aindex keeps track of the number of entries in each row of A to be updated
  Aindex <- rep(0, r)
  for (i in 1:r){
    if (i <= q) Aindex[i] <- i-1
    else Aindex[i] <- q
  }
  # the draws group by r (beta, A, tau2) and by q (U and sig2)
  # draws are stored in 3d arrays of the dimension rows by cols by iters
  keepers_beta <- array(0, c(p, r, iters))
  # colnames(keepers_beta) <- c(paste0("r", 1:r))
  # rownames(keepers_beta) <- c(paste0("p", 1:p))
  keepers_A <- array(0, c(r, q, iters))
  # colnames(keepers_A) <- c(paste0("q", 1:q))
  # rownames(keepers_A) <- c(paste0("r", 1:r))
  for (i in 1:r){
    if (i <= q) keepers_A[i,i,] <- 1    
  }
  keepers_U <- array(0, c(n, q, iters))
  # colnames(keepers_U) <- c(paste0("q", 1:q))
  # rownames(keepers_U) <- c(paste0("n", 1:n))
  # keepers_varz contains sigma2Q, tau2R, lambdaU
  keepers_varz <- matrix(0, nrow = iters, ncol = q+r+1)
  # colnames(keepers_varz) <- c(paste0("sig2_q", 1:q), paste0("tau2_r", 1:r), "lambda")
  keepers_varz[1,] <- 0.5
  att <- acc <- 0
  # priors and initial values
  MH <- 0.1
  if (is.null(W)){
    W <- eigen(Q)$values[num1:num2] 
  }
  cov_B_inv <- diag(1/100, p)
  a_tau <- b_tau <- a_sig2 <- b_sig2 <- .1
  # logpdf1 <- 0
  # logpdf2 <- 0
  Amat <- keepers_A[,,1]
  # Amat <- matrix(c(1, 0, 0.5, 1), byrow = T, nrow = 2)
  # Umat <- Umat
  Umat <- keepers_U[,,1]
  betaz <- keepers_beta[,,1]
  sig2 <- keepers_varz[1, 1:q]
  tau2 <- keepers_varz[1, (q+1):(q+r)]
  lambda <- keepers_varz[1, q+r+1]
  
  for (iter in 2:iters){
    print(paste0("This is the ", iter, "th iterations."))
    for (i in 1:r){
      # update beta_r
      KtX <- sweep(t(X), 2, 1/tau2[i], "*")
      A <- KtX %*% X + cov_B_inv
      # solA <- solve(A)
      B <- KtX%*%(Y[,i]-Umat%*%Amat[i,])
      # betaz[,i] <- mvrnorm(1, solA%*%B, solA)
      betaz[,i] <- spam::rmvnorm.canonical(1, B, A)
      # update tau2_r
      tau2[i] <- rinvgamma(1, n/2 + a_tau, 0.5*sum((Y[,i]-X%*%betaz[,i]-Umat%*%Amat[i,])^2) + b_tau)
      # update A_r
      z <- Aindex[i]
      if (z == 0) {
        next
      }
      # tUK <- sweep(t(Umat[,1:z]), 2, 1/tau2[i], "*")
      # AA <- tUK %*% Umat[,1:z] + diag(1/.25, z)
      tUK <- sweep(t(Umat[,]), 2, 1/tau2[i], "*")
      AA <- tUK %*% Umat[,] + diag(1/.25, q)
      # solAA <- solve(AA)
      BB <- tUK %*% (Y[,i]-X%*%betaz[,i])
      # Amat[i,1:z] <- mvrnorm(1, solAA%*%BB, solAA)[1:z]
      Amat[i,1:z] <- spam::rmvnorm.canonical(1, BB, AA)[1:z]
    }
    for (j in 1:q){
      # update sig2Q
      Kinv_q <- makeKinv(Umat[,j], lambda, W)
      sig2[j] <- rinvgamma(1, n/2+a_sig2, Kinv_q/2+b_sig2)
      # update Uq
      Qindex <- setdiff(1:q, j)
      sumR <- rep(0, n)
      sumArqdivTaur <- 0
      for (w in 1:r){
        sumP <- rep(0, n)
        for (qprime in Qindex){
          sumP <- sumP + Umat[,qprime] * Amat[w,qprime]
        }
        sumArqdivTaur <- sumArqdivTaur + ((Amat[w,j])^2/tau2[w])
        sumR <- sumR + ((Y[,w]-X%*%betaz[,w]-sumP)*Amat[w,j]/tau2[w])
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
    for (m in 1:q){
      curK <- makeK(sig2[m], lambda, W)
      canK <- makeK(sig2[m], can, W)
      curKinv <- makeKinv(Umat[,m], lambda, W)
      canKinv <- makeKinv(Umat[,m], can, W)
      curlp <- curlp + log_like(curK, curKinv, sig2[m])
      canlp <- canlp + log_like(canK, canKinv, sig2[m])
    }
    R <- canlp - curlp + dnorm(qnorm(can), log = TRUE) - dnorm(qnorm(lambda), log = TRUE)
    if (!is.na(R) & log(runif(1)) < R){
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
    keepers_beta[,,iter] <- betaz
    # keepers_A[,,iter] <- Amat
    # keepers_U[,,iter] <- Umat
    # keepers_varz[iter,] <- c(sig2, tau2, lambda)
    # calculate posterior of the likelihood 
    # yhat <- X %*% betaz + Umat%*%t(Amat)
    
    # curll <- matrix(0, nrow = n, ncol = r)
    # for (i in 1:r){
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
  out <- list(beta = keepers_beta[,,(burn+1):iters], acc_rate = acc/att, time = tok - tik, MH = MH)
  return(out)
}


###############################################################
############## Simulation Functions ###########################
###############################################################

# 8 competing models
# use uncorrelated B functions
# use extra smoothing for Z
# no standard normal prior
sim_func4_1_tuneL_6 <- function(n, nz, E, R, L, K, Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, betas, B1, B2, rho1, rho2, iters, burn, num1, num2, n1, bXZ, bw, i){
  set.seed(i)
  library(MASS)
  library(plyr)
  library(splines2)
  library(MCMCpack)
  EI <- eigen(Qmat)
  Gamma <- EI$vectors
  W <- EI$values
  B <- bSpline(W,df=L,intercept=T)
  B20 <- bSpline(W,df = 20, intercept = T)
  # B30 <- bSpline(W,df = 30, intercept = T)
  set <- create_set7(n, nz, E, R, Qmat, sig2X, sig2Z, sig2theta, tau2r, lambdaZ, lambda1, lambda2, betas, B1, B2, rho1, rho2, n1, bXZ, bw)
  Ystar <- t(Gamma)%*%set$Y
  Xstar <- t(Gamma)%*%set$X
  # mean, first, median, third, sd, time
  results <- array(0, c(n, E, R, 7, 8))
  test <- list()
  # full model, HC prior
  test[[1]] <- MCMC_model4_8_v7(Ystar, Xstar, Qmat = Qmat, E = E, n = n, R = R, Q = R, L, K, B = B, iters = iters, burn = burn, full_result = T, W=W)
  test[[2]] <- MCMC_model4_8_v7(Ystar, Xstar, Qmat = Qmat, E = E, n = n, R = R, Q = R, L=20, K, B = B20, iters = iters, burn = burn, full_result = T, W=W)
  # full model, IG prior
  test[[3]] <- MCMC_model4_8_IVG(Ystar, Xstar, Qmat = Qmat, E = E, n = n, R = R, Q = R, L, K, B = B, iters = iters, burn = burn)
  test[[4]] <- MCMC_model4_8_IVG(Ystar, Xstar, Qmat = Qmat, E = E, n = n, R = R, Q = R, L=20, K, B=B20, iters = iters, burn = burn)
  for (j in 1:n){
    for (k in 1:E){
      for (l in 1:R){
        for (m in 1:4){
          results[j,k,l,1,m] <- mean(test[[m]]$beta[j,k,l,])
          results[j,k,l,2,m] <- quantile(test[[m]]$beta[j,k,l,], probs = c(0.025))
          results[j,k,l,3,m] <- quantile(test[[m]]$beta[j,k,l,], probs = c(0.5))
          results[j,k,l,4,m] <- quantile(test[[m]]$beta[j,k,l,], probs = c(0.975))
          results[j,k,l,5,m] <- sd(test[[m]]$beta[j,k,l,])
          # results[j,k,l,6,m] <- test[[m]]$time[3]
          results[j,k,l,7,m] <- effectiveSize(test[[m]]$beta[j,k,l,])
        }
      }
    }
  }
  # slopes do not vary by frequency
  test[[5]] <- MCMC_model3_2(Ystar, Xstar, Q = Qmat, p = E, n = n, r = R, q = R, iters = iters, burn = burn, 1, n)
  # one-step spatial plus
  test[[6]] <- MCMC_model3_2(Ystar, Xstar, Q = Qmat, p = E, n = n, r = R, q = R, iters = iters, burn = burn, num1, num2)
  test[[7]] <- MCMC_model3_2(Ystar, Xstar, Q = Qmat, p = E, n = n, r = R, q = R, iters = iters, burn = burn, num1, num2/2)
  for (k in 1:E){
    for (l in 1:R){
      for (m in 5:7){
        results[1,k,l,1,m] <- mean(test[[m]]$beta[k,l,])
        results[1,k,l,2,m] <- quantile(test[[m]]$beta[k,l,], probs = c(0.025))
        results[1,k,l,3,m] <- quantile(test[[m]]$beta[k,l,], probs = c(0.5))
        results[1,k,l,4,m] <- quantile(test[[m]]$beta[k,l,], probs = c(0.975))
        results[1,k,l,5,m] <- sd(test[[m]]$beta[k,l,])
        results[1,k,l,6,m] <- test[[m]]$time[3]
        results[1,k,l,7,m] <- effectiveSize(test[[m]]$beta[k,l,])
      }
    }
  }
  # the next three models are univariate models with horse prior
  # note the WAIC for each R is different
  # should sum the WAIC for each R 
  tol_time <- 0
  for (r in 1:R){
    test[[r+7]] <- MCMC_model4_9_v2(Ystar[,r], Xstar, Qmat = Qmat, E = E, n = n, L, K, B = B, iters = iters, burn = burn)
    tol_time <- tol_time + test[[r+7]]$time[3]
  }
  for (j in 1:n){
    for (k in 1:E){
      for (l in 1:R){
        for (m in 8:8){
          results[j,k,l,1,m] <- mean(test[[(l+7)]]$beta[j,k,])
          results[j,k,l,2,m] <- quantile(test[[(l+7)]]$beta[j,k,], probs = c(0.025))
          results[j,k,l,3,m] <- quantile(test[[(l+7)]]$beta[j,k,], probs = c(0.5))
          results[j,k,l,4,m] <- quantile(test[[(l+7)]]$beta[j,k,], probs = c(0.975))
          results[j,k,l,5,m] <- sd(test[[(l+7)]]$beta[j,k,])
          results[j,k,l,6,m] <- tol_time
          results[j,k,l,7,m] <- effectiveSize(test[[(l+7)]]$beta[j,k,])
        }
      }
    }
  }
  return(results)
}

#####################################################################
######### Functions for metrics, bias, and results ##################
#####################################################################


# return scaled score for each M,E,R
# need to separate 0 and non-0 coefficients
metrics_over_all <- function(sim, true, M, groupE){
  E <- nrow(true) 
  R <- ncol(true) 
  MSE1 <- bias1 <- cov1 <- sd1 <- array(0, c(M,groupE,R))
  MSE1[,1,] <- bias1[,1,] <- cov1[,1,] <- sd1[,1,] <- NA
  MSE0 <- bias0 <- cov0 <- sd0 <- array(0, c(M,(E-groupE),R))
  for (m in 1:M){
    for (i in 2:groupE){
      for (j in 1:R){
        MSE1[m,i,j] <- mean((sim[1,i,j,1,m,] - true[i,j])^2)
        bias1[m,i,j] <- mean((sim[1,i,j,1,m,] - true[i,j]))
        cov1[m,i,j] <- mean((sim[1,i,j,2,m,]<true[i,j])*(sim[1,i,j,4,m,]>true[i,j]))
        sd1[m,i,j] <- mean(sim[1,i,j,5,m,])
      }
    }
  }
  for (m in 1:M){
    for (i in (groupE+1):E){
      for (j in 1:R){
        MSE0[m,(i-groupE),j] <- mean((sim[1,i,j,1,m,] - true[i,j])^2)
        bias0[m,(i-groupE),j] <- mean((sim[1,i,j,1,m,] - true[i,j]))
        cov0[m,(i-groupE),j] <- mean((sim[1,i,j,2,m,]<true[i,j])*(sim[1,i,j,4,m,]>true[i,j]))
        sd0[m,(i-groupE),j] <- mean(sim[1,i,j,5,m,])
      }
    }
  }
  MSE_ave1 <- rowMeans(MSE1, na.rm = TRUE, dims = 1)
  MSE_ave0 <- rowMeans(MSE0, na.rm = TRUE, dims = 1)  
  bias_ave1 <- rowMeans(abs(bias1), na.rm = TRUE, dims = 1)
  bias_ave0 <- rowMeans(abs(bias0), na.rm = TRUE, dims = 1)
  cov_ave1 <- rowMeans(cov1, na.rm = TRUE, dims = 1)
  cov_ave0 <- rowMeans(cov0, na.rm = TRUE, dims = 1)
  sd_ave1 <- rowMeans(sd1, na.rm = TRUE, dims = 1)
  sd_ave0 <- rowMeans(sd0, na.rm = TRUE, dims = 1)
  out <- list(MSE_ave1 = MSE_ave1, MSE_ave0 = MSE_ave0, bias_ave1 = bias_ave1, bias_ave0 = bias_ave0, cov_ave1 = cov_ave1, cov_ave0 = cov_ave0, sd_ave1 = sd_ave1, sd_ave0 = sd_ave0)
  return(out)
}

# return scaled score for each M,E,R
# separate 0 and non-0 coefficients
# added length
# mean absolute deviation
metrics_over_all2 <- function(sim, true, M, groupE){
  E <- nrow(true) 
  R <- ncol(true) 
  MSE1 <- bias1 <- cov1 <- sd1 <- length1 <- array(0, c(M,groupE,R))
  MSE1[,1,] <- bias1[,1,] <- cov1[,1,] <- sd1[,1,] <- length1[,1,] <- NA
  MSE0 <- bias0 <- cov0 <- sd0 <- length0 <- array(0, c(M,(E-groupE),R))
  for (m in 1:M){
    for (i in 2:groupE){
      for (j in 1:R){
        MSE1[m,i,j] <- mean((sim[1,i,j,1,m,] - true[i,j])^2)
        bias1[m,i,j] <- mean((sim[1,i,j,1,m,] - true[i,j]))
        cov1[m,i,j] <- mean((sim[1,i,j,2,m,]<true[i,j])*(sim[1,i,j,4,m,]>true[i,j]))
        sd1[m,i,j] <- mean(sim[1,i,j,5,m,])
        length1[m,i,j] <- mean(sim[1,i,j,4,m,] - sim[1,i,j,2,m,])
      }
    }
  }
  for (m in 1:M){
    for (i in (groupE+1):E){
      for (j in 1:R){
        MSE0[m,(i-groupE),j] <- mean((sim[1,i,j,1,m,] - true[i,j])^2)
        bias0[m,(i-groupE),j] <- mean((sim[1,i,j,1,m,] - true[i,j]))
        cov0[m,(i-groupE),j] <- mean((sim[1,i,j,2,m,]<true[i,j])*(sim[1,i,j,4,m,]>true[i,j]))
        sd0[m,(i-groupE),j] <- mean(sim[1,i,j,5,m,])
        length0[m,(i-groupE),j] <- mean(sim[1,i,j,4,m,] - sim[1,i,j,2,m,])
      }
    }
  }
  MSE_ave1 <- rowMeans(MSE1, na.rm = TRUE, dims = 1)
  MSE_ave0 <- rowMeans(MSE0, na.rm = TRUE, dims = 1)  
  # bias_ave1 <- rowMeans(abs(bias1), na.rm = TRUE, dims = 1)
  # bias_ave0 <- rowMeans(abs(bias0), na.rm = TRUE, dims = 1)
  bias_ave1 <- rowMeans(bias1, na.rm = TRUE, dims = 1)
  bias_ave0 <- rowMeans(bias0, na.rm = TRUE, dims = 1)
  cov_ave1 <- rowMeans(cov1, na.rm = TRUE, dims = 1)
  cov_ave0 <- rowMeans(cov0, na.rm = TRUE, dims = 1)
  sd_ave1 <- rowMeans(sd1, na.rm = TRUE, dims = 1)
  sd_ave0 <- rowMeans(sd0, na.rm = TRUE, dims = 1)
  length_ave1 <- rowMeans(length1, na.rm = TRUE, dims = 1)
  length_ave0 <- rowMeans(length0, na.rm = TRUE, dims = 1)
  out <- list(MSE_ave1 = MSE_ave1, MSE_ave0 = MSE_ave0, bias_ave1 = bias_ave1, bias_ave0 = bias_ave0, cov_ave1 = cov_ave1, cov_ave0 = cov_ave0, sd_ave1 = sd_ave1, sd_ave0 = sd_ave0, length_ave1 = length_ave1, length_ave0 = length_ave0)
  return(out)
}

# return scaled score for each M,E,R
# need to separate 0 and non-0 coefficients
se_over_all <- function(sim, true, M, groupE){
  E <- nrow(true) 
  R <- ncol(true) 
  dim(sim[[1]])[6]
  MSE_se1 <- bias_se1 <- cov_se1 <- sd_se1 <- array(0, c(M,groupE,R))
  MSE_se1[,1,] <- bias_se1[,1,] <- cov_se1[,1,] <- sd_se1[,1,] <- NA
  MSE_se0 <- bias_se0 <- cov_se0 <- sd_se0 <- array(0, c(M,(E-groupE),R))
  for (m in 1:M){
    for (i in 2:groupE){
      for (j in 1:R){
        MSE_se1[m,i,j] <- sd((sim[1,i,j,1,m,] - true[i,j])^2)
        bias_se1[m,i,j] <- sd((sim[1,i,j,1,m,] - true[i,j]))
        cov_se1[m,i,j] <- sd((sim[1,i,j,2,m,]<true[i,j])*(sim[1,i,j,4,m,]>true[i,j]))
        sd_se1[m,i,j] <- sd(sim[1,i,j,5,m,])
      }
    }
  }
  for (m in 1:M){
    for (i in (groupE+1):E){
      for (j in 1:R){
        MSE_se0[m,(i-groupE),j] <- sd((sim[1,i,j,1,m,] - true[i,j])^2)
        bias_se0[m,(i-groupE),j] <- sd((sim[1,i,j,1,m,] - true[i,j]))
        cov_se0[m,(i-groupE),j] <- sd((sim[1,i,j,2,m,]<true[i,j])*(sim[1,i,j,4,m,]>true[i,j]))
        sd_se0[m,(i-groupE),j] <- sd(sim[1,i,j,5,m,])
      }
    }
  }
  MSE_se_ave1 <- rowMeans(MSE_se1, na.rm = TRUE, dims = 1)/sqrt(reps)
  MSE_se_ave0 <- rowMeans(MSE_se0, na.rm = TRUE, dims = 1)/sqrt(reps)  
  bias_se_ave1 <- rowMeans(bias_se1, na.rm = TRUE, dims = 1)/sqrt(reps)
  bias_se_ave0 <- rowMeans(bias_se0, na.rm = TRUE, dims = 1)/sqrt(reps)
  cov_se_ave1 <- rowMeans(cov_se1, na.rm = TRUE, dims = 1)/sqrt(reps)
  cov_se_ave0 <- rowMeans(cov_se0, na.rm = TRUE, dims = 1)/sqrt(reps)
  sd_se_ave1 <- rowMeans(sd_se1, na.rm = TRUE, dims = 1)/sqrt(reps)
  sd_se_ave0 <- rowMeans(sd_se0, na.rm = TRUE, dims = 1)/sqrt(reps)
  out <- list(MSE_se_ave1 = MSE_se_ave1, MSE_se_ave0 = MSE_se_ave0, bias_se_ave1 = bias_se_ave1, bias_se_ave0 = bias_se_ave0, cov_se_ave1 = cov_se_ave1, cov_se_ave0 = cov_se_ave0, sd_se_ave1 = sd_se_ave1, sd_se_ave0 = sd_se_ave0)
  return(out)
}

# return scaled score for each M,E,R
# need to separate 0 and non-0 coefficients
# added length
se_over_all2 <- function(sim, true, M, groupE){
  E <- nrow(true) 
  R <- ncol(true) 
  dim(sim[[1]])[6]
  MSE_se1 <- bias_se1 <- cov_se1 <- sd_se1 <- length_se1 <- array(0, c(M,groupE,R))
  MSE_se1[,1,] <- bias_se1[,1,] <- cov_se1[,1,] <- sd_se1[,1,] <- length_se1[,1,] <- NA
  MSE_se0 <- bias_se0 <- cov_se0 <- sd_se0 <- length_se0 <- array(0, c(M,(E-groupE),R))
  for (m in 1:M){
    for (i in 2:groupE){
      for (j in 1:R){
        MSE_se1[m,i,j] <- sd((sim[1,i,j,1,m,] - true[i,j])^2)
        bias_se1[m,i,j] <- sd((sim[1,i,j,1,m,] - true[i,j]))
        cov_se1[m,i,j] <- sd((sim[1,i,j,2,m,]<true[i,j])*(sim[1,i,j,4,m,]>true[i,j]))
        sd_se1[m,i,j] <- sd(sim[1,i,j,5,m,])
        length_se1[m,i,j] <- sd(sim[1,i,j,4,m,] - sim[1,i,j,2,m,])
      }
    }
  }
  for (m in 1:M){
    for (i in (groupE+1):E){
      for (j in 1:R){
        MSE_se0[m,(i-groupE),j] <- sd((sim[1,i,j,1,m,] - true[i,j])^2)
        bias_se0[m,(i-groupE),j] <- sd((sim[1,i,j,1,m,] - true[i,j]))
        cov_se0[m,(i-groupE),j] <- sd((sim[1,i,j,2,m,]<true[i,j])*(sim[1,i,j,4,m,]>true[i,j]))
        sd_se0[m,(i-groupE),j] <- sd(sim[1,i,j,5,m,])
        length_se0[m,(i-groupE),j] <- sd(sim[1,i,j,4,m,] - sim[1,i,j,2,m,])
      }
    }
  }
  MSE_se_ave1 <- rowMeans(MSE_se1, na.rm = TRUE, dims = 1)/sqrt(reps)
  MSE_se_ave0 <- rowMeans(MSE_se0, na.rm = TRUE, dims = 1)/sqrt(reps)  
  bias_se_ave1 <- rowMeans(bias_se1, na.rm = TRUE, dims = 1)/sqrt(reps)
  bias_se_ave0 <- rowMeans(bias_se0, na.rm = TRUE, dims = 1)/sqrt(reps)
  cov_se_ave1 <- rowMeans(cov_se1, na.rm = TRUE, dims = 1)/sqrt(reps)
  cov_se_ave0 <- rowMeans(cov_se0, na.rm = TRUE, dims = 1)/sqrt(reps)
  sd_se_ave1 <- rowMeans(sd_se1, na.rm = TRUE, dims = 1)/sqrt(reps)
  sd_se_ave0 <- rowMeans(sd_se0, na.rm = TRUE, dims = 1)/sqrt(reps)
  length_se_ave1 <- rowMeans(length_se1, na.rm = TRUE, dims = 1)/sqrt(reps)
  length_se_ave0 <- rowMeans(length_se0, na.rm = TRUE, dims = 1)/sqrt(reps)
  out <- list(MSE_se_ave1 = MSE_se_ave1, MSE_se_ave0 = MSE_se_ave0, bias_se_ave1 = bias_se_ave1, bias_se_ave0 = bias_se_ave0, cov_se_ave1 = cov_se_ave1, cov_se_ave0 = cov_se_ave0, sd_se_ave1 = sd_se_ave1, sd_se_ave0 = sd_se_ave0, length_se_ave1 = length_se_ave1, length_se_ave0 = length_se_ave0)
  return(out)
}

combine_sim_res <- function(metric, m, namez){
  res <- met[[1]][[metric]]
  for (i in 2:m){
    res <- cbind(res, met[[i]][[metric]])
  }
  colnames(res) <- paste0("setting ", 1:m)
  rownames(res) <- namez
  return(res)
}

combine_sim <- function(list, metric, m, namez){
  res <- list[[1]][[metric]]
  for (i in 2:m){
    res <- cbind(res, list[[i]][[metric]])
  }
  colnames(res) <- paste0("setting ", 1:m)
  rownames(res) <- namez
  return(res)
}

flatten <- function(sims, reps, m, n, dim, E, R){
  results <- array(0, c(n, E, R, dim, m, reps))
  for (i in 1:reps){
    results[,,,,,i] <- sims[[i]]
  }
  return(results)
}
