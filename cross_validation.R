# full MV method with shrinkage prior 
# allow for covariates (Z) and treatments (X)
# cross validation
MCMC_model4_8_XZ_cv <- function(Y, X, Z, Qmat, E, n, R, Q, L, K, B, iters = 1000, burn = 50, full_result = FALSE){
  tik <- proc.time()
  # Aindex keeps track of the number of entries in each row of A to be updated
  Aindex <- rep(0, R)
  for (r in 1:R){
    if (r <= Q) Aindex[r] <- r-1
    else Aindex[r] <- Q
  }
  keepers_beta <- array(0, c(n, E, R, iters))
  p <- ncol(Z)
  keepers_alpha <- array(0, c(p, R, iters))
  Amat <- matrix(0, nrow = R, ncol = Q)
  for (r in 1:R){
    if (r <= Q) Amat[r,r] <- 1    
  }
  att <- acc <- 0
  # priors and initial values
  MH <- 0.1
  W<- eigen(Qmat)$values
  a_tau <- .5
  b_tau <- .005
  Umat <- matrix(0, nrow = n, ncol = Q)
  Tl <- matrix(0, nrow = L, ncol = K)
  Te <- matrix(0, nrow = E, ncol = K)
  Tr <- matrix(0, nrow = R, ncol = K)
  sig2 <- rep(0.5, Q)
  tau2 <- rep(0.5, R)
  lambda <- 0.5
  keepers_varz <- matrix(0.5, nrow = iters, ncol = R+1)
  keepers_AU <- array(0, c(n, R, iters))
  lambdaKHC <- rep(1, K)
  lambdaEHC <- rep(1, E)
  lambdaRHC <- rep(1, R)
  nuRHC <- rep(1, R)
  tauKHC <- tauEHC <- tauHC <- tauQHC <- 1
  xiHC <- 1
  X_tensor <- array(X, c(n, E, R))
  estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
  alpha <- keepers_alpha[,,1]
  
  # fill in values for test set 
  miss    <- which(is.na(Y[,1]))
  nm      <- length(miss)
  for (r in 1:ncol(Y)){
    Y[miss,r] <- rnorm(nm) 
  }
  # running mean and var for predicted Y test set
  Moment1 <- Moment2 <- 0
  
  for (iter in 2:iters){
    # update tensor margins for beta
    AU <- Umat%*%t(Amat)
    Zalpha <- Z%*%alpha
    for (k in 1:K){
      for (l in 1:L){
        Cl <- B[,l] * (X %*% outer(Te[,k], Tr[,k]))
        Vl <- sum(colSums(Cl^2)/tau2) + 1/(lambdaKHC[k])
        # Vl <- sum(colSums(Cl^2)/tau2) + 1
        Tl_prev <- Tl[l,k]
        Ml <- sum(colSums(Cl*(Y - AU - Zalpha - estimate + Tl_prev*Cl))/tau2)
        Tl[l,k] <- rnorm(1, Ml/Vl, 1/sqrt(Vl))
        # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
        estimate <- estimate + (Tl[l,k] - Tl_prev)*Cl 
      }
      for (e in 1:E){
        Ce <- X[,e] * (B %*% outer(Tl[,k], Tr[,k]))
        Ve <- sum(colSums(Ce^2)/tau2) + 1/(lambdaEHC[e])
        # Ve <- sum(colSums(Ce^2)/tau2) + 1
        Te_prev <- Te[e,k]
        Me <- sum(colSums(Ce*(Y - AU - Zalpha - estimate + Te_prev*Ce))/tau2)
        Te[e,k] <- rnorm(1, Me/Ve, 1/sqrt(Ve))
        # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
        estimate <- estimate + (Te[e,k] - Te_prev)*Ce
      }
      for (r in 1:R){
        Cr <- rowSums( X * (B %*% outer(Tl[,k], Te[,k])))
        Vr <- sum(Cr^2)/tau2[r] + 1/(lambdaRHC[r]*tauHC*tau2[r])
        # Vr <- sum(Cr^2)/tau2[r] + 1
        Tr_prev <- Tr[r,k]
        Mr <- sum(Cr*(Y[,r] - AU[,r] - Zalpha[,r] - estimate[,r] + Tr_prev*Cr)/tau2[r])
        Tr[r,k] <- rnorm(1, Mr/Vr, 1/sqrt(Vr))
        # estimate <- find_estimate(X_tensor, B, Tl, Te, Tr, n, R)
        estimate[,r] <- estimate[,r] + (Tr[r,k] - Tr_prev)*Cr
      }
      # update lambdaKHC[k]
      # sum_tauHC[k] <- sum(Tr[,k]^2/(2*lambdaRHC*tau2))
      lambdaKHC[k] <- rinvgamma(1, (L+1)/2, sum(Tl[,k]^2)/2 + 1/tauKHC)
    }
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
      tau2[i] <- rinvgamma(1, (n+K)/2 + a_tau, 0.5*sum((Y[,i]-Zalpha[,i]-estimate[,i]-Umat%*%Amat[i,])^2) + sumbyr[i]/(lambdaRHC[i]*tauHC) + b_tau)
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
    
    AU <- Umat%*%t(Amat)
    Zalpha <- Z%*%alpha
    yhat <- estimate + AU + Zalpha
    
    # missing values
    if(nm>0){
      for (r in 1:R){
        Y[miss,r] <- rnorm(nm, yhat[miss,r], sqrt(tau2[r])) 
      }
    }
    
    if(iter>burn){
      Moment1 <- Moment1 + Y/(iters-burn)
      Moment2 <- Moment2 + Y*Y/(iters-burn)
    }
    
    keepers_beta[,,,iter] <- bbygamma(B, find_tol_gamma(Tl, Te, Tr), n)
    keepers_alpha[,,iter] <- alpha
    keepers_AU[,,iter] <- AU
    keepers_varz[iter,] <- c(tau2, lambda)
    
  }

  tok <- proc.time()
  if (full_result){
    out <- list(beta = keepers_beta[,,,(burn+1):iters], alpha = keepers_alpha[,,(burn+1):iters], varz = keepers_varz, AU = keepers_AU, acc_rate = acc/att, time = tok - tik, MH = MH, y_mean = Moment1, y_var  = Moment2-Moment1^2)  
  }
  else{
    out <- list(beta = keepers_beta[1,,,(burn+1):iters], alpha = keepers_alpha[,,(burn+1):iters], varz = keepers_varz, AU = keepers_AU, acc_rate = acc/att, time = tok - tik, MH = MH, y_mean = Moment1, y_var  = Moment2-Moment1^2)  
  }
  
  return(out)
}


# full model with inverse gamma prior for beta
# lambda1 \sim invgamma(1/2, 1/gamma1) # this is lambdaKHC
# lambda2 \sim invgamma(1/2, 1/gamma2) # this is lambdaEHC
# lambda3 (lambdaRHC) is the same as tauHC
# allow for covariates (Z) and treatments (X)
MCMC_model4_8_3_XZ_cv <- function(Y, X, Z, Qmat, E, n, R, Q, L, K, B, iters = 1000, burn = 50, W, thin = 2){
  tik <- proc.time()
  # Aindex keeps track of the number of entries in each row of A to be updated
  Aindex <- rep(0, R)
  for (r in 1:R){
    if (r <= Q) Aindex[r] <- r-1
    else Aindex[r] <- Q
  }
  iters_thin <- (iters-burn)%/%thin
  p <- ncol(Z)
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
  
  AU <- Umat%*%t(Amat)
  Zalpha <- Z%*%alpha
  
  # fill in values for test set 
  miss    <- which(is.na(Y[,1]))
  nm      <- length(miss)
  for (r in 1:ncol(Y)){
    Y[miss,r] <- rnorm(nm) 
  }
  # running mean and var for predicted Y test set
  Moment1 <- Moment2 <- 0
  
  for (iter in 2:iters){
    print(paste0("This is the ", iter, "th iterations."))
    # update tensor margins for beta
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
    
    AU <- Umat%*%t(Amat)
    Zalpha <- Z%*%alpha
    yhat <- estimate + AU + Zalpha
    
    # missing values
    if(nm>0){
      for (r in 1:R){
        Y[miss,r] <- rnorm(nm, yhat[miss,r], sqrt(tau2[r])) 
      }
    }
    
    if(iter>burn){
      Moment1 <- Moment1 + Y/(iters-burn)
      Moment2 <- Moment2 + Y*Y/(iters-burn)
    }
    
    if (iter > burn){
      if ((iter-burn)%%thin == 0){
        # keepers_beta[,,,(iter-burn)%/%thin] <- bbygamma(B, find_tol_gamma(Tl, Te, Tr), n)
        keepers_Tl[,,(iter-burn)%/%thin] <- Tl
        keepers_Te[,,(iter-burn)%/%thin] <- Te
        keepers_Tr[,,(iter-burn)%/%thin] <- Tr
      }
    }
  }
  
  tok <- proc.time()
  out <- list(Tl = keepers_Tl, Te = keepers_Te, Tr = keepers_Tr, acc_rate = acc/att, time = tok - tik, MH = MH, y_mean = Moment1, y_var  = Moment2-Moment1^2)  
  return(out)
}

