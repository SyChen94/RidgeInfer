
## DDL
estimate_coefficients <- function(X, Y, rho = 0.5) {
    # perform Single Val Decomp on X
    UDV_list <- svd(X)
    U <- UDV_list$u
    D <- UDV_list$d
    V <- t(UDV_list$v)

    # forms QX and QY(STEP 1)
    tau <- quantile(D, rho)
    Dtilde <- pmin(D, tau)
    Q <- diag(nrow(X)) - U %*% diag(1 - Dtilde / D) %*% t(U)
    Xtilde <- Q %*% X
    Ytilde <- Q %*% Y

    # perform a lasso fit for trimmed QX and QY(STEP 2)
    fit <- glmnet::cv.glmnet(x = Xtilde, y = Ytilde)

    # B-init
    betahat <- as.matrix((coef(fit, S = fit$lambda.min)[-1]))

    # determine bhat
    res <- Y - X %*% betahat
    fit <- glmnet::cv.glmnet(x = U %*% diag(D), y = res, alpha = 1)
    bhat <- t(V) %*% coef(fit, s = fit$lambda.min)[-1]

    return_listcoeff <- list("betahat" = betahat, "bhat" = bhat)
    return(return_listcoeff)
}
estimate_coefficients1 <- function(X, Y, rho = 0.5) {
  # perform Single Val Decomp on X
  UDV_list <- svd(X)
  U <- UDV_list$u
  D <- UDV_list$d
  V <- t(UDV_list$v)
  
  # forms QX and QY(STEP 1)
  tau <- quantile(D, rho)
  Dtilde <- pmin(D, tau)
  #Q <- diag(nrow(X)) - U %*% diag(1 - Dtilde / D) %*% t(U)
  Q=diag(nrow(X))
  Xtilde <- Q %*% X
  Ytilde <- Q %*% Y
  # perform a lasso fit for trimmed QX and QY(STEP 2)
  fit <- glmnet::cv.glmnet(x = Xtilde, y = Ytilde)
  
  # B-init
  betahat <- as.matrix((coef(fit, S = fit$lambda.min)[-1]))
  
  # determine bhat
  res <- Y - X %*% betahat
  fit <- glmnet::cv.glmnet(x = U %*% diag(D), y = res, alpha = 1)
  bhat <- t(V) %*% coef(fit, s = fit$lambda.min)[-1]
  
  return_listcoeff <- list("betahat" = betahat, "bhat" = bhat)
  return(return_listcoeff)
}


find_z <- function(X, index) {
    n <- dim(X)[1]
    p <- dim(X)[2]

    # Xj and X-j
    X_j <- X[, index]
    X_negj <- X[, -index]

    # regress X-j on xj, use least min lambda to estimate gamma(Step 4)
    cvfit <- glmnet::cv.glmnet(x = X_negj, y = X_j)
    gamma <- coef(cvfit, s = cvfit$lambda.min)[-1]
    # eq. (8), residuals(Step 5)
    z <- n^-0.5 * (X_j - X_negj %*% gamma)

    # variation from eq. (23) with 25% increase(read 3.6)
    V <- 1.25 * n^0.5 * norm(z, type = "2") / (t(z) %*% X_j)

    # take first z whose variance is at most 25% larger than for the CV lambda
    for (lam in cvfit$glmnet.fit$lambda) {
        gamma <- coef(cvfit, s = lam)[-1]
        z <- n^(-0.5) * (X_j - X_negj %*% matrix(gamma, p - 1, 1))
        if (n^0.5 * (norm(z, type = "2") / (t(z) %*% X_j)) > V) {
            break
        }
    }

    # normalize Z with 2 norm
    z <- z / norm(z, type = "2")
    return(z)
}



estimate_sigma <- function(X, Y, rho = 0.5, alt = FALSE, active_set_scaling = FALSE) {
    # uses both fitting estimators to create an unbiased estimator of sigma
    est_coef <- estimate_coefficients(X, Y, rho)
    betahat <- est_coef$betahat
    bhat <- est_coef$bhat

    UDV_list <- svd(X)
    U <- UDV_list$u
    D <- UDV_list$d
    V <- t(UDV_list$v)

    tau <- quantile(D, rho)
    Dtilde <- pmin(D, tau)
    # The length of D has length min(n,p). In order to avoid the over-shrinkage to the last n-p eigenvector
    # when n>p, the following is necessary.
    Q <- diag(nrow(X)) - U %*% diag(1 - Dtilde / D) %*% t(U)
    Xtilde <- Q %*% X
    Ytilde <- Q %*% Y

    # Step 7, two methods of sigmahat computation, the first has shown more robust results
    divisor <- sum(diag(Q %*% Q))
    error <- (norm(Ytilde - Xtilde %*% betahat, type = "2"))^2
    sigmahat <- (error / divisor)^0.5
    if (active_set_scaling) {
        sigmahat <- (nrow(X) / (nrow(X) - min(nrow(X) / 2, Matrix::nnzero(betahat))))^0.5 * sigmahat
    }
    if (alt) {
        residuals <- Y - X %*% betahat - X %*% bhat
        sigmahat <- mean(residuals^2)^0.5
    }

    return(sigmahat)
}


estimate_sigma1 <- function(X, Y, rho = 0.5, alt = FALSE, active_set_scaling = FALSE) {
  # uses both fitting estimators to create an unbiased estimator of sigma
  est_coef1 <- estimate_coefficients1(X, Y, rho)
  betahat <- est_coef1$betahat
  bhat <- est_coef1$bhat
  
  # UDV_list <- svd(X)
  # U <- UDV_list$u
  # D <- UDV_list$d
  # V <- t(UDV_list$v)
  
  # tau <- quantile(D, rho)
  # Dtilde <- pmin(D, tau)
  # The length of D has length min(n,p). In order to avoid the over-shrinkage to the last n-p eigenvector
  # when n>p, the following is necessary.
  #Q <- diag(nrow(X)) - U %*% diag(1 - Dtilde / D) %*% t(U)
  Q <- diag(nrow(X))
  Xtilde <- Q %*% X
  Ytilde <- Q %*% Y
  
  # Step 7, two methods of sigmahat computation, the first has shown more robust results
  divisor <- sum(diag(Q %*% Q))
  error <- (norm(Ytilde - Xtilde %*% betahat, type = "2"))^2
  sigmahat <- (error / divisor)^0.5
  if (active_set_scaling) {
    sigmahat <- (nrow(X) / (nrow(X) - min(nrow(X) / 2, Matrix::nnzero(betahat))))^0.5 * sigmahat
  }
  if (alt) {
    residuals <- Y - X %*% betahat - X %*% bhat
    sigmahat <- mean(residuals^2)^0.5
  }
  
  return(sigmahat)
}


DDL1 <- function(X, Y, index, rho , rhop ,sig_level = 0.05) {
    # determines parameters
    n <- dim(X)[1]
    p <- dim(X)[2]
    dblasso <- rep(NA, length(index))
    stddev <- rep(NA, length(index))
    lower <- rep(NA, length(index))
    upper <- rep(NA, length(index))
    liein <- rep(NA, length(index))
    bias <- rep(NA, length(index))
    cover0 <- rep(NA, length(index))
    est <- estimate_coefficients(X, Y, rho)
    betahat <- est$betahat
    i=1
    for (i in seq(length(index))) {
        X_negj <- X[, -index[i]]
        # single value decomposition of X (Trim transform)
        UDV_list <- svd(X_negj)
        U <- UDV_list$u
        D <- UDV_list$d
        V <- t(UDV_list$v)

        # Reduce D, then make P trim
        Dtilde <- pmin(D, quantile(D, rhop))
        P <- diag(nrow(X_negj)) - U %*% diag(1 - Dtilde / D) %*% t(U)
        P_X <- P %*% X

        # determine projection direction then estimate betahat and bhat(the z here has been multiplied with P)
        z <- find_z(P_X, index[i])
        # bhat = est$bhat

        # eq. (12) for a point estimation of Bj(Step 6)
        dblasso[i] <- t(z) %*% P %*% (Y - X_negj %*% betahat[-index[i]]) / (t(z) %*% P %*% X[, index[i]])

        # eq. (23) for
        Variance <- (t(z) %*% (P^2) %*% z) / (t(z) %*% P %*% X[, index[i]])^2

        # (Step 7)
        sigmahat <- estimate_sigma(X, Y, rho)
        stddev[i] <- sigmahat * Variance^.5
        # 2-sided CI critical value for sig_level
        quantile <- qnorm(1 - sig_level / 2, mean = 0, sd = 1, lower.tail = TRUE)
        # #determine bounds of interval, from eq. (13) (Step 8)
        lower[i] <- dblasso[i] - quantile * stddev[i]
        upper[i] <- dblasso[i] + quantile * stddev[i]
        # bias[i] <- abs(dblasso[i] - true[i])
        # liein[i] <- ifelse((true[i] > lower[i]) && (true[i] < upper[i]), 1, 0)
        # B_b=t(z) %*% P_X %*% b/(t(z) %*% P_X[,index]*stddev)
        cover0[i] <- ifelse(lower[i] * upper[i] > 0, 0, 1)
        # B_beta= t(z) %*% P_X[,-index] %*% (beta[-index]-betahat[-index])/(t(z) %*% P_X[,index]*stddev)
    }
    rs <- list(
        est_ddl = dblasso,
        se = stddev,
        sigmahat = sigmahat,
        # bi = bias,
        len = upper - lower,
         CI=c(lower,upper),
        # li = liein,
        co0 = cover0
    )
    return(rs)
}

### debias without hidden
DLasso <- function(X, Y, index, rho , rhop ,sig_level = 0.05) {
  # determines parameters
  n <- dim(X)[1]
  p <- dim(X)[2]
  dblasso <- rep(NA, length(index))
  stddev <- rep(NA, length(index))
  lower <- rep(NA, length(index))
  upper <- rep(NA, length(index))
  liein <- rep(NA, length(index))
  bias <- rep(NA, length(index))
  cover0 <- rep(NA, length(index))
  est <- estimate_coefficients1(X, Y, rho)
  betahat <- est$betahat
  i=1
  for (i in seq(length(index))) {
    X_negj <- X[, -index[i]]
    # single value decomposition of X (Trim transform)
    UDV_list <- svd(X_negj)
    U <- UDV_list$u
    D <- UDV_list$d
    V <- t(UDV_list$v)
    
    Dtilde <- pmin(D, quantile(D, rhop))
    #P <- diag(nrow(X_negj)) - U %*% diag(1 - Dtilde / D) %*% t(U)
    P=diag(nrow(X))
    P_X <- P %*% X
    
    # determine projection direction then estimate betahat and bhat(the z here has been multiplied with P)
    z <- find_z(P_X, index[i])
    # bhat = est$bhat
    
    # eq. (12) for a point estimation of Bj(Step 6)
    dblasso[i] <- t(z) %*% P %*% (Y - X_negj %*% betahat[-index[i]]) / (t(z) %*% P %*% X[, index[i]])
    
    # eq. (23) for
    Variance <- (t(z) %*% (P^2) %*% z) / (t(z) %*% P %*% X[, index[i]])^2
    
    # (Step 7)
    sigmahat <- estimate_sigma1(X, Y, rho)
    stddev[i] <- sigmahat * Variance^.5
    # 2-sided CI critical value for sig_level
    quantile <- qnorm(1 - sig_level / 2, mean = 0, sd = 1, lower.tail = TRUE)
    # #determine bounds of interval, from eq. (13) (Step 8)
    lower[i] <- dblasso[i] - quantile * stddev[i]
    upper[i] <- dblasso[i] + quantile * stddev[i]
    # bias[i] <- abs(dblasso[i] - true[i])
    # liein[i] <- ifelse((true[i] > lower[i]) && (true[i] < upper[i]), 1, 0)
    # B_b=t(z) %*% P_X %*% b/(t(z) %*% P_X[,index]*stddev)
    cover0[i] <- ifelse(lower[i] * upper[i] > 0, 0, 1)
    # B_beta= t(z) %*% P_X[,-index] %*% (beta[-index]-betahat[-index])/(t(z) %*% P_X[,index]*stddev)
  }
  rs <- list(
    est_ddl = dblasso,
    se = stddev,
    sigmahat = sigmahat,
    # bi = bias,
    len = upper - lower,
    CI=c(lower,upper),
    # li = liein,
    co0 = cover0
  )
  return(rs)
}

### DLassoAUC
DLassoAUC<-function(X, Y, index, K.factors, rho, rhop, sig_level = 0.05){
  output = farm.res(X,K.factors=K.factors,robust=FALSE) #default options
  U=output$factors
  X=cbind(X,U)
  resu=DLasso(X, Y, index, rho , rhop ,sig_level = 0.05)
  return(resu)
}
