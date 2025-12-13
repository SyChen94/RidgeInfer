library(glmnet)



myDLasso0222 <- function(X, Y, index, rho , rhop ,sig_level = 0.05) {
  # determines parameters
  n <- dim(X)[1]
  p <- dim(X)[2]
  dblasso <- rep(NA, length(index))
  stddev <- rep(NA, length(index))
  lower <- rep(NA, length(index))
  upper <- rep(NA, length(index))
  cover0 <- rep(NA, length(index))
  # est <- my_estglmnet(X, Y)
  ###################################################################
  ###################################################################
  ###################################################################
  ###################################################################
  lam_min = 0.1*max(abs(t(X)%*%Y)/(n*p))
  lam_max = 40*max(abs(t(X)%*%Y)/(n*p))
  the_lambdas = seq(lam_min,lam_max,(lam_max - lam_min)/20)
  fit <- glmnet::cv.glmnet(x = X, y = Y,alpha=1,lambda=the_lambdas)
  betahat <- as.matrix((coef(fit, S = fit$lambda.min)[-1]))
  # betahat <- as.matrix((coef(fit, S = fit$lambda.1se)[-1]))
  # fit <- glmnet::cv.glmnet(x = X, y = Y,alpha=1)
  # betahat <- as.matrix((coef(fit, S = fit$lambda.min)[-1]))

  ###################################################################
  ###################################################################
  ###################################################################
  ###################################################################
  i=1
  for (i in seq(length(index))) {
    X_negj <- X[, -index[i]]
    # z <- myfind_z(P_X, index[i])
    ######################################################
    ######################################################
    
    # Xj and X-j
    X_j <- X[, index]
    X_negj <- X[, -index]

    # regress X-j on xj, use least min lambda to estimate gamma(Step 4)
    lam_min = 0.1*max(abs(t(X_negj)%*%X_j)/(n*p))
    lam_max = 40*max(abs(t(X_negj)%*%X_j)/(n*p))
    cvfit <- glmnet::cv.glmnet(x = X_negj, y = X_j,alpha=1)
    # gamma <- coef(cvfit, s = cvfit$lambda.1se)[-1]
    gamma <- coef(cvfit, s = cvfit$lambda.min)[-1]
    z <- n^-0.5 * (X_j - X_negj %*% gamma)
    
    # V <- 1.25 * n^0.5 * norm(z, type = "2") / (t(z) %*% X_j)
    ######################################################
    ######################################################
    # eq. (12) for a point estimation of Bj(Step 6)
    dblasso[i] <- t(z)  %*% (Y - X_negj %*% betahat[-index[i]]) / (t(z) %*%  X[, index[i]])
    Variance <- (t(z) %*% z) / (t(z)  %*% X[, index[i]])^2
    ### calculate sigmahat
    if(1){
      fit <- glmnet::cv.glmnet(x = X, y = Y,alpha=1)
      betahat <- as.matrix((coef(fit, S = fit$lambda.min)[-1]))
      error <- (norm(Y - X %*% betahat, type = "2"))^2
      sigmahat <- (error / n)^0.5
    }
    stddev[i] <- sigmahat * Variance^.5
    quantile <- qnorm(1 - sig_level / 2, mean = 0, sd = 1, lower.tail = TRUE)
    lower[i] <- dblasso[i] - quantile * stddev[i]
    upper[i] <- dblasso[i] + quantile * stddev[i]
    cover0[i] <- ifelse(lower[i] * upper[i] > 0, 0, 1)
  }
  rs <- list(
    est_ddl = dblasso,
    se = stddev,
    sigmahat = sigmahat,
    len = upper - lower,
    CI=c(lower,upper),
    co0 = cover0
  )
  return(rs)
}

### DLassoAUC
