library(glmnet)

library(ncvreg)

myDLassoncvreg <- function(X, Y, index, rho , rhop ,sig_level = 0.05) {
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
  if(0){
    lam_min = 0.1*max(abs(t(X)%*%Y)/(n*p))
    lam_max = 10*max(abs(t(X)%*%Y)/(n*p))
    the_lambdas = seq(lam_min,lam_max,(lam_max - lam_min)/10)
    fit <- glmnet::cv.glmnet(x = X, y = Y,alpha=1,lambda=the_lambdas)
    betahat <- as.matrix((coef(fit, S = fit$lambda.min)[-1]))
  }

  if(1){
    fit0 = cv.ncvreg(X,Y,penalty = "lasso")
    fit0$lambda.min
    names(fit0)
    names(fit0$fit)
    betahat0 = fit0$fit$beta[,fit0$min]
    betahat = betahat0[-1]
  }
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
    if(0){
        lam_min = 0.1*max(abs(t(X_negj)%*%X_j)/(n*p))
        lam_max = 10*max(abs(t(X_negj)%*%X_j)/(n*p))
        cvfit <- glmnet::cv.glmnet(x = X_negj, y = X_j,alpha=1)
        # gamma <- coef(cvfit, s = cvfit$lambda.1se)[-1]
        gamma <- coef(cvfit, s = cvfit$lambda.min)[-1]
        z <- n^-0.5 * (X_j - X_negj %*% gamma)
    }
    if(1){
          fit0 = cv.ncvreg(X_negj,X_j,penalty = "lasso")
          # fit0$lambda.min
          # names(fit0)
          # names(fit0$fit)
          gamma0 = fit0$fit$beta[,fit0$min]
          gamma = gamma0[-1]
          z <- n^-0.5 * (X_j - X_negj %*% gamma)
    }
    # V <- 1.25 * n^0.5 * norm(z, type = "2") / (t(z) %*% X_j)
    ######################################################
    ######################################################
    # eq. (12) for a point estimation of Bj(Step 6)
    dblasso[i] <- t(z)  %*% (Y - X_negj %*% betahat[-index[i]]) / (t(z) %*%  X[, index[i]])
    Variance <- (t(z) %*% z) / (t(z)  %*% X[, index[i]])^2
    ## sigmahat
    {
      # UDV_list <- svd(X)
      # U <- UDV_list$u
      # D <- UDV_list$d
      # V <- t(UDV_list$v)
      
      # # forms QX and QY(STEP 1)
      # tau <- quantile(D, rho)
      # Dtilde <- pmin(D, tau)
      #Q <- diag(nrow(X)) - U %*% diag(1 - Dtilde / D) %*% t(U)
      # Q=diag(nrow(X))

      fit <- cv.ncvreg(X, Y,penalty = "lasso")
      betahat0 = fit$fit$beta[,fit$min]
      betahat = betahat0[-1]
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
