library(glmnet)
library(ncvreg)
### DLassoAUC
myDLassoAUC2ncvreg <- function(X, Y, index, K.factors, rho, rhop, sig_level = 0.05){
  output = farm.res(X,K.factors=K.factors,robust=FALSE) #default options
  U=output$factors
  XX=cbind(X,U)

  # resu=myDLasso0222(XX, Y, index, rho , rhop ,sig_level = 0.05)
  ####################################
  ####################################
  n <- dim(XX)[1]
  p <- dim(XX)[2]
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
    lam_min = 0.1*max(abs(t(XX)%*%Y)/(n*p))
    lam_max = 20*max(abs(t(XX)%*%Y)/(n*p))
    # the_lambdas = seq(lam_min,lam_max,(lam_max - lam_min)/10)
    the_lambdas = seq(lam_max,lam_min,(lam_min - lam_max)/10)

    fit0 <- glmnet::cv.glmnet(x = XX, y = Y,alpha=1, lambda=the_lambdas,penalty.factor = c(rep(1,ncol(X)),rep(0,ncol(U))))
    fit0 <- glmnet::cv.glmnet(x = XX, y = Y,alpha=1, penalty.factor = c(rep(1,ncol(X)),rep(0,ncol(U))))
    # fit0$lambda
    # fit0 <- glmnet::cv.glmnet(x = XX, y = Y,alpha=1,penalty.factor = c(rep(1,ncol(X)),rep(0,ncol(U))))
    # fit0$lambda
    # fit0$lambda.min
    fit = glmnet::glmnet(x=XX,y=Y,lambda=fit0$lambda.1se,penalty.factor = c(rep(1,ncol(X)),rep(0,ncol(U))))
    betahat <- as.matrix((coef(fit)[-1]))
  }
  if(1){
    fit0 = cv.ncvreg(XX,Y,penalty = "lasso",penalty.factor = c(rep(1,ncol(X)),rep(0,ncol(U))))
    fit0$lambda.min
    names(fit0)
    names(fit0$fit)
    betahat0 = fit0$fit$beta[,fit0$min]
    betahat = betahat0[-1]
    }
  

  ###################################################################
  ###################################################################
  ###################################################################
  ###################################################################
  i=1
  for (i in seq(length(index))) {
    X_negj <- XX[, -index[i]]
    # z <- myfind_z(P_X, index[i])
    ######################################################
    ######################################################
    # Xj and X-j
    X_j <- XX[, index]
    X_negj <- XX[, -index]

    if(0){
      # regress X-j on xj, use least min lambda to estimate gamma(Step 4)
      lam_min = 0.1*max(abs(t(X_negj)%*%X_j)/(n*p))
      lam_max = 20*max(abs(t(X_negj)%*%X_j)/(n*p))
      the_lambdas = seq(lam_max,lam_min,(lam_min - lam_max)/10)
      cvfit0 <- glmnet::cv.glmnet(x = X_negj, y = X_j,alpha=1,lambda=the_lambdas, penalty.factor = c(rep(1,(ncol(X_negj) - ncol(U))  ),rep(0,ncol(U))))
      # cvfit0 <- glmnet::cv.glmnet(x = X_negj, y = X_j,alpha=1, penalty.factor = c(rep(1,(ncol(X_negj) - ncol(U))  ),rep(0,ncol(U))))
      # gamma <- coef(cvfit, s = cvfit$lambda.1se)[-1]
      cvfit = glmnet::glmnet(x=X_negj,y=X_j,lambda=cvfit0$lambda.min, penalty.factor = c(rep(1,(ncol(X_negj) - ncol(U))  ),rep(0,ncol(U)))             )
      gamma <- coef(cvfit, s = cvfit$lambda.min)[-1]
      # z <- n^(-0.5) * (X_j - X_negj %*% gamma)
      z <- n^(-0.5) * (X_j - X_negj %*% gamma)
    }
    if(1){
        fit0 = cv.ncvreg(X_negj,X_j,penalty = "lasso",penalty.factor = c(rep(1,(ncol(X_negj) - ncol(U))  ),rep(0,ncol(U))))
        fit0$lambda.min
        names(fit0)
        names(fit0$fit)
        gamma = fit0$fit$beta[,fit0$min]
        z <- n^(-0.5) * (X_j - X_negj %*% gamma[-1])
    }

    dblasso[i] <- t(z)  %*% (Y - X_negj %*% betahat[-index[i]]) / (t(z) %*%  XX[, index[i]])
    Variance <- (t(z) %*% z) / (t(z)  %*% XX[, index[i]])^2
    ### sigmahat <- estimate_sigma1(XX, Y, rho)
    {
      fit <- cv.ncvreg(XX, Y,penalty = "lasso",penalty.factor = c(rep(1,(ncol(XX) - ncol(U))),rep(0,ncol(U))))
      
      betahat0 = fit$fit$beta[,fit$min]
      betahat = betahat0[-1]
      error <- (norm(Y - XX %*% betahat, type = "2"))^2
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