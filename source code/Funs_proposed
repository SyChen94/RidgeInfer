

## CV

CV_lam1<-function(X,Y,Lam,rhop2){
  N=nrow(X)
  res=rep(0,length(Lam))
  for (k in 1:length(Lam)) {
    err=0
    for (h in 1:5) {
      idx=((h-1)*N/5+1):(h*N/5)
      estt =mypro(X[-idx,], Y[-idx], index=1,rhop=rhop2,eta=Lam[k],sig_level = 0.05,alpha2=0.1)
      err=err+mean((Y[idx,]-X[idx,1] * estt$est_ddl)^2)
    }
    res[k]=err
  }
  idx=which.min(res)
  return(Lam[idx])
}
CV_lam<-function(X,Y,Lam,rhop2){
  N=nrow(X)
  res=matrix(0,nrow=100,ncol=(length(Lam)))
  for (k in 1:length(Lam)) {
    for (l in 1:100) {
      idx=sample(100,100)
      res[l,k] =mypro(X[idx,], Y[idx], index=1,rhop=rhop2,eta=Lam[k],sig_level = 0.05,alpha2=0.1)$est_ddl
    }
  }
  idx=which.min(apply(res,2,function(x) sum((x-signal_strength)^2/100)))
  return(Lam[idx])
}

## DDL
myridgesig <- function(y, x, alpha2=0.1,rhop1) {
  n1 <- dim(x)[1]
  p1 <- dim(x)[2]
  UDV_list = svd(x)
  U = UDV_list$u
  D = UDV_list$d
  V = t(UDV_list$v)
  #Reduce D, then make P trim
  Dtilde = pmin(D, quantile(D,rhop1))
  #P = diag(nrow(X_negj)) - U %*% diag(1 - Dtilde / D) %*% t(U)
  p = U %*% diag(Dtilde / D) %*% t(U)
  px = p %*% x
  py = p %*% y
  Eta <- alpha2 * max(t(px) %*% py) / (n1 * p1)
  Sn2 <- 1/n1*t(px)%*%px
  X2_re  = chol2inv(chol(px %*% t(px) / n1 + Eta * diag(nrow(px)))) 
  Sn_inv2= 1 / Eta * (diag(ncol(px)) - t(px) %*% X2_re %*%px  / n1)
  A <- px %*% Sn_inv2 %*% t(px) / n1
  #RES = py-px[,index[i]] * dblasso[i]
  RES = py
  sigmahat=1/n1*t(RES)%*%(diag(n1)-A)%*%RES/(1-1/n1*tr(A))
  return(sigmahat)
}
myridgesig1 <- function(y, x, alpha2=0.1,rhop1) {
  n1 <- dim(x)[1]
  p1 <- dim(x)[2]
  UDV_list = svd(x)
  U = UDV_list$u
  D = UDV_list$d
  V = t(UDV_list$v)
  #Reduce D, then make P trim
  Dtilde = pmin(D, quantile(D,rhop1))
  p = diag(n1)
  px = p %*% x
  py = p %*% y
  Eta <- alpha2 * max(t(px) %*% py) / (n1 * p1)
  Sn2 <- 1/n1*t(px)%*%px
  X2_re  = chol2inv(chol(px %*% t(px) / n1 + Eta * diag(nrow(px)))) 
  Sn_inv2= 1 / Eta * (diag(ncol(px)) - t(px) %*% X2_re %*%px  / n1)
  A <- px %*% Sn_inv2 %*% t(px) / n1
  #RES = py-px[,index[i]] * dblasso[i]
  RES = py
  sigmahat=1/n1*t(RES)%*%(diag(n1)-A)%*%RES/(1-1/n1*tr(A))
  return(sigmahat)
}

## proposed method
RidgeInfer = function(X,Y,index,rhop,sig_level = 0.05,alpha1=0.1,alpha2=0.1){
  #determines parameters
  n = dim(X)[1]
  p = dim(X)[2]
  dblasso = rep(NA,length(index))
  stddev = rep(NA,length(index))
  lower = rep(NA,length(index))
  upper = rep(NA,length(index))
  liein = rep(NA,length(index))
  bias = rep(NA,length(index))
  cover0=rep(NA,length(index))
  i=1
  for (i in seq(length(index))){
    X_negj=X[,-index[i]]
    #single value decomposition of X (Trim transform)
    UDV_list = svd(X_negj)
    U = UDV_list$u
    D = UDV_list$d
    V = t(UDV_list$v)
    #Reduce D, then make P trim
    Dtilde = pmin(D, quantile(D,rhop))
    #P = diag(nrow(X_negj)) - U %*% diag(1 - Dtilde / D) %*% t(U)
    P = diag(nrow(X_negj))
    P_X = P %*% X
    P_Y = P %*% Y
    PX_negj = P %*% X_negj
    eta <- alpha1 * max(t(P_X) %*% P_Y) / (n * p)
    ### a point estimation of Bj
    Sn2 <- 1/n*t(PX_negj)%*%PX_negj
    X2_re  = chol2inv(chol(PX_negj %*% t(PX_negj) / n + eta * diag(nrow(PX_negj)))) 
    Sn_inv2= 1 / eta * (diag(ncol(PX_negj)) - t(PX_negj) %*% X2_re %*%PX_negj  / n)
    A <- PX_negj %*% Sn_inv2 %*% t(PX_negj) / n
    
    dblasso[i]=t(P_X[,index[i]])%*%(diag(n)-A)%*%P_Y/(t(P_X[,index[i]])%*%(diag(n)-A)%*%P_X[,index[i]])
    dblasso[i]
    ## variance
    Variance=(t(P_X[,index[i]])%*%t((diag(n)-A))%*%(diag(n)-A)%*%P_X[,index[i]])/(t(P_X[,index[i]])%*%(diag(n)-A)%*%P_X[,index[i]])^2
    
    
    sigmahat=myridgesig1(Y,X,alpha2=alpha2,rhop1=rhop)
    
    stddev[i] = (sigmahat*Variance)^(0.5)
    
    #2-sided CI critical value for alpha
    qquantile = qnorm(1-sig_level/2,mean=0,sd = 1,lower.tail=TRUE)
    # #determine bounds of interval, from eq. (13) (Step 8)
    lower[i] = dblasso[i] - qquantile * stddev[i]
    upper[i] = dblasso[i] + qquantile * stddev[i]
    # bias[i] = dblasso[i] - true[i]
    # liein[i] = ifelse((true[i]>lower[i]) && (true[i]<upper[i]),1,0)
    cover0[i]= ifelse(lower[i]*upper[i]>0,0,1)
  }
  rs=list(est_ddl= dblasso,
       se = stddev[i],
       sigmahat=sigmahat,
      #  bi=bias,
       len=upper-lower,
       CI=c(lower,upper),
      #  li=liein,
       co0=cover0
  )
  return(rs)
}

#mypro(X,Y,index,rhop=0.5,eta=1e-5,sig_level = 0.05)
#mypro(X, Y, index, rhop = the_rhop, eta = the_eta,alpha2=0.1)
## proposed method
mypro_rawinv = function(X,Y,index,rhop=0.5,eta=1e-5,sig_level = 0.05){
  #determines parameters
  n = dim(X)[1]
  p = dim(X)[2]
  dblasso = rep(NA,length(index))
  stddev = rep(NA,length(index))
  lower = rep(NA,length(index))
  upper = rep(NA,length(index))
  liein = rep(NA,length(index))
  bias = rep(NA,length(index))
  cover0=rep(NA,length(index))
  i=1
  for (i in seq(length(index))){
    X_negj=X[,-index[i]]
    #single value decomposition of X (Trim transform)
    UDV_list = svd(X_negj)
    U = UDV_list$u
    D = UDV_list$d
    V = t(UDV_list$v)
    #Reduce D, then make P trim
    Dtilde = pmin(D, quantile(D,rhop))
    #P = diag(nrow(X_negj)) - U %*% diag(1 - Dtilde / D) %*% t(U)
    P = U %*% diag(Dtilde / D) %*% t(U)
    P_X = P %*% X
    P_Y = P %*% Y
    PX_negj = P %*% X_negj
    ### a point estimation of Bj
      Sn2 <- 1/n*t(PX_negj)%*%PX_negj
      X2_re  = chol2inv(chol(PX_negj %*% t(PX_negj) / n + eta * diag(nrow(PX_negj)))) 
      Sn_inv2= 1 / eta * (diag(ncol(PX_negj)) - t(PX_negj) %*% X2_re %*%PX_negj  / n)
      A <- PX_negj %*% Sn_inv2 %*% t(PX_negj) / n
    
    dblasso[i]=t(P_X[,index[i]])%*%(diag(n)-A)%*%P_Y/(t(P_X[,index[i]])%*%(diag(n)-A)%*%P_X[,index[i]])

    ## variance
    Variance=(t(P_X[,index[i]])%*%t((diag(n)-A))%*%(diag(n)-A)%*%P_X[,index[i]])/
              (t(P_X[,index[i]])%*%(diag(n)-A)%*%P_X[,index[i]])^2
    

    #RES = P_Y-P_X[,index[i]] * dblasso[i]
   # sigmahat=myridgesig(RES,PX_negj,eta_alpha=alpha2)
    sigmahat = 1
    
    
    stddev[i] = (sigmahat*Variance)^(0.5)
    
    #2-sided CI critical value for sig_level
    quantile = qnorm(1-sig_level/2,mean=0,sd = 1,lower.tail=TRUE)
    # #determine bounds of interval, from eq. (13) (Step 8)
    lower[i] = dblasso[i] - quantile * stddev[i]
    upper[i] = dblasso[i] + quantile * stddev[i]
    # bias[i] = dblasso[i] - true[i]
    # liein[i] = ifelse((true[i]>lower[i]) && (true[i]<upper[i]),1,0)
    cover0[i]= ifelse(lower[i]*upper[i]>0,0,1)
  }
  rs=list(est_ddl= dblasso,
       se = stddev[i],
       sigmahat=sigmahat,
      #  bi=bias,
       len=upper-lower,
       CI=c(lower,upper),
      #  li=liein,
       co0=cover0
  )
  return(rs)
}

## ridge without hidden

RidgeInfer_wh= function(X,Y,index,rhop,sig_level = 0.05,alpha1=0.1,alpha2=0.1){
  #determines parameters
  n = dim(X)[1]
  p = dim(X)[2]
  dblasso = rep(NA,length(index))
  stddev = rep(NA,length(index))
  lower = rep(NA,length(index))
  upper = rep(NA,length(index))
  liein = rep(NA,length(index))
  bias = rep(NA,length(index))
  cover0=rep(NA,length(index))
  i=1
  for (i in seq(length(index))){
    X_negj=X[,-index[i]]
    #single value decomposition of X (Trim transform)
    UDV_list = svd(X_negj)
    U = UDV_list$u
    D = UDV_list$d
    V = t(UDV_list$v)
    #Reduce D, then make P trim
    Dtilde = pmin(D, quantile(D,rhop))
    #P = diag(nrow(X_negj)) - U %*% diag(1 - Dtilde / D) %*% t(U)
    P = U %*% diag(Dtilde / D) %*% t(U)
    P_X = P %*% X
    P_Y = P %*% Y
    PX_negj = P %*% X_negj
    eta <- alpha1 * max(t(P_X) %*% P_Y) / (n * p)
    ### a point estimation of Bj
    Sn2 <- 1/n*t(PX_negj)%*%PX_negj
    X2_re  = chol2inv(chol(PX_negj %*% t(PX_negj) / n + eta * diag(nrow(PX_negj)))) 
    Sn_inv2= 1 / eta * (diag(ncol(PX_negj)) - t(PX_negj) %*% X2_re %*%PX_negj  / n)
    A <- PX_negj %*% Sn_inv2 %*% t(PX_negj) / n
    
    dblasso[i]=t(P_X[,index[i]])%*%(diag(n)-A)%*%P_Y/(t(P_X[,index[i]])%*%(diag(n)-A)%*%P_X[,index[i]])
    dblasso[i]
    ## variance
    Variance=(t(P_X[,index[i]])%*%t((diag(n)-A))%*%(diag(n)-A)%*%P_X[,index[i]])/(t(P_X[,index[i]])%*%(diag(n)-A)%*%P_X[,index[i]])^2
    
    
    sigmahat=myridgesig(Y,X,alpha2=alpha2,rhop1=rhop)
    
    stddev[i] = (sigmahat*Variance)^(0.5)
    
    #2-sided CI critical value for sig_level
    quantile = qnorm(1-sig_level/2,mean=0,sd = 1,lower.tail=TRUE)
    # #determine bounds of interval, from eq. (13) (Step 8)
    lower[i] = dblasso[i] - quantile * stddev[i]
    upper[i] = dblasso[i] + quantile * stddev[i]
    # bias[i] = dblasso[i] - true[i]
    # liein[i] = ifelse((true[i]>lower[i]) && (true[i]<upper[i]),1,0)
    cover0[i]= ifelse(lower[i]*upper[i]>0,0,1)
  }
  rs=list(est_ddl= dblasso,
          se = stddev[i],
          sigmahat=sigmahat,
          #  bi=bias,
          len=upper-lower,
          CI=c(lower,upper),
          #  li=liein,
          co0=cover0
  )
  return(rs)
}


