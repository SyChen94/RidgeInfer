generate_dataset = function(n,p, q=3, sigma=2, pert=1){
  #create H and Gamma with N(0,1) values and of appropriate size. H can be tuned with pert
  H = pert*matrix(rnorm(n*q,mean=0,sd=1),n,q,byrow = TRUE)
  Gamma = matrix(rnorm(q*p,mean=0,sd=1),q,p,byrow = TRUE)
  ###
  sigmaE=matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigmaE[i,j]=0.5^(abs(i-j))
    }
  }
  E = matrix(mvrnorm(n,rep(0,p),sigmaE),n,p,byrow = TRUE)
  
  #defined in eq. (2), high-dimensional measured covariates
  X = E + H %*% Gamma
  
  delta = matrix(rnorm(q*1,mean=0,sd=1),q,1,byrow = TRUE)
  
  #px1 matrix, creates beta with 1s in the first s entries and the remaining p-s as 0s
  beta = true
  
  #nx1 matrix with values of mean 0 and SD of sigma, error in Y independent of X
  nu = matrix(rnorm(n*1,mean=0,sd=sigma),n,1,byrow = TRUE)
  
  #eq. (1), the response of the Structural Equation Model
  Y = X %*% beta + H %*% delta + nu
  return_list = list("X"= X,"Y"= Y)
  return(return_list)
}
## liu's method
generate_dataset = function(n, p, q=3, sigma=2, pert=1){
  #create H and Gamma with N(0,1) values and of appropriate size. H can be tuned with pert
  H = pert*matrix(rnorm(n*q,mean=0,sd=1),n,q,byrow = TRUE)
  Gamma = matrix(rnorm(q*p,mean=0,sd=1),q,p,byrow = TRUE)
  ###
  T=
    rho=runif(n,0,1)
  sigmaE=matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      if(abs(j-i)<T){
        for (k in 1:(T-abs(j-i))) {
          sigmaE[i,j]=sigmaE[i,j]+rho[k]*rho[k+abs(j-i)]
        }
      }
    }
  }
  E = matrix(mvrnorm(n,rep(0,p),sigmaE),n,p,byrow = TRUE)
  #defined in eq. (2), high-dimensional measured covariates
  X = E + H %*% Gamma
  
  delta = matrix(rnorm(q*1,mean=0,sd=1),q,1,byrow = TRUE)
  
  #px1 matrix, creates beta with 1s in the first s entries and the remaining p-s as 0s
  beta = true
  
  #nx1 matrix with values of mean 0 and SD of sigma, error in Y independent of X
  nu = matrix(rnorm(n*1,mean=0,sd=sigma),n,1,byrow = TRUE)
  
  #eq. (1), the response of the Structural Equation Model
  Y = X %*% beta + H %*% delta + nu
  return_list = list("X"= X,"Y"= Y)
  return(return_list)
}
###等相关
generate_dataset = function(n, p, q=3, sigma=2, pert=1){
  #create H and Gamma with N(0,1) values and of appropriate size. H can be tuned with pert
  H = pert*matrix(rnorm(n*q,mean=0,sd=1),n,q,byrow = TRUE)
  Gamma = matrix(rnorm(q*p,mean=0,sd=1),q,p,byrow = TRUE)
  ###
  sigmaE=matrix(0.4,p,p)
  for (i in 1:p) {
    sigmaE[i,i]=1
  }
  E = matrix(mvrnorm(n,rep(0,p),sigmaE),n,p,byrow = TRUE)
  
  #defined in eq. (2), high-dimensional measured covariates
  X = E + H %*% Gamma
  
  delta = matrix(rnorm(q*1,mean=0,sd=1),q,1,byrow = TRUE)
  
  #px1 matrix, creates beta with 1s in the first s entries and the remaining p-s as 0s
  beta = true
  
  #nx1 matrix with values of mean 0 and SD of sigma, error in Y independent of X
  nu = matrix(rnorm(n,mean=0,sd=sigma),n,1,byrow = TRUE)
  
  #eq. (1), the response of the Structural Equation Model
  Y = X %*% beta + H %*% delta + nu
  return_list = list("X"= X,"Y"= Y)
  return(return_list)
}
