## Golbal Args
#setwd("C:/Users/ELSAID/Dropbox/Elsayed_Mun/Second Project/Real Data/Numerical Study/Numerical study_cross validation/Realdata_HKP")
rm(list=ls())
M0 <- 2 # number of mixture component
p0 <- 2 # number of predictor
## BMPARMC and BMPBUTTO with corr = 0.81
maxit_em0 <- 200
epsilon_em0 <- 1e-5 
k_fold <- 5
CN0 <- 2000        ##{500,2000}
Nsim <- 1000         ##{200,500,2000}
Nsample0 <- c(60,100)       # sample size in estimation 
n0t <- Nsample0[1]

## Reading data

## BMPARMC and BMPBUTTO with corr = 0.81

load("BMPARMC_BMPBUTTO.RData")
x1 <- BMPARMC_BMPBUTTO$X1
x2 <- BMPARMC_BMPBUTTO$X2
Y0 <- BMPARMC_BMPBUTTO$Y

Beta01 <-  c(0.004,0.006)
Beta02 <-  c(0.011,0.005)
Beta0 <- rbind(Beta01,Beta02)
pi0 <- c(0.6,0.4)
sigma0 <- c(0.0081,0.0144)

##*********************************************

## Starting point as a little bit different with true parameter

Beta.01 <- Beta0+0.001
Beta.02 <- Beta0+0.001
pi.01 <- c(0.7,0.3)
pi.02 <- c(0.7,0.3)
sigma.01 <- sigma0+0.001
sigma.02 <- sigma0+0.001

##*********************************************

# Data Matrix

Design_Matrix <- cbind(Y0,x1,x2)

##*********************************************


##*********************************************

## Generates New design matrix data

## Args:
## n: sample size
## p: number of predictors
## rho: correlation between first 2 predictors
## phi: correlation between last 2 predictors

## Returns:
## Matrix of size (n*p)

# Design_Matrix <- function(n,p,rho,phi)
# {
#   x <- matrix(NA, nrow=n, ncol=(p-1))
#   wr <- rnorm(n*(p-1),mean=0,sd=1)
#   wr <- matrix(wr, nrow=n, ncol=p-1)
#   wf <- rnorm(n,0,1)
#   for (i in 1:n) {
#     for (j in 1:(p-1)) {
#       if (j==1||j==2) {
#         x[i,j]= (1-phi^2)*wr[i,j]+phi*wf[i]
#       }
#       else {
#         x[i,j]= (1-rho^2)*wr[i,j]+rho*wf[i]
#       }
#     }
#   }
#   x <- cbind(1,x)
#   return(x)
# }

#x <- Design_Matrix(n=n0t,p=p0,rho=rho0t,phi=phi0t)
# cor(X[,-1])

##*********************************************

## Generate vector of response variable  

## Args:
## x: Data Matrix (n by p)
## beta: true coefficients (M by p)
## pi: mixing proportion (of length M)
## M: number of mixture component
## sigma: vector of variances of length M

## Returns:
## y: a matrix of responses and memberships of size (n*(M+1))
## the first column is y values and (2:M+1) are memberships

# Design_Response_Mixture <- function(x, beta,pi,M,sigma)
# {
#   n <- dim(x)[1]
#   L <- length(pi)
#   m <- sample(x=1:M, size=n, replace=TRUE, prob=pi)
#   y <- rep(NA,n)
#   class_Assigned=matrix(data=0, nrow=n, ncol=L)
#   for(i in 1:n)
#   {
#     if (m[i]==1)
#     {
#       y[i]<-x[i,]%*%beta[1,]+rnorm(1,0,sqrt(sigma[1]))
#     }
#     
#     else y[i]<-x[i,]%*%beta[2,]+rnorm(1,0,sqrt(sigma[2]))
#     
#     class_Assigned[i,m[i]]<-1
#   }
#   res<-cbind(y,class_Assigned)
#   return(res)
# }

# x <- Design_Matrix(n=n0t,p=p0,rho=rho0t,phi=phi0t)
# Yres <- Design_Response_Mixture(x=x,beta=Beta0,pi=pi0,M=M0,sigma=sigma0)
# y=Yres[,1]

##***************************************************

Design_Response_prediction <- function(x, beta,pi,M)
{
  n <- dim(x)[1]
  L <- length(pi)
  m <- sample(x=1:M, size=n, replace=TRUE, prob=pi)
  y <- rep(NA,n)
  class_Assigned=matrix(data=0, nrow=n, ncol=L)
  for(i in 1:n)
  {
    if (m[i]==1)
    {
      y[i]<-x[i,]%*%beta[1,]
    }
    
    else y[i]<-x[i,]%*%beta[2,]
    
    class_Assigned[i,m[i]]<-1
  }
  res<-cbind(y,class_Assigned)
  return(res)
}

##***************************************************
## Generate E_Step
## Generate tau
## Args:
## x: Data Matrix
## y: Response variable
## beta: coefficients (updated from previous iteration) (M by p)
## Pi: Probability vector according to number of mixture component (Pi)
## sigma: variance vector according to number of mixture component (Pi)

## Returns:
## tau: a Matix of size (n*k) whose each row corresponds probability vector for each instance

Design_Tau <- function(x,y,beta,pi,sigma)
{
  n <- dim(x)[1]
  L <- length(pi)
  tau_Matrix <- matrix(data=NA, nrow=n, ncol=L)
  for (j in 1:L) 
  {
    for (i in 1:n) 
    {
      P_j <- dnorm(y[i],x[i,]%*%beta[j,],sqrt(sigma[j]))
      tau_ij <- pi[j]*P_j
      tau_Matrix[i,j] <- tau_ij
    }
  }
  dev <- rowSums(tau_Matrix)
  tau <- apply(tau_Matrix, 2, "/", dev)
  return(tau)
}

#tau <- Design_Tau(x,y,beta=Beta0,pi=pi0,sigma=sigma0)
# apply(tau,1,sum)

##*********************************************

Z_EM <- function(tau)
{
  n <- dim(tau)[1]
  L <- dim(tau)[2]
  z <- matrix(data=0, nrow=n, ncol=L)
  for (i in 1:n)
  {
    z[i,] <- t(rmultinom(n=1,size=1,prob=c(tau[i,])))
  }
  return(z)
}

#Z <- Z_EM(tau=tau)
# apply(Z,2,sum)

##*********************************************

# second step for E_step (Calculation Z_i)

Genration_Z_i <- function(tau)
{
  n <- dim(tau)[1]
  L <- dim(tau)[2]
  z <- matrix(data=0, nrow=n, ncol=L)
  for (i in 1:n)
  {
    ind=which(tau[i,]==max(tau[i,]))
    z[i,ind] <- 1
  }
  return(z)
}

#Z <- Genration_Z_i(tau=tau)
# apply(Z,2,sum)

##*********************************************

## Liklihood function

## Args:
## x: Data Matrix
## Y: Response variable
## beta: coefficients (updated from previous iteration) (M by p)
## Pi: Probability vector according to number of mixture component (Pi)
## sigma: vector of variances according to number of mixture component

## Returns:
## Log_Liklihood result 

Log_Liklihood <- function(x,y,pi,beta,sigma)
{
  n <- dim(x)[1]
  L <- length(pi)
  LL_Matrix <- matrix(data=NA, nrow=n, ncol=L)
  for (j in 1:L) 
  {
    for (i in 1:n) 
    {
      P_j <- dnorm(y[i],x[i,]%*%beta[j,],sqrt(sigma[j]))
      Mat_ij <- pi[j]*P_j
      LL_Matrix[i,j] <- Mat_ij
    }
  }
  Sum_Rows <- rowSums(LL_Matrix)
  Log_Sum_Rows <- sum(log(Sum_Rows))
  return(Log_Sum_Rows)
}

#Log_Liklihood(x,y,pi=pi0,beta=Beta0,sigma=sigma0)

##*********************************************

## calculate Sigma for EM
## Args:
## tau : probability vector for j Regression
## X: A numeric data matrix
## y: Response vector
## beta: vetor of parameter according to j th regression


Calc_Sigma<- function(tau,x,y,beta)
{
  n <- dim(x)[1]
  part_1<- c()
  
  for(i in 1:n)
  {
    part_1[i]=tau[i]*(y[i]-x[i,]%*%beta)^2
    
  }
  S=sum(part_1)/sum(tau)
  return(S)
}

#S<- Calc_Sigma(tau[,1],x,y,Beta0[1,])

##*********************************************

## Computes the K_hkb

## Args:
## p: number of predictors
## sig2ls : the estimate of sigma^2 by lS 
## Bhat: the LSE of Beta parameters

## Returns:
## an integer K_HKB 

K.hkb <- function(p,sig2ls,Bhat)
{
  k=c()
  n=length(sig2ls)
  for (i in 1:n)
  {
    k[i] <- c(p*sig2ls[i]/(t(Bhat[i,])%*%Bhat[i,]))
    
  }
  return(k)
}

##**********************************************

## the estimate of dOpt
## Args:
## sig2hat: The final estimate variance of ridge
## k: The first [arameter of liu function from  K.hkb 
## alpha: obtained from canolical
## lamda: lamda from canolical

.dOpt <- function(sig2hat,k,alpha,lamda)
{
  la <- diag(lamda)
  u <- 0
  l <- 0
  p <- length(alpha)
  for (i in 1:p) 
  {
    u <- u + ((la[i]*(sig2hat-(k*(alpha[i]^2))))/((la[i]+k)^3))
    l <- l + ((la[i]*(la[i]*alpha[i]^2+sig2hat))/((la[i]+k)^4))
  }
  d <- u/l
  return(d)
}


# ##*********************************************
# ## calculate canolical form for ridge to calculate d for liu
# canolical <- function(x,y,w,k,sig2hat)
# {
#   np <- dim(x)
#   n  <- np[1]
#   p  <- np[2]
#   W <- as.numeric(w)
#   XtX <- t(x)%*%((diag(W)))%*%x
#   e <- eigen(XtX)
#   eigenvalue <- e$values
#   eigenvector <- e$vectors
#   z <- x%*%eigenvector
#   lamda <- t(z)%*%((diag(W)))%*%z
#   ll <- diag(lamda)
#   laKI <- lamda+(k*diag(nrow(lamda)))
#   Zty  <- t(z)%*%((diag(W)))%*%y
#   alpha <- solve(laKI)%*%Zty
#   d <- .dOpt(sig2hat,k,alpha,lamda)
#   return (d)
# }

##*********************************************
## calculate canolical form for ridge to calculate d for liu
canolical_Liu <- function(x,y,w,k,sig2hat,bR)
{
  np <- dim(x)
  n  <- np[1]
  p  <- np[2]
  W <- as.numeric(w)
  XtX <- t(x)%*%((diag(W)))%*%x
  e <- eigen(XtX)
  eigenvalue <- e$values
  eigenvector <- e$vectors
  z <- x%*%eigenvector
  lamda <- t(z)%*%((diag(W)))%*%z
  ll <- diag(lamda)
  laKI <- lamda+(k*diag(nrow(lamda)))
  Zty  <- t(z)%*%((diag(W)))%*%y
#  alpha <- solve(laKI)%*%Zty
  alpha <- t(eigenvector)%*%bR
  d <- .dOpt(sig2hat,k,alpha,lamda)
  return (d)
}

##*********************************************

## calculate the updated Beta for EM
## Args:
## X: A numeric data matrix
## y: Response vector
## tau : probability vector for j Regression
## beta.k : one previous iteration of beta
## beta.old : two previous iteration of beta

update_Beta <- function(x,y,tau,beta.k,beta.old)
{
  b_old2 <- beta.old
  W <- as.numeric(tau)
  XtX <- crossprod(x,diag(W)%*%x)
  detV <- det(XtX)
  if(is.na(detV)) return(b_old2)
  if(detV==0) return(b_old2)
  eXtX <- eigen(XtX)$values
  CN <- sqrt(abs(max(eXtX)/min(eXtX)))
  if(is.nan(CN)) return(c(b_old2))
  Xty <- crossprod(x,W*y)
  if(CN >= CN0)  return(c(b_old2))     ## singularity
  else          beta.k <- solve(XtX,Xty)
  if(any(is.nan(beta.k))) return(c(b_old2))
  return(c(beta.k))
}

##*********************************************

update_Beta_Ridge <- function(x,y,tau,beta.k,beta.old,La)
{
  pp <- ncol(x)
  IPP <- diag(pp)
  b_old2 <- beta.old
  W <- as.numeric(tau)
  XtX <- crossprod(x,diag(W)%*%x)
  V=XtX+La*IPP
  detV <- det(V)
  if(is.na(detV)) return(b_old2)
  if(detV==0) return(b_old2)
  ev <- eigen(V)$values
  CN <- sqrt(abs(max(ev)/min(ev)))
  if(is.nan(CN)) return(c(b_old2))
  Xty <- crossprod(x,W*y)
  if(CN >= CN0)  return(c(b_old2))     ## singularity
  else          beta.k <- solve(V,Xty)
  if(any(is.nan(beta.k))) return(c(b_old2))
  return(c(beta.k))
}


##*********************************************

update_Beta_Liu <- function(x,y,tau,beta.k,beta.old,La,dFix)
{
  pp <- ncol(x)
  IPP <- diag(pp)
  b_old2 <- beta.old
  W <- as.numeric(tau)
  XtX <- crossprod(x,diag(W)%*%x)
  V=XtX+La*IPP
  detV <- det(V)
  if(is.na(detV)) return(b_old2)
  if(detV==0) return(b_old2)
  ev <- eigen(V)$values
  CN <- sqrt(abs(max(ev)/min(ev)))
  if(is.nan(CN)) return(c(b_old2))
  if(CN >= CN0)  return(c(b_old2))     ## singularit
  Xty <- crossprod(x,W*y)
  B_Hat <- solve(V,Xty)
  if(any(is.nan(B_Hat))) return(b_old2)
  V1=XtX-(dFix*IPP)
  beta.k=solve(V,V1%*%B_Hat)
  if(any(is.nan(beta.k))) return(b_old2)
  return(c(beta.k))
}

##*********************************************

## Estimator : EM
## Args:
## epsilon: Numeric tolarence parameter regarding EM algorithm (Stopping rule).
## X: A numeric data matrix
## y: Response vector
## beta.0 and pi.0 and sigma.0 intial values not true parameters
## maxit_EM: Integer maximum number of iterations regarding EM algorithm.

## Returns:
## a list of two elements
## $pi.new: mixing proportion of length M
## $beta.new: Matrix of predictors of size (M by p)
## $sigma: vector of variance of length M 


em_ml <-function(epsilon,x,y,beta.0,beta.old,pi.0,pi.old,sigma.0,sigma.old,maxit_EM,printout=F)
{
  pp <- dim(x)[2]
  M <- length(pi.0)  
  Beta_Mle <- matrix(data=NA, nrow=M, ncol=pp)
  convergence<-F
  beta.k <- beta.0
  pi.k  <- pi.0
  sigma.k <- sigma.0
  sigma.MLe <-c()
  for(k in 1:maxit_EM)
  {
    #print(k)
    if(any(is.na(beta.k)))
    {
      #MSE <- .msepi_one_iteration(pi.old[1],pi0[1])
      #MSE_Beta <- .msebeta_one_iteration(beta.old,Beta0)
      #pi_vec <- c(pi.old,MSE,MSE_Beta,k,101)
      #ret <- list(k=k,pi=round(pi_vec,3),beta=beta.old)
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(pi.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(sigma.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    tau <- Design_Tau(x,y,beta.k,pi.k,sigma.k) 
    if(any(is.na(tau)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    Z <- Z_EM(tau=tau)
    
    if(any(is.na(Z)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    for (j in 1:ncol(Z)) 
    {
      beta_mle <- update_Beta(x,y,tau[,j],beta.k[j,],beta.old[j,])
      #beta_mle <- solve(t(x)%*%diag(tau[,j])%*%x)%*%(t(x)%*%diag(tau[,j])%*%y)
      Beta_Mle[j,] <- beta_mle
      sigma.MLe[j]<- Calc_Sigma(tau[,j],x,y,Beta_Mle[j,])
    }
    pi.new <- colMeans(tau)
    beta.new <- Beta_Mle
    sigma.new <- sigma.MLe
    
    #if(printout==T) print(cbind(pi.new,beta.new),digits=3)
    
    changeValuell <- abs(Log_Liklihood(x,y,pi.new,beta.new,sigma.new)-Log_Liklihood(x,y,pi.k,beta.k,sigma.k)) 
    #print(changeValuell)
  
    if(changeValuell < epsilon)    break
    pi.old <- pi.k
    pi.k <- pi.new 
    beta.old <- beta.k
    beta.k <- beta.new
    sigma.old <- sigma.k
    sigma.k <- sigma.new
  } 
  
  ret <- ret <- list(k=k,pi=pi.new,sigma=sigma.new,beta=beta.new)
  return(ret)
}

# EM_Result<- em_ml(epsilon=epsilon_em0,x,y,beta.0=Beta.01,beta.old=Beta.02,
#                   pi.0=pi.01,pi.old=pi.02,sigma.0=sigma.01
#                   ,sigma.old=sigma.02,maxit_EM=maxit_em0)
# 
# EM_Result$k
# EM_Result$pi
# EM_Result$sigma
# EM_Result$beta

##**********************************************************


em_ridge <-function(epsilon,x,y,beta.0,beta.old,pi.0,pi.old,sigma.0,sigma.old,maxit_EM,k_ridge,printout=F)
{
  pp <- dim(x)[2]
  M <- length(pi.0)  
  Beta_Ridge <- matrix(data=NA, nrow=M, ncol=pp)
  convergence<-F
  beta.k <- beta.0
  pi.k  <- pi.0
  sigma.k <- sigma.0
  sigma.Ridge <-c()
  for(k in 1:maxit_EM)
  {
    #print(k)
    if(any(is.na(beta.k)))
    {
      #MSE <- .msepi_one_iteration(pi.old[1],pi0[1])
      #MSE_Beta <- .msebeta_one_iteration(beta.old,Beta0)
      #pi_vec <- c(pi.old,MSE,MSE_Beta,k,101)
      #ret <- list(k=k,pi=round(pi_vec,3),beta=beta.old)
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(pi.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(sigma.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    tau <- Design_Tau(x,y,beta.k,pi.k,sigma.k) 
    if(any(is.na(tau)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    Z <- Z_EM(tau=tau)
    
    if(any(is.na(Z)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    for (j in 1:ncol(Z)) 
    {
      beta_ridge <- update_Beta_Ridge(x,y,tau[,j],beta.k[j,],beta.old[j,],k_ridge[j])
      #beta_mle <- solve(t(x)%*%diag(tau[,j])%*%x)%*%(t(x)%*%diag(tau[,j])%*%y)
      Beta_Ridge[j,] <- beta_ridge
      sigma.Ridge[j]<- Calc_Sigma(tau[,j],x,y,Beta_Ridge[j,])
    }
    pi.new <- colMeans(tau)
    beta.new <- Beta_Ridge
    sigma.new <- sigma.Ridge
    
    #if(printout==T) print(cbind(pi.new,beta.new),digits=3)
    
    changeValuell <- abs(Log_Liklihood(x,y,pi.new,beta.new,sigma.new)-Log_Liklihood(x,y,pi.k,beta.k,sigma.k)) 
    #print(changeValuell)
    
    if(changeValuell < epsilon)    break
    pi.old <- pi.k
    pi.k <- pi.new 
    beta.old <- beta.k
    beta.k <- beta.new
    sigma.old <- sigma.k
    sigma.k <- sigma.new
    
  } 
  
  ret <- ret <- list(k=k,pi=pi.new,sigma=sigma.new,beta=beta.new)
  return(ret)
}

# k_ridge <- K.hkb(p0,EM_Result$sigma,EM_Result$beta)
# 
# EM_Result_ridge<- em_ridge(epsilon=epsilon_em0,x,y,beta.0=Beta.01,beta.old=Beta.02,
#                    pi.0=pi.01,pi.old=pi.02,sigma.0=sigma.01
#                    ,sigma.old=sigma.02,maxit_EM=maxit_em0,k_ridge)
# 
# EM_Result_ridge$k
# EM_Result_ridge$pi
# EM_Result_ridge$sigma
# EM_Result_ridge$beta

##***********************************************************

##**********************************************************


em_Liu <-function(epsilon,x,y,beta.0,beta.old,pi.0,pi.old,sigma.0,sigma.old,maxit_EM,k_Liu,Sigma_Ridge_EM,printout=F)
{
  pp <- dim(x)[2]
  M <- length(pi.0)  
  Beta_Liu <- matrix(data=NA, nrow=M, ncol=pp)
  convergence<-F
  beta.k <- beta.0
  pi.k  <- pi.0
  sigma.k <- sigma.0
  sigma.Liu <-c()
  for(k in 1:maxit_EM)
  {
    #print(k)
    if(any(is.na(beta.k)))
    {
      #MSE <- .msepi_one_iteration(pi.old[1],pi0[1])
      #MSE_Beta <- .msebeta_one_iteration(beta.old,Beta0)
      #pi_vec <- c(pi.old,MSE,MSE_Beta,k,101)
      #ret <- list(k=k,pi=round(pi_vec,3),beta=beta.old)
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(pi.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(sigma.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    tau <- Design_Tau(x,y,beta.k,pi.k,sigma.k) 
    if(any(is.na(tau)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    Z <- Z_EM(tau=tau)
    
    if(any(is.na(Z)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    for (j in 1:ncol(Z)) 
    {
      BR  <-   update_Beta_Ridge(x,y,tau[,j],beta.k[j,],beta.old[j,],k_Liu[j])
      d_j <- canolical_Liu(x,y,tau[,j],k_Liu[j],Sigma_Ridge_EM[j],BR)
      beta_Liu <- update_Beta_Liu(x,y,tau[,j],beta.k[j,],beta.old[j,],k_Liu[j],d_j)
      #beta_mle <- solve(t(x)%*%diag(tau[,j])%*%x)%*%(t(x)%*%diag(tau[,j])%*%y)
      Beta_Liu[j,] <- beta_Liu
      sigma.Liu[j]<- Calc_Sigma(tau[,j],x,y,Beta_Liu[j,])
    }
    pi.new <- colMeans(tau)
    beta.new <- Beta_Liu
    sigma.new <- sigma.Liu
    
    #if(printout==T) print(cbind(pi.new,beta.new),digits=3)
    
    changeValuell <- abs(Log_Liklihood(x,y,pi.new,beta.new,sigma.new)-Log_Liklihood(x,y,pi.k,beta.k,sigma.k)) 
    #print(changeValuell)
    
    if(changeValuell < epsilon)    break
    pi.old <- pi.k
    pi.k <- pi.new 
    beta.old <- beta.k
    beta.k <- beta.new
    sigma.old <- sigma.k
    sigma.k <- sigma.new
  } 
  
  ret <- ret <- list(k=k,pi=pi.new,sigma=sigma.new,beta=beta.new)
  return(ret)
}

# k_Liu <- K.hkb(p0,EM_Result_ridge$sigma,EM_Result_ridge$beta)
# 
# EM_Result_Liu<- em_Liu(epsilon=epsilon_em0,x,y,beta.0=Beta.01,beta.old=Beta.02,
#                    pi.0=pi.01,pi.old=pi.02,sigma.0=sigma.01
#                    ,sigma.old=sigma.02,maxit_EM=maxit_em0,k_Liu,EM_Result_ridge$sigma)
# 
# EM_Result_Liu$k
# EM_Result_Liu$pi
# EM_Result_Liu$sigma
# EM_Result_Liu$beta

##***********************************************************



## Estimator : CEM for ML method

CEM_ml <-function(epsilon,x,y,beta.0,beta.old,pi.0,pi.old,sigma.0,sigma.old,maxit_EM,printout=F)
{
  pp <- dim(x)[2]
  M <- length(pi.0)  
  Beta_Mle <- matrix(data=NA, nrow=M, ncol=pp)
  convergence<-F
  beta.k <- beta.0
  pi.k  <- pi.0
  sigma.k <- sigma.0
  sigma.MLe <-c()
  for(k in 1:maxit_EM)
  {
    #print(k)
    if(any(is.na(beta.k)))
    {
      #MSE <- .msepi_one_iteration(pi.old[1],pi0[1])
      #MSE_Beta <- .msebeta_one_iteration(beta.old,Beta0)
      #pi_vec <- c(pi.old,MSE,MSE_Beta,k,101)
      #ret <- list(k=k,pi=round(pi_vec,3),beta=beta.old)
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(pi.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(sigma.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    tau <- Design_Tau(x,y,beta.k,pi.k,sigma.k) 
    if(any(is.na(tau)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    Z <- Genration_Z_i(tau=tau)
    
    if(any(is.na(Z)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    for (j in 1:ncol(Z)) 
    {
      zInd <- Z[,j]==1
      if(sum(zInd)==1||sum(zInd)==0)
      {
        ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
        return(ret)
      }
      beta_mle <- update_Beta(x=x[zInd,],y=y[zInd],tau=tau[zInd,j],beta.k[j,],beta.old[j,])
      #beta_mle <- solve(t(x)%*%diag(tau[,j])%*%x)%*%(t(x)%*%diag(tau[,j])%*%y)
      Beta_Mle[j,] <- beta_mle
      sigma.MLe[j]<- Calc_Sigma(tau=tau[zInd,j],x=x[zInd,],y=y[zInd],Beta_Mle[j,])
    }
    pi.new <- colMeans(Z)
    beta.new <- Beta_Mle
    sigma.new <- sigma.MLe
    
    #if(printout==T) print(cbind(pi.new,beta.new),digits=3)
    
    changeValuell <- abs(Log_Liklihood(x,y,pi.new,beta.new,sigma.new)-Log_Liklihood(x,y,pi.k,beta.k,sigma.k)) 
    #print(changeValuell)
    
    if(changeValuell < epsilon)    break
    pi.old <- pi.k
    pi.k <- pi.new 
    beta.old <- beta.k
    beta.k <- beta.new
    sigma.old <- sigma.k
    sigma.k <- sigma.new
  } 
  
  ret <- ret <- list(k=k,pi=pi.new,sigma=sigma.new,beta=beta.new)
  return(ret)
}

# CEM_Result_ml<- CEM_ml(epsilon=epsilon_em0,x,y,beta.0=Beta.01,beta.old=Beta.02,
#                   pi.0=pi.01,pi.old=pi.02,sigma.0=sigma.01
#                   ,sigma.old=sigma.02,maxit_EM=maxit_em0)
# 
# CEM_Result_ml$k
# CEM_Result_ml$pi
# CEM_Result_ml$sigma
# CEM_Result_ml$beta

##*********************************************

CEM_ridge <-function(epsilon,x,y,beta.0,beta.old,pi.0,pi.old,sigma.0,sigma.old,maxit_EM,k_ridge,printout=F)
{
  pp <- dim(x)[2]
  M <- length(pi.0)  
  Beta_Ridge <- matrix(data=NA, nrow=M, ncol=pp)
  convergence<-F
  beta.k <- beta.0
  pi.k  <- pi.0
  sigma.k <- sigma.0
  sigma.Ridge <-c()
  for(k in 1:maxit_EM)
  {
    #print(k)
    if(any(is.na(beta.k)))
    {
      #MSE <- .msepi_one_iteration(pi.old[1],pi0[1])
      #MSE_Beta <- .msebeta_one_iteration(beta.old,Beta0)
      #pi_vec <- c(pi.old,MSE,MSE_Beta,k,101)
      #ret <- list(k=k,pi=round(pi_vec,3),beta=beta.old)
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(pi.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(sigma.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    tau <- Design_Tau(x,y,beta.k,pi.k,sigma.k) 
    if(any(is.na(tau)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    Z <- Genration_Z_i(tau=tau)
    
    if(any(is.na(Z)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    for (j in 1:ncol(Z)) 
    {
      zInd <- Z[,j]==1
      if(sum(zInd)==1||sum(zInd)==0)
      {
        ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
        return(ret)
      }
    
      beta_ridge <- update_Beta_Ridge(x=x[zInd,],y=y[zInd],tau=tau[zInd,j],beta.k[j,],beta.old[j,],k_ridge[j])
      Beta_Ridge[j,] <- beta_ridge
      sigma.Ridge[j]<- Calc_Sigma(tau=tau[zInd,j],x=x[zInd,],y=y[zInd],Beta_Ridge[j,])
    }
    pi.new <- colMeans(Z)
    beta.new <- Beta_Ridge
    sigma.new <- sigma.Ridge
    
    #if(printout==T) print(cbind(pi.new,beta.new),digits=3)
    
    changeValuell <- abs(Log_Liklihood(x,y,pi.new,beta.new,sigma.new)-Log_Liklihood(x,y,pi.k,beta.k,sigma.k)) 
    #print(changeValuell)
    
    if(changeValuell < epsilon)    break
    pi.old <- pi.k
    pi.k <- pi.new 
    beta.old <- beta.k
    beta.k <- beta.new
    sigma.old <- sigma.k
    sigma.k <- sigma.new
  } 
  
  ret <- ret <- list(k=k,pi=pi.new,sigma=sigma.new,beta=beta.new)
  return(ret)
}

# k_ridge_CEM <- K.hkb(p0,CEM_Result_ml$sigma,CEM_Result_ml$beta)
# 
# CEM_Result_ridge<- CEM_ridge(epsilon=epsilon_em0,x,y,beta.0=Beta.01,beta.old=Beta.02,
#                    pi.0=pi.01,pi.old=pi.02,sigma.0=sigma.01
#                    ,sigma.old=sigma.02,maxit_EM=maxit_em0,k_ridge_CEM)
# 
# CEM_Result_ridge$k
# CEM_Result_ridge$pi
# CEM_Result_ridge$sigma
# CEM_Result_ridge$beta

##***********************************************************

CEM_Liu <-function(epsilon,x,y,beta.0,beta.old,pi.0,pi.old,sigma.0,sigma.old,maxit_EM,k_Liu,Sigma_Ridge_CEM,printout=F)
{
  pp <- dim(x)[2]
  M <- length(pi.0)  
  Beta_Liu <- matrix(data=NA, nrow=M, ncol=pp)
  convergence<-F
  beta.k <- beta.0
  pi.k  <- pi.0
  sigma.k <- sigma.0
  sigma.Liu <-c()
  for(k in 1:maxit_EM)
  {
    #print(k)
    if(any(is.na(beta.k)))
    {
      #MSE <- .msepi_one_iteration(pi.old[1],pi0[1])
      #MSE_Beta <- .msebeta_one_iteration(beta.old,Beta0)
      #pi_vec <- c(pi.old,MSE,MSE_Beta,k,101)
      #ret <- list(k=k,pi=round(pi_vec,3),beta=beta.old)
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(pi.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(sigma.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    tau <- Design_Tau(x,y,beta.k,pi.k,sigma.k) 
    if(any(is.na(tau)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    Z <- Genration_Z_i(tau=tau)
    
    if(any(is.na(Z)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    for (j in 1:ncol(Z)) 
    {
      zInd <- Z[,j]==1
      if(sum(zInd)==1||sum(zInd)==0)
      {
        ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
        return(ret)
      }
      BR <- update_Beta_Ridge(x=x[zInd,],y=y[zInd],tau=tau[zInd,j],beta.k[j,],beta.old[j,],k_Liu[j])
      d_j <- canolical_Liu(x=x[zInd,],y=y[zInd],tau[zInd,j],k_Liu[j],Sigma_Ridge_CEM[j],BR)
      beta_Liu <- update_Beta_Liu(x=x[zInd,],y=y[zInd],tau=tau[zInd,j],beta.k[j,],beta.old[j,],k_Liu[j],d_j)
      Beta_Liu[j,] <- beta_Liu
      sigma.Liu[j]<- Calc_Sigma(tau=tau[zInd,j],x=x[zInd,],y=y[zInd],Beta_Liu[j,])
    }
    pi.new <- colMeans(Z)
    beta.new <- Beta_Liu
    sigma.new <- sigma.Liu
    
    #if(printout==T) print(cbind(pi.new,beta.new),digits=3)
    
    changeValuell <- abs(Log_Liklihood(x,y,pi.new,beta.new,sigma.new)-Log_Liklihood(x,y,pi.k,beta.k,sigma.k)) 
    #print(changeValuell)
    
    if(changeValuell < epsilon)    break
    pi.old <- pi.k
    pi.k <- pi.new 
    beta.old <- beta.k
    beta.k <- beta.new
    sigma.old <- sigma.k
    sigma.k <- sigma.new
  } 
  
  ret <- ret <- list(k=k,pi=pi.new,sigma=sigma.new,beta=beta.new)
  return(ret)
}



##*******************************************************************

## Estimator : SEM for ML method

SEM_ml <-function(epsilon,x,y,beta.0,beta.old,pi.0,pi.old,sigma.0,sigma.old,maxit_EM,printout=F)
{
  pp <- dim(x)[2]
  M <- length(pi.0)  
  Beta_Mle <- matrix(data=NA, nrow=M, ncol=pp)
  convergence<-F
  beta.k <- beta.0
  pi.k  <- pi.0
  sigma.k <- sigma.0
  sigma.MLe <-c()
  for(k in 1:maxit_EM)
  {
    #print(k)
    if(any(is.na(beta.k)))
    {
      #MSE <- .msepi_one_iteration(pi.old[1],pi0[1])
      #MSE_Beta <- .msebeta_one_iteration(beta.old,Beta0)
      #pi_vec <- c(pi.old,MSE,MSE_Beta,k,101)
      #ret <- list(k=k,pi=round(pi_vec,3),beta=beta.old)
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(pi.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(sigma.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    tau <- Design_Tau(x,y,beta.k,pi.k,sigma.k) 
    if(any(is.na(tau)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    Z <- Z_EM(tau=tau)
    
    if(any(is.na(Z)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    for (j in 1:ncol(Z)) 
    {
      zInd <- Z[,j]==1
      if(sum(zInd)==1||sum(zInd)==0)
      {
        ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
        return(ret)
      }
      beta_mle <- update_Beta(x=x[zInd,],y=y[zInd],tau=tau[zInd,j],beta.k[j,],beta.old[j,])
      #beta_mle <- solve(t(x)%*%diag(tau[,j])%*%x)%*%(t(x)%*%diag(tau[,j])%*%y)
      Beta_Mle[j,] <- beta_mle
      sigma.MLe[j]<- Calc_Sigma(tau=tau[zInd,j],x=x[zInd,],y=y[zInd],Beta_Mle[j,])
    }
    pi.new <- colMeans(Z)
    beta.new <- Beta_Mle
    sigma.new <- sigma.MLe
    
    #if(printout==T) print(cbind(pi.new,beta.new),digits=3)
    
    changeValuell <- abs(Log_Liklihood(x,y,pi.new,beta.new,sigma.new)-Log_Liklihood(x,y,pi.k,beta.k,sigma.k)) 
    #print(changeValuell)
    
    if(changeValuell < epsilon)    break
    pi.old <- pi.k
    pi.k <- pi.new 
    beta.old <- beta.k
    beta.k <- beta.new
    sigma.old <- sigma.k
    sigma.k <- sigma.new
  } 
  
  ret <- ret <- list(k=k,pi=pi.new,sigma=sigma.new,beta=beta.new)
  return(ret)
}

# SEM_Result_ml<- SEM_ml(epsilon=epsilon_em0,x,y,beta.0=Beta.01,beta.old=Beta.02,
#                   pi.0=pi.01,pi.old=pi.02,sigma.0=sigma.01
#                   ,sigma.old=sigma.02,maxit_EM=maxit_em0)
# 
# SEM_Result_ml$k
# SEM_Result_ml$pi
# SEM_Result_ml$sigma
# SEM_Result_ml$beta

##**********************************************************************

SEM_ridge <-function(epsilon,x,y,beta.0,beta.old,pi.0,pi.old,sigma.0,sigma.old,maxit_EM,k_ridge,printout=F)
{
  pp <- dim(x)[2]
  M <- length(pi.0)  
  Beta_Ridge <- matrix(data=NA, nrow=M, ncol=pp)
  convergence<-F
  beta.k <- beta.0
  pi.k  <- pi.0
  sigma.k <- sigma.0
  sigma.Ridge <-c()
  for(k in 1:maxit_EM)
  {
    #print(k)
    if(any(is.na(beta.k)))
    {
      #MSE <- .msepi_one_iteration(pi.old[1],pi0[1])
      #MSE_Beta <- .msebeta_one_iteration(beta.old,Beta0)
      #pi_vec <- c(pi.old,MSE,MSE_Beta,k,101)
      #ret <- list(k=k,pi=round(pi_vec,3),beta=beta.old)
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(pi.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(sigma.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    tau <- Design_Tau(x,y,beta.k,pi.k,sigma.k) 
    if(any(is.na(tau)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    Z <- Z_EM(tau=tau)
    
    if(any(is.na(Z)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    for (j in 1:ncol(Z)) 
    {
      zInd <- Z[,j]==1
      if(sum(zInd)==1||sum(zInd)==0)
      {
        ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
        return(ret)
      }
      
      beta_ridge <- update_Beta_Ridge(x=x[zInd,],y=y[zInd],tau=tau[zInd,j],beta.k[j,],beta.old[j,],k_ridge[j])
      Beta_Ridge[j,] <- beta_ridge
      sigma.Ridge[j]<- Calc_Sigma(tau=tau[zInd,j],x=x[zInd,],y=y[zInd],Beta_Ridge[j,])
    }
    pi.new <- colMeans(Z)
    beta.new <- Beta_Ridge
    sigma.new <- sigma.Ridge
    
    #if(printout==T) print(cbind(pi.new,beta.new),digits=3)
    
    changeValuell <- abs(Log_Liklihood(x,y,pi.new,beta.new,sigma.new)-Log_Liklihood(x,y,pi.k,beta.k,sigma.k)) 
    #print(changeValuell)
    
    if(changeValuell < epsilon)    break
    pi.old <- pi.k
    pi.k <- pi.new 
    beta.old <- beta.k
    beta.k <- beta.new
    sigma.old <- sigma.k
    sigma.k <- sigma.new
  } 
  
  ret <- ret <- list(k=k,pi=pi.new,sigma=sigma.new,beta=beta.new)
  return(ret)
}

# k_ridge_SEM <- K.hkb(p0,SEM_Result_ml$sigma,SEM_Result_ml$beta)
# 
# SEM_Result_ridge<- CEM_ridge(epsilon=epsilon_em0,x,y,beta.0=Beta.01,beta.old=Beta.02,
#                    pi.0=pi.01,pi.old=pi.02,sigma.0=sigma.01
#                    ,sigma.old=sigma.02,maxit_EM=maxit_em0,k_ridge_SEM)
# 
# SEM_Result_ridge$k
# SEM_Result_ridge$pi
# SEM_Result_ridge$sigma
# SEM_Result_ridge$beta

##**********************************************************************


SEM_Liu <-function(epsilon,x,y,beta.0,beta.old,pi.0,pi.old,sigma.0,sigma.old,maxit_EM,k_Liu,Sigma_Ridge_SEM,printout=F)
{
  pp <- dim(x)[2]
  M <- length(pi.0)  
  Beta_Liu <- matrix(data=NA, nrow=M, ncol=pp)
  convergence<-F
  beta.k <- beta.0
  pi.k  <- pi.0
  sigma.k <- sigma.0
  sigma.Liu <-c()
  for(k in 1:maxit_EM)
  {
    #print(k)
    if(any(is.na(beta.k)))
    {
      #MSE <- .msepi_one_iteration(pi.old[1],pi0[1])
      #MSE_Beta <- .msebeta_one_iteration(beta.old,Beta0)
      #pi_vec <- c(pi.old,MSE,MSE_Beta,k,101)
      #ret <- list(k=k,pi=round(pi_vec,3),beta=beta.old)
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(pi.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    if(any(is.na(sigma.k)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    tau <- Design_Tau(x,y,beta.k,pi.k,sigma.k) 
    if(any(is.na(tau)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    Z <- Z_EM(tau=tau)
    
    if(any(is.na(Z)))
    {
      ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
      return(ret)
    }
    
    for (j in 1:ncol(Z)) 
    {
      zInd <- Z[,j]==1
      if(sum(zInd)==1||sum(zInd)==0)
      {
        ret <- list(k=k,pi=pi.old,sigma=sigma.old,beta=beta.old)
        return(ret)
      }
      BR <- update_Beta_Ridge(x=x[zInd,],y=y[zInd],tau=tau[zInd,j],beta.k[j,],beta.old[j,],k_Liu[j])
      d_j <- canolical_Liu(x=x[zInd,],y=y[zInd],tau[zInd,j],k_Liu[j],Sigma_Ridge_SEM[j],BR)
      beta_Liu <- update_Beta_Liu(x=x[zInd,],y=y[zInd],tau=tau[zInd,j],beta.k[j,],beta.old[j,],k_Liu[j],d_j)
      Beta_Liu[j,] <- beta_Liu
      sigma.Liu[j]<- Calc_Sigma(tau=tau[zInd,j],x=x[zInd,],y=y[zInd],Beta_Liu[j,])
    }
    pi.new <- colMeans(Z)
    beta.new <- Beta_Liu
    sigma.new <- sigma.Liu
    
    #if(printout==T) print(cbind(pi.new,beta.new),digits=3)
    
    changeValuell <- abs(Log_Liklihood(x,y,pi.new,beta.new,sigma.new)-Log_Liklihood(x,y,pi.k,beta.k,sigma.k)) 
    #print(changeValuell)
    
    if(changeValuell < epsilon)    break
    pi.old <- pi.k
    pi.k <- pi.new 
    beta.old <- beta.k
    beta.k <- beta.new
    sigma.old <- sigma.k
    sigma.k <- sigma.new
  } 
  
  ret <- ret <- list(k=k,pi=pi.new,sigma=sigma.new,beta=beta.new)
  return(ret)
}



##**********************************************************************
## Computes the MSE of Beta hat (over Nsim) 

## Args:
## Mhat: Matrix of beta hat values

## Returns:
## data frame containing Mean and Median and 5% and 95% confidence interval

.mseBeta <- function(Mhat,Btrue)
{
  #dif1 <- sqrt((Mhat-Btrue)^2)
  jpn <- dim(Mhat)
  j <- jpn[1]
  p <- jpn[2]
  n <- jpn[3]
  dif1 <- ((Mhat-Btrue)^2)
  dif2 <- apply(dif1, 3, sum)
  dif3 <- sqrt(dif2/(j*p))
  meanBeta <- round(mean(dif3),4)
  quanBeta <- round((quantile(dif3, probs=c(0.025,0.5,0.975))),4)
  quanlen <- quanBeta[3]-quanBeta[1]
  tab <- c(meanBeta,quanBeta[2],quanBeta[1],quanBeta[3],quanlen)
  return(tab)
}

##*********************************************
## Computes the MSE of pi (over Nsim) 

## Args:
## P_hat: Matrix of beta hat values (Nsim*M0)
## P_true: true pi (1*M0)

## Returns:
## data frame containing Mean and Median and 5% and 95% confidence interval

.msepi <- function(P_hat,P_true)
{
  #dif1 <- sqrt((P_hat-P_true)^2)
  dif1 <- sqrt((P_hat-P_true)^2)
  mean_pi <- round(mean(dif1),4)
  quan_pi <- round((quantile(dif1, probs=c(0.025,0.5,0.975))),4)
  quanlen <- quan_pi[3]-quan_pi[1]
  tab <- c(mean_pi,quan_pi[2],quan_pi[1],quan_pi[3],quanlen)
  return(tab)
}

##*********************************************

## Computes the MSE of sigma (over Nsim) 

## Args:
## sigma_hat: Matrix of sigma values (Nsim*M0)
## sigma_true: true sigma (1*M0)

## Returns:
## data frame containing Mean and Median and 5% and 95% confidence interval

.msesigma <- function(sigma_hat,sigma_true)
{
  #dif1 <- sqrt((sigma_hat-sigma_true)^2)
  dif1 <- ((t(sigma_hat)-sigma_true)^2)
  dif2 <- colSums(dif1)
  dif3 <- sqrt(dif2/2)
  mean_sigma <- round(mean(dif3),4)
  quan_sigma <- round((quantile(dif3, probs=c(0.025,0.5,0.975))),4)
  quanlen <- quan_sigma[3]-quan_sigma[1]
  tab <- c(mean_sigma,quan_sigma[2],quan_sigma[1],quan_sigma[3],quanlen)
  return(tab)
}

##*********************************************


## Computes the bias of Beta hat (over Nsim) 

## Args:
## Mhat: Matrix of beta hat values

## Returns:
## data frame containing Mean and Median and 5% and 95% confidence interval

.biasBeta <- function(Mhat,Btrue)
{
  dif1 <- Mhat-Btrue
  dif2 <- apply(dif1, 3, sum)
  meanBeta <- round(mean(dif2),3)
  #quanBeta <- round(quantile(dif2, probs=c(0.025,0.5,0.975)),3)
  #quanlen <- quanBeta[3]-quanBeta[1]
  #tab <- c(meanBeta,quanBeta[2],quanBeta[1],quanBeta[3],quanlen)
  return(meanBeta)
}

##*********************************************

## Computes the bias of pi (over Nsim) 

## Args:
## P_hat: Matrix of beta hat values (Nsim*M0)
## P_true: true pi (1*M0)

## Returns:
## data frame containing Mean and Median and 5% and 95% confidence interval

.biaspi <- function(P_hat,P_true)
{
  #dif1 <- sqrt((P_hat-P_true)^2)
  dif1 <- P_hat-P_true
  mean_pi <- round(mean(dif1),3)
  #quan_pi <- round(quantile(dif1, probs=c(0.025,0.5,0.975)),3)
  #quanlen <- quan_pi[3]-quan_pi[1]
  #tab <- c(mean_pi,quan_pi[2],quan_pi[1],quan_pi[3],quanlen)
  return(mean_pi)
}

##*********************************************

## Computes the bias of sigma (over Nsim) 

## Args:
## sigma_hat: Matrix of sigma values (Nsim*M0)
## sigma_true: true sigma (1*M0)

## Returns:
## data frame containing Mean and Median and 5% and 95% confidence interval

.biassigma <- function(sigma_hat,sigma_true)
{
  dif1 <- t(sigma_hat)-sigma_true
  dif2 <- colSums(dif1)
  mean_sigma <- round(mean(dif2),3)
  #quan_pi <- round(quantile(dif1, probs=c(0.025,0.5,0.975)),3)
  #quanlen <- quan_pi[3]-quan_pi[1]
  #tab <- c(mean_pi,quan_pi[2],quan_pi[1],quan_pi[3],quanlen)
  return(mean_sigma)
}

##*********************************************


## simulator

Simulation <- function(p,epsilon,M,Beta,beta.0,beta.old,pi,pi.0,pi.old,sigma,sigma.0,sigma.old,maxit_EM,printout=F)
{
  ## intialization for prediction (p+1) means the number of indep variables + response variable
  
  Data <- array(NA,dim=c(n0t,p+1,Nsim))
  
  ## intialization for ML Method
  resArray_LS_EM <-resArray_LS_CEM <- resArray_LS_SEM <- array(NA,dim=c(2,p,Nsim))
  MO_LS_EM_pi <- MO_LS_EM_sigma <- matrix(NA, nrow = Nsim, ncol = length(pi))
  MO_LS_CEM_pi <- MO_LS_CEM_sigma <- matrix(NA, nrow = Nsim, ncol = length(pi))
  MO_LS_SEM_pi <- MO_LS_SEM_sigma <- matrix(NA, nrow = Nsim, ncol = length(pi))
  k_LS_EM <- k_LS_CEM <- k_LS_SEM <- c()
  
  ## intialization for Ridge Method
  
  resArray_ridge_EM <-resArray_ridge_CEM <- resArray_ridge_SEM <- array(NA,dim=c(2,p,Nsim))
  MO_Ridge_EM_pi <- MO_Ridge_EM_sigma <- matrix(NA, nrow = Nsim, ncol = length(pi))
  MO_Ridge_CEM_pi <- MO_Ridge_CEM_sigma <- matrix(NA, nrow = Nsim, ncol = length(pi))
  MO_Ridge_SEM_pi <- MO_Ridge_SEM_sigma <- matrix(NA, nrow = Nsim, ncol = length(pi))
  k_Ridge_EM <- k_Ridge_CEM <- k_Ridge_SEM <- c()
  
  ## intialization for Liu Method
  
  resArray_Liu_EM <-resArray_Liu_CEM <- resArray_Liu_SEM <- array(NA,dim=c(2,p,Nsim))
  MO_Liu_EM_pi <- MO_Liu_EM_sigma <- matrix(NA, nrow = Nsim, ncol = length(pi))
  MO_Liu_CEM_pi <- MO_Liu_CEM_sigma <- matrix(NA, nrow = Nsim, ncol = length(pi))
  MO_Liu_SEM_pi <- MO_Liu_SEM_sigma <- matrix(NA, nrow = Nsim, ncol = length(pi))
  k_Liu_EM <- k_Liu_CEM <- k_Liu_SEM <- c()
  
  #MO_LS_pi <- MO_Ridge_pi <-  MO_Liu_pi <- matrix(NA, nrow = Nsim, ncol = length(pi))
  for (j in 1:Nsim)
  {
    cat("j = ",j,'\n')
    ## Estimation  Procedure 
    Ind  <- sample(nrow(Design_Matrix),size = n0t,replace = T)
    Data_train <- Design_Matrix[Ind,]
    
    x <- Data_train[,-1]
    y <- Data_train[,1]
    
    ## combine data with response variables
    
    #data <- cbind(y,x)
    Data[,,j] <- Data_train  
    
    ## Methods of EM
    
    ## ML Method
    
    EM_Result_LS <- em_ml(epsilon=epsilon,x=x,y=y,beta.0=beta.0,
                          beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                          sigma.old=sigma.old,maxit_EM=maxit_EM,printout=F)
    
    K_LS_EM <- EM_Result_LS$k
    P_ml_EM <- EM_Result_LS$pi
    Sigma_LS_EM <- EM_Result_LS$sigma
    B_ml_EM <- EM_Result_LS$beta
    
    resArray_LS_EM[,,j] <- B_ml_EM
    MO_LS_EM_pi[j,] <- P_ml_EM
    MO_LS_EM_sigma[j,] <-Sigma_LS_EM
    k_LS_EM[j] <-K_LS_EM
    
    
    ## Ridge Method
    
    k_ridge <- K.hkb(p,Sigma_LS_EM,B_ml_EM)
    
    EM_Result_ridge <- em_ridge(epsilon=epsilon,x=x,y=y,beta.0=beta.0,
                          beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                          sigma.old=sigma.old,maxit_EM=maxit_EM,k_ridge,printout=F)
    
    K_Ridge_EM <- EM_Result_ridge$k
    P_Ridge_EM <- EM_Result_ridge$pi
    Sigma_Ridge_EM <- EM_Result_ridge$sigma
    B_Ridge_EM <- EM_Result_ridge$beta
    
    resArray_ridge_EM[,,j] <- B_Ridge_EM
    MO_Ridge_EM_pi[j,] <- P_Ridge_EM
    MO_Ridge_EM_sigma[j,] <-Sigma_Ridge_EM
    k_Ridge_EM[j] <-K_Ridge_EM
    
    ## Liu Method
    
    k_Liu <- K.hkb(p,Sigma_Ridge_EM,B_Ridge_EM )

    EM_Result_Liu <- em_Liu(epsilon=epsilon,x=x,y=y,beta.0=beta.0,
                                beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                                sigma.old=sigma.old,maxit_EM=maxit_EM,k_Liu,Sigma_Ridge_EM,printout=F)

    K_Liu_EM <- EM_Result_Liu$k
    P_Liu_EM <- EM_Result_Liu$pi
    Sigma_Liu_EM <- EM_Result_Liu$sigma
    B_Liu_EM <- EM_Result_Liu$beta

    resArray_Liu_EM[,,j] <- B_Liu_EM
    MO_Liu_EM_pi[j,] <- P_Liu_EM
    MO_Liu_EM_sigma[j,] <-Sigma_Liu_EM
    k_Liu_EM[j] <-K_Liu_EM
    
    
    ## Methods of CEM
    
    ## ML Method
    
    CEM_Result_LS <- CEM_ml(epsilon=epsilon,x=x,y=y,beta.0=beta.0,
                          beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                          sigma.old=sigma.old,maxit_EM=maxit_EM,printout=F)
    
    K_ml_CEM <- CEM_Result_LS$k
    P_ml_CEM <- CEM_Result_LS$pi
    Sigma_ml_CEM <- CEM_Result_LS$sigma
    B_ml_CEM <- CEM_Result_LS$beta
    
    resArray_LS_CEM[,,j] <- B_ml_CEM
    MO_LS_CEM_pi[j,] <- P_ml_CEM
    MO_LS_CEM_sigma[j,] <-Sigma_ml_CEM
    k_LS_CEM[j] <-K_ml_CEM
    
    ## Ridge Method
    
    k_ridge <- K.hkb(p,Sigma_ml_CEM,B_ml_CEM)
    
    CEM_Result_ridge <- CEM_ridge(epsilon=epsilon,x=x,y=y,beta.0=beta.0,
                                beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                                sigma.old=sigma.old,maxit_EM=maxit_EM,k_ridge,printout=F)
    
    K_Ridge_CEM <- CEM_Result_ridge$k
    P_Ridge_CEM <- CEM_Result_ridge$pi
    Sigma_Ridge_CEM <- CEM_Result_ridge$sigma
    B_Ridge_CEM <- CEM_Result_ridge$beta
    
    resArray_ridge_CEM[,,j] <- B_Ridge_CEM
    MO_Ridge_CEM_pi[j,] <- P_Ridge_CEM
    MO_Ridge_CEM_sigma[j,] <-Sigma_Ridge_CEM
    k_Ridge_CEM[j] <-K_Ridge_CEM
    
    
    ## Liu Method
    
    k_Liu <- K.hkb(p,Sigma_Ridge_CEM,B_Ridge_CEM )
    
    CEM_Result_Liu <- CEM_Liu(epsilon=epsilon,x=x,y=y,beta.0=beta.0,
                            beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                            sigma.old=sigma.old,maxit_EM=maxit_EM,k_Liu,Sigma_Ridge_CEM,printout=F)
    
    K_Liu_CEM <- CEM_Result_Liu$k
    P_Liu_CEM <- CEM_Result_Liu$pi
    Sigma_Liu_CEM <- CEM_Result_Liu$sigma
    B_Liu_CEM <- CEM_Result_Liu$beta
    
    resArray_Liu_CEM[,,j] <- B_Liu_CEM
    MO_Liu_CEM_pi[j,] <- P_Liu_CEM
    MO_Liu_CEM_sigma[j,] <-Sigma_Liu_CEM
    k_Liu_CEM[j] <-K_Liu_CEM
    
    
    ## Methods of SEM
    
    ## ML Method
    
    SEM_Result_LS <- SEM_ml(epsilon=epsilon,x=x,y=y,beta.0=beta.0,
                            beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                            sigma.old=sigma.old,maxit_EM=maxit_EM,printout=F)
    
    K_ml_SEM <- SEM_Result_LS$k
    P_ml_SEM <- SEM_Result_LS$pi
    Sigma_ml_SEM <- SEM_Result_LS$sigma
    B_ml_SEM <- SEM_Result_LS$beta
    
    resArray_LS_SEM[,,j] <- B_ml_SEM
    MO_LS_SEM_pi[j,] <- P_ml_SEM
    MO_LS_SEM_sigma[j,] <-Sigma_ml_SEM
    k_LS_SEM[j] <-K_ml_SEM
    
    ## Ridge Method
    
    k_ridge <- K.hkb(p,Sigma_ml_SEM,B_ml_SEM)
    
    SEM_Result_ridge <- SEM_ridge(epsilon=epsilon,x=x,y=y,beta.0=beta.0,
                                  beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                                  sigma.old=sigma.old,maxit_EM=maxit_EM,k_ridge,printout=F)
    
    K_Ridge_SEM <- SEM_Result_ridge$k
    P_Ridge_SEM <- SEM_Result_ridge$pi
    Sigma_Ridge_SEM <- SEM_Result_ridge$sigma
    B_Ridge_SEM <- SEM_Result_ridge$beta
    
    resArray_ridge_SEM[,,j] <- B_Ridge_SEM
    MO_Ridge_SEM_pi[j,] <- P_Ridge_SEM
    MO_Ridge_SEM_sigma[j,] <-Sigma_Ridge_SEM
    k_Ridge_SEM[j] <-K_Ridge_SEM
    
    ## Liu Method
    
    k_Liu <- K.hkb(p,Sigma_Ridge_SEM,B_Ridge_SEM )
    
    SEM_Result_Liu <- SEM_Liu(epsilon=epsilon,x=x,y=y,beta.0=beta.0,
                              beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                              sigma.old=sigma.old,maxit_EM=maxit_EM,k_Liu,Sigma_Ridge_SEM,printout=F)
    
    K_Liu_SEM <- SEM_Result_Liu$k
    P_Liu_SEM <- SEM_Result_Liu$pi
    Sigma_Liu_SEM <- SEM_Result_Liu$sigma
    B_Liu_SEM <- SEM_Result_Liu$beta
    
    resArray_Liu_SEM[,,j] <- B_Liu_SEM
    MO_Liu_SEM_pi[j,] <- P_Liu_SEM
    MO_Liu_SEM_sigma[j,] <-Sigma_Liu_SEM
    k_Liu_SEM[j] <-K_Liu_SEM
    
    
  }
  
  return(list(Beta_LS_EM=resArray_LS_EM,Pi_LS_EM=MO_LS_EM_pi,
              Sigma_LS_EM=MO_LS_EM_sigma,K_LS_EM=k_LS_EM,
              Beta_LS_CEM=resArray_LS_CEM,Pi_LS_CEM=MO_LS_CEM_pi,
              Sigma_LS_CEM=MO_LS_CEM_sigma,K_LS_CEM=k_LS_CEM,
              Beta_LS_SEM=resArray_LS_SEM,Pi_LS_SEM=MO_LS_SEM_pi,
              Sigma_LS_SEM=MO_LS_SEM_sigma,K_LS_SEM=k_LS_SEM,
              
              Beta_Ridge_EM=resArray_ridge_EM,Pi_Ridge_EM=MO_Ridge_EM_pi,
              Sigma_Ridge_EM=MO_Ridge_EM_sigma,K_Ridge_EM=k_Ridge_EM,
              
              Beta_Ridge_CEM=resArray_ridge_CEM,Pi_Ridge_CEM=MO_Ridge_CEM_pi,
              Sigma_Ridge_CEM=MO_Ridge_CEM_sigma,K_Ridge_CEM=k_Ridge_CEM,
              
              Beta_Ridge_SEM=resArray_ridge_SEM,Pi_Ridge_SEM=MO_Ridge_SEM_pi,
              Sigma_Ridge_SEM=MO_Ridge_SEM_sigma,K_Ridge_SEM=k_Ridge_SEM,
              
              Beta_Liu_EM=resArray_Liu_EM,Pi_Liu_EM=MO_Liu_EM_pi,
              Sigma_Liu_EM=MO_Liu_EM_sigma,K_Liu_EM=k_Liu_EM,
              
              Beta_Liu_CEM=resArray_Liu_CEM,Pi_Liu_CEM=MO_Liu_CEM_pi,
              Sigma_Liu_CEM=MO_Liu_CEM_sigma,K_Liu_CEM=k_Liu_CEM,
              
              Beta_Liu_SEM=resArray_Liu_SEM,Pi_Liu_SEM=MO_Liu_SEM_pi,
              Sigma_Liu_SEM=MO_Liu_SEM_sigma,K_Liu_SEM=k_Liu_SEM,Data=Data))
}


## Results

Sim_Result=Simulation(p=p0,epsilon=epsilon_em0,M=M0,
                      Beta=Beta0,beta.0=Beta.01,beta.old=Beta.02,
                      pi=pi0,pi.0=pi.01,pi.old=pi.02,sigma=sigma0,sigma.0=sigma.01
                      ,sigma.old=sigma.02,maxit_EM=maxit_em0)

## LS results for 3 methods
Beta_LS_EM <-Sim_Result$Beta_LS_EM
Pi_LS_EM <-Sim_Result$Pi_LS_EM
Sigma_LS_EM <- Sim_Result$Sigma_LS_EM
K_LS_EM <- Sim_Result$K_LS_EM

Beta_LS_CEM <-Sim_Result$Beta_LS_CEM
Pi_LS_CEM <-Sim_Result$Pi_LS_CEM
Sigma_LS_CEM <- Sim_Result$Sigma_LS_CEM
K_LS_CEM <- Sim_Result$K_LS_CEM

Beta_LS_SEM <-Sim_Result$Beta_LS_SEM
Pi_LS_SEM <-Sim_Result$Pi_LS_SEM
Sigma_LS_SEM <- Sim_Result$Sigma_LS_SEM
K_LS_SEM <- Sim_Result$K_LS_SEM

## Ridge results for 3 methods

Beta_Ridge_EM <-Sim_Result$Beta_Ridge_EM
Pi_Ridge_EM <-Sim_Result$Pi_Ridge_EM
Sigma_Ridge_EM <- Sim_Result$Sigma_Ridge_EM
K_Ridge_EM <- Sim_Result$K_Ridge_EM

Beta_Ridge_CEM <-Sim_Result$Beta_Ridge_CEM
Pi_Ridge_CEM <-Sim_Result$Pi_Ridge_CEM
Sigma_Ridge_CEM <- Sim_Result$Sigma_Ridge_CEM
K_Ridge_CEM <- Sim_Result$K_Ridge_CEM

Beta_Ridge_SEM <-Sim_Result$Beta_Ridge_SEM
Pi_Ridge_SEM <-Sim_Result$Pi_Ridge_SEM
Sigma_Ridge_SEM <- Sim_Result$Sigma_Ridge_SEM
K_Ridge_SEM <- Sim_Result$K_Ridge_SEM

## Liu results for 3 methods

Beta_Liu_EM <-Sim_Result$Beta_Liu_EM
Pi_Liu_EM <-Sim_Result$Pi_Liu_EM
Sigma_Liu_EM <- Sim_Result$Sigma_Liu_EM
K_Liu_EM <- Sim_Result$K_Liu_EM

Beta_Liu_CEM <-Sim_Result$Beta_Liu_CEM
Pi_Liu_CEM <-Sim_Result$Pi_Liu_CEM
Sigma_Liu_CEM <- Sim_Result$Sigma_Liu_CEM
K_Liu_CEM <- Sim_Result$K_Liu_CEM

Beta_Liu_SEM <-Sim_Result$Beta_Liu_SEM
Pi_Liu_SEM <-Sim_Result$Pi_Liu_SEM
Sigma_Liu_SEM <- Sim_Result$Sigma_Liu_SEM
K_Liu_SEM <- Sim_Result$K_Liu_SEM

Data <- Sim_Result$Data

##*********************************************


##*********************************************

# Wrapping
## Args:
## Beta_LS_EM: Matrix of beta hat values for LS with dimension (2*3*Nsim)
## Pi_LS_EM: Matrix of Pi-hat for LS with dimension (Nsim*M0)
## Sigma_LS_EM: Matrix of sigma-hat values for Ls with dimension (Nsim*M0)
## Beta_LS_CEM: Matrix of beta hat values for LS with dimension (2*3*Nsim) (EM algorithm)
## Pi_LS_CEM: Matrix of Pi-hat for LS with dimension (Nsim*M0) (CEML algorithm)
## Sigma_LS_CEM: Matrix of sigma-hat values for Ls with dimension (Nsim*M0) (CEM algorithm)
## Beta_LS_SEM: Matrix of beta hat values for LS with dimension (2*3*Nsim) (SEM algorithm)
## Pi_LS_SEM: Matrix of Pi-hat for LS with dimension (Nsim*M0) (SEM algorithm)
## Sigma_LS_SEM: Matrix of sigma-hat values for Ls with dimension (Nsim*M0) (SEM algorithm)
## B_true: Matrix of true values with dimension (2*3)
## Pi_true: vector of true values of pi with dimension (1*2)
## sigma: vector of true values of sigma with dimension (1*2)

## Returns:
## data frame for MSE for each Beta and pi

Wrapping <- function(Beta_LS_EM,Pi_LS_EM,Sigma_LS_EM,
                     Beta_LS_CEM,Pi_LS_CEM,Sigma_LS_CEM,
                     Beta_LS_SEM,Pi_LS_SEM,Sigma_LS_SEM,
                     Beta_Ridge_EM,Pi_Ridge_EM,Sigma_Ridge_EM,
                     Beta_Ridge_CEM,Pi_Ridge_CEM,Sigma_Ridge_CEM,
                     Beta_Ridge_SEM,Pi_Ridge_SEM,Sigma_Ridge_SEM,
                     Beta_Liu_EM,Pi_Liu_EM,Sigma_Liu_EM,
                     Beta_Liu_CEM,Pi_Liu_CEM,Sigma_Liu_CEM,
                     Beta_Liu_SEM,Pi_Liu_SEM,Sigma_Liu_SEM
                     ,B_true,Pi_true,sigma_true)
{
  #naIndLS <- which(!is.na(Beta_LS_EM[1,1,]))
  #Beta_LS_EM <- Beta_LS_EM[,,naIndLS]
  
  Barray <- array(B_true,dim=c(M0,p0,Nsim))
  
  ## MSE calculation EM
  MSE_LS_Beta_EM   <- .mseBeta(Mhat =  Beta_LS_EM,Btrue = Barray)
  MSE_LS_Pi_EM     <- .msepi(P_hat =  Pi_LS_EM[,1],P_true = Pi_true[1])
  MSE_LS_sigma_EM  <- .msesigma(sigma_hat =  Sigma_LS_EM,sigma_true = sigma_true)
  
  MSE_Ridge_Beta_EM   <- .mseBeta(Mhat =  Beta_Ridge_EM,Btrue = Barray)
  MSE_Ridge_Pi_EM     <- .msepi(P_hat =  Pi_Ridge_EM[,1],P_true = Pi_true[1])
  MSE_Ridge_sigma_EM  <- .msesigma(sigma_hat =  Sigma_Ridge_EM,sigma_true = sigma_true)
  
  MSE_Liu_Beta_EM   <- .mseBeta(Mhat =  Beta_Liu_EM,Btrue = Barray)
  MSE_Liu_Pi_EM     <- .msepi(P_hat =  Pi_Liu_EM[,1],P_true = Pi_true[1])
  MSE_Liu_sigma_EM  <- .msesigma(sigma_hat =  Sigma_Liu_EM,sigma_true = sigma_true)
  
  ## MSE calculation CEM
  
  MSE_LS_Beta_CEM   <- .mseBeta(Mhat =  Beta_LS_CEM,Btrue = Barray)
  MSE_LS_Pi_CEM     <- .msepi(P_hat =  Pi_LS_CEM[,1],P_true = Pi_true[1])
  MSE_LS_sigma_CEM  <- .msesigma(sigma_hat =  Sigma_LS_CEM,sigma_true = sigma_true)
  
  MSE_Ridge_Beta_CEM   <- .mseBeta(Mhat =  Beta_Ridge_CEM,Btrue = Barray)
  MSE_Ridge_Pi_CEM     <- .msepi(P_hat =  Pi_Ridge_CEM[,1],P_true = Pi_true[1])
  MSE_Ridge_sigma_CEM  <- .msesigma(sigma_hat =  Sigma_Ridge_CEM,sigma_true = sigma_true)
  
  MSE_Liu_Beta_CEM   <- .mseBeta(Mhat =  Beta_Liu_CEM,Btrue = Barray)
  MSE_Liu_Pi_CEM     <- .msepi(P_hat =  Pi_Liu_CEM[,1],P_true = Pi_true[1])
  MSE_Liu_sigma_CEM  <- .msesigma(sigma_hat =  Sigma_Liu_CEM,sigma_true = sigma_true)
  
  ## MSE calculation SEM
  
  MSE_LS_Beta_SEM   <- .mseBeta(Mhat =  Beta_LS_SEM,Btrue = Barray)
  MSE_LS_Pi_SEM     <- .msepi(P_hat =  Pi_LS_SEM[,1],P_true = Pi_true[1])
  MSE_LS_sigma_SEM  <- .msesigma(sigma_hat =  Sigma_LS_SEM,sigma_true = sigma_true)
  
  MSE_Ridge_Beta_SEM   <- .mseBeta(Mhat =  Beta_Ridge_SEM,Btrue = Barray)
  MSE_Ridge_Pi_SEM     <- .msepi(P_hat =  Pi_Ridge_SEM[,1],P_true = Pi_true[1])
  MSE_Ridge_sigma_SEM  <- .msesigma(sigma_hat =  Sigma_Ridge_SEM,sigma_true = sigma_true)
  
  MSE_Liu_Beta_SEM   <- .mseBeta(Mhat =  Beta_Liu_SEM,Btrue = Barray)
  MSE_Liu_Pi_SEM     <- .msepi(P_hat =  Pi_Liu_SEM[,1],P_true = Pi_true[1])
  MSE_Liu_sigma_SEM  <- .msesigma(sigma_hat =  Sigma_Liu_SEM,sigma_true = sigma_true)
  
  ## bias calculation for EM
  Bias_LS_Beta_EM   <- .biasBeta(Mhat =  Beta_LS_EM,Btrue = Barray)
  Bias_LS_Pi_EM     <- .biaspi(P_hat =  Pi_LS_EM[,1],P_true = Pi_true[1])
  Bias_LS_sigma_EM  <- .biassigma(sigma_hat =  Sigma_LS_EM,sigma_true = sigma_true)
  
  Bias_Ridge_Beta_EM   <- .biasBeta(Mhat =  Beta_Ridge_EM,Btrue = Barray)
  Bias_Ridge_Pi_EM     <- .biaspi(P_hat =  Pi_Ridge_EM[,1],P_true = Pi_true[1])
  Bias_Ridge_sigma_EM  <- .biassigma(sigma_hat =  Sigma_Ridge_EM,sigma_true = sigma_true)
  
  Bias_Liu_Beta_EM   <- .biasBeta(Mhat =  Beta_Liu_EM,Btrue = Barray)
  Bias_Liu_Pi_EM     <- .biaspi(P_hat =  Pi_Liu_EM[,1],P_true = Pi_true[1])
  Bias_Liu_sigma_EM  <- .biassigma(sigma_hat =  Sigma_Liu_EM,sigma_true = sigma_true)
  
  ## bias calculation for CEM
  
  Bias_LS_Beta_CEM   <- .biasBeta(Mhat =  Beta_LS_CEM,Btrue = Barray)
  Bias_LS_Pi_CEM     <- .biaspi(P_hat =  Pi_LS_CEM[,1],P_true = Pi_true[1])
  Bias_LS_sigma_CEM  <- .biassigma(sigma_hat =  Sigma_LS_CEM,sigma_true = sigma_true)
  
  Bias_Ridge_Beta_CEM   <- .biasBeta(Mhat =  Beta_Ridge_CEM,Btrue = Barray)
  Bias_Ridge_Pi_CEM     <- .biaspi(P_hat =  Pi_Ridge_CEM[,1],P_true = Pi_true[1])
  Bias_Ridge_sigma_CEM  <- .biassigma(sigma_hat =  Sigma_Ridge_CEM,sigma_true = sigma_true)
  
  Bias_Liu_Beta_CEM   <- .biasBeta(Mhat =  Beta_Liu_CEM,Btrue = Barray)
  Bias_Liu_Pi_CEM     <- .biaspi(P_hat =  Pi_Liu_CEM[,1],P_true = Pi_true[1])
  Bias_Liu_sigma_CEM  <- .biassigma(sigma_hat =  Sigma_Liu_CEM,sigma_true = sigma_true)
  
  ## bias calculation for SEM
  
  Bias_LS_Beta_SEM   <- .biasBeta(Mhat =  Beta_LS_SEM,Btrue = Barray)
  Bias_LS_Pi_SEM     <- .biaspi(P_hat =  Pi_LS_SEM[,1],P_true = Pi_true[1])
  Bias_LS_sigma_SEM  <- .biassigma(sigma_hat =  Sigma_LS_SEM,sigma_true = sigma_true)
  
  Bias_Ridge_Beta_SEM   <- .biasBeta(Mhat =  Beta_Ridge_SEM,Btrue = Barray)
  Bias_Ridge_Pi_SEM     <- .biaspi(P_hat =  Pi_Ridge_SEM[,1],P_true = Pi_true[1])
  Bias_Ridge_sigma_SEM  <- .biassigma(sigma_hat =  Sigma_Ridge_SEM,sigma_true = sigma_true)
  
  Bias_Liu_Beta_SEM   <- .biasBeta(Mhat =  Beta_Liu_SEM,Btrue = Barray)
  Bias_Liu_Pi_SEM     <- .biaspi(P_hat =  Pi_Liu_SEM[,1],P_true = Pi_true[1])
  Bias_Liu_sigma_SEM  <- .biassigma(sigma_hat =  Sigma_Liu_SEM,sigma_true = sigma_true)
  
  tabMSE <- rbind(MSE_LS_Beta_EM,MSE_LS_Pi_EM,MSE_LS_sigma_EM,
                  MSE_Ridge_Beta_EM,MSE_Ridge_Pi_EM,MSE_Ridge_sigma_EM,
                  MSE_Liu_Beta_EM,MSE_Liu_Pi_EM,MSE_Liu_sigma_EM,
                  
                  MSE_LS_Beta_CEM,MSE_LS_Pi_CEM,MSE_LS_sigma_CEM,
                  MSE_Ridge_Beta_CEM,MSE_Ridge_Pi_CEM,MSE_Ridge_sigma_CEM,
                  MSE_Liu_Beta_CEM,MSE_Liu_Pi_CEM,MSE_Liu_sigma_CEM,
                  
                  MSE_LS_Beta_SEM,MSE_LS_Pi_SEM,MSE_LS_sigma_SEM,
                  MSE_Ridge_Beta_SEM,MSE_Ridge_Pi_SEM,MSE_Ridge_sigma_SEM,
                  MSE_Liu_Beta_SEM,MSE_Liu_Pi_SEM,MSE_Liu_sigma_SEM)
  
  colnames(tabMSE) <- c("Mean","Median","2.5%","97.5%","length")
  
  rownames(tabMSE) <- c("Beta_LS_EM","Pi_LS_EM","Sigma_LS_EM",
                        "Beta_Ridge_EM","Pi_Ridge_EM","Sigma_Ridge_EM",
                        "Beta_Liu_EM","Pi_Liu_EM","Sigma_Liu_EM",
                        
                        "Beta_LS_CEM","Pi_LS_CEM","Sigma_LS_CEM",
                        "Beta_Ridge_CEM","Pi_Ridge_CEM","Sigma_Ridge_CEM",
                        "Beta_Liu_CEM","Pi_Liu_CEM","Sigma_Liu_CEM",
                        
                        "Beta_LS_SEM","Pi_LS_SEM","Sigma_LS_SEM",
                        "Beta_Ridge_SEM","Pi_Ridge_SEM","Sigma_Ridge_SEM",
                        "Beta_Liu_SEM","Pi_Liu_SEM","Sigma_Liu_SEM")
  
  
  
  
  
  
  
   # MSE_LS <- c(MSE_LS_Beta_EM,MSE_LS_Pi_EM,MSE_LS_sigma_EM,
   #             MSE_LS_Beta_CEM,MSE_LS_Pi_CEM,MSE_LS_sigma_CEM,
   #             MSE_LS_Beta_SEM,MSE_LS_Pi_SEM,MSE_LS_sigma_SEM)
   # 
   # Bias_LS <- c(Bias_LS_Beta_EM,Bias_LS_Pi_EM,Bias_LS_sigma_EM,
   #           Bias_LS_Beta_CEM,Bias_LS_Pi_CEM,Bias_LS_sigma_CEM,
   #           Bias_LS_Beta_SEM,Bias_LS_Pi_SEM,Bias_LS_sigma_SEM)
   # 
   # MSE_Ridge <- c(MSE_Ridge_Beta_EM,MSE_Ridge_Pi_EM,MSE_Ridge_sigma_EM,
   #             MSE_Ridge_Beta_CEM,MSE_Ridge_Pi_CEM,MSE_Ridge_sigma_CEM,
   #             MSE_Ridge_Beta_SEM,MSE_Ridge_Pi_SEM,MSE_Ridge_sigma_SEM)
   # 
   # Bias_Ridge <- c(Bias_Ridge_Beta_EM,Bias_Ridge_Pi_EM,Bias_Ridge_sigma_EM,
   #           Bias_Ridge_Beta_CEM,Bias_Ridge_Pi_CEM,Bias_Ridge_sigma_CEM,
   #           Bias_Ridge_Beta_SEM,Bias_Ridge_Pi_SEM,Bias_Ridge_sigma_SEM)
   # 
   # MSE_Liu <- c(MSE_Liu_Beta_EM,MSE_Liu_Pi_EM,MSE_Liu_sigma_EM,
   #              MSE_Liu_Beta_CEM,MSE_Liu_Pi_CEM,MSE_Liu_sigma_CEM,
   #              MSE_Liu_Beta_SEM,MSE_Liu_Pi_SEM,MSE_Liu_sigma_SEM)
   # 
   # Bias_Liu <- c(Bias_Liu_Beta_EM,Bias_Liu_Pi_EM,Bias_Liu_sigma_EM,
   #               Bias_Liu_Beta_CEM,Bias_Liu_Pi_CEM,Bias_Liu_sigma_CEM,
   #               Bias_Liu_Beta_SEM,Bias_Liu_Pi_SEM,Bias_Liu_sigma_SEM)
   # 
   # tabMSE <- cbind(MSE_LS,Bias_LS,MSE_Ridge,Bias_Ridge,MSE_Liu,Bias_Liu )
   # 
   # colnames(tabMSE) <- c("MSE_LS","Bias_LS","MSE_Ridge","Bias_Ridge","MSE_Liu","Bias_Liu")
   # rownames(tabMSE) <- c("Beta_EM","Pi_EM","Sigma_EM","Beta_CEM","Pi_CEM","Sigma_CEM","Beta_SEM","Pi_SEM","Sigma_SEM")
   # 

  #****************************************************************** 
  
  #  Beta_vec  <- c(MSE_LS_Beta_EM,MSE_Ridge_Beta_EM,MSE_LS_Beta_CEM,MSE_Ridge_Beta_CEM,MSE_LS_Beta_SEM,MSE_Ridge_Beta_SEM,
  #                 Bias_LS_Beta_EM,Bias_Ridge_Beta_EM,Bias_LS_Beta_CEM,Bias_Ridge_Beta_CEM,Bias_LS_Beta_SEM,Bias_Ridge_Beta_SEM)
  # 
  #  Pi_vec    <- c(MSE_LS_Pi_EM,MSE_Ridge_Pi_EM,MSE_LS_Pi_CEM,MSE_Ridge_Pi_CEM,MSE_LS_Pi_SEM,MSE_Ridge_Pi_SEM,
  #                 Bias_LS_Pi_EM,Bias_Ridge_Pi_EM,Bias_LS_Pi_CEM,Bias_Ridge_Pi_CEM,Bias_LS_Pi_SEM,Bias_Ridge_Pi_SEM)
  # 
  #  Sigma_vec <- c(MSE_LS_sigma_EM,MSE_Ridge_sigma_EM,MSE_LS_sigma_CEM,MSE_Ridge_sigma_CEM,MSE_LS_sigma_SEM,MSE_Ridge_sigma_SEM,
  #                 Bias_LS_sigma_EM,Bias_Ridge_sigma_EM,Bias_LS_sigma_CEM,Bias_Ridge_sigma_CEM,Bias_LS_sigma_SEM,Bias_Ridge_sigma_SEM)
  # 
  # tabMSE <- rbind(Beta_vec,Pi_vec,Sigma_vec)
  # colnames(tabMSE) <- c("MSE_EM","MSE_Ridge_EM","MSE_CEM","MSE_Ridge_CEM","MSE_SEM","MSE_Ridge_SEM",
  #                       "Bias_EM","Bias_Ridge_EM","Bias_CEM","Bias_Ridge_CEM","Bias_SEM","Bias_Ridge_SEM")
  # rownames(tabMSE) <- c("Beta","Pi","Sigma")
  
  return(tabMSE)
}


MSE_Result <-Wrapping(Beta_LS_EM,Pi_LS_EM,Sigma_LS_EM,
                      Beta_LS_CEM,Pi_LS_CEM,Sigma_LS_CEM,
                      Beta_LS_SEM,Pi_LS_SEM,Sigma_LS_SEM,
                      Beta_Ridge_EM,Pi_Ridge_EM,Sigma_Ridge_EM,
                      Beta_Ridge_CEM,Pi_Ridge_CEM,Sigma_Ridge_CEM,
                      Beta_Ridge_SEM,Pi_Ridge_SEM,Sigma_Ridge_SEM,
                      Beta_Liu_EM,Pi_Liu_EM,Sigma_Liu_EM,
                      Beta_Liu_CEM,Pi_Liu_CEM,Sigma_Liu_CEM,
                      Beta_Liu_SEM,Pi_Liu_SEM,Sigma_Liu_SEM,
                      Beta0,pi0,sigma0)
#MSE_Result

# res <- Estimate(Beta_LS_EM,Beta_LS_CEM,Beta_LS_SEM,Pi_LS_EM,Pi_LS_CEM,Pi_LS_SEM,
#                 Sigma_LS_EM,Sigma_LS_CEM,Sigma_LS_SEM,Beta0,pi0,sigma0)
# res

# cat("The mean number of iteration required for EM for ML is:",mean(K_LS_EM),'\n')
# cat("The mean number of iteration required for CEM for ML is:",mean(K_LS_CEM),'\n')
# cat("The mean number of iteration required for SEM for ML is:",mean(K_LS_SEM),'\n')
# cat("The mean number of iteration required for EM for Ridge is:",mean(K_Ridge_EM),'\n')
# cat("The mean number of iteration required for CEM for Ridge is:",mean(K_Ridge_CEM),'\n')
# cat("The mean number of iteration required for SEM for Ridge is:",mean(K_Ridge_SEM),'\n')

df_Iterantion<-data.frame(EM=c(mean(K_LS_EM),mean(K_Ridge_EM),mean(K_Liu_EM)),CEM=c(mean(K_LS_CEM),mean(K_Ridge_CEM),mean(K_Liu_CEM)),
                          SEM=c(mean(K_LS_SEM),mean(K_Ridge_SEM),mean(K_Liu_SEM)),row.names = c("ML","Ridge","Liu"))


##**************************************************************************

Calc_RMSEP <- function(Y_test,Y_Hat)
{
  dif1 <- ((Y_test-Y_Hat)^2)
  diff2 <- mean(dif1)
  diff3 <- sqrt(diff2)
  return(diff3)
}

Calc_MRSEP <- function(RMSEP)
{
  Mean_MRSEP <- round(mean(RMSEP),3)
  quanMRSEP <- round(quantile(RMSEP, probs=c(0.025,0.5,0.975)),3)
  quanlen <- quanMRSEP[3]-quanMRSEP[1]
  tab <- c(Mean_MRSEP,quanMRSEP[2],quanMRSEP[1],quanMRSEP[3],quanlen)
  return(tab)
}

##**************************************************************************


## Predictions

## Args:
## Data: Matrix witch contains the independent variables and response variables for NSim
## n: sample size
## K: Number of cross validation
## Nsim: the number of simulation
## epsilon <- the error for stopping rule
## beta.0,beta.old,pi.0,pi.old,sigma.0,sigma.old: values of starting point


Prediction <- function(Data,k,n,Nsim,p,epsilon,M,beta.0,beta.old,pi.0,pi.old,sigma.0,sigma.old,maxit_EM)
{
  RMSEP_LS_EM <- RMSEP_Ridge_EM <- RMSEP_Liu_EM <- c()
  RMSEP_LS_CEM <- RMSEP_Ridge_CEM <- RMSEP_Liu_CEM <- c()
  RMSEP_LS_SEM <- RMSEP_Ridge_SEM <- RMSEP_Liu_SEM <- c()
  
  for (i in 1:Nsim)
  {
    Mydata=Data[,,i]
    folds <- cut(seq(1,nrow(Mydata)),breaks = k,labels=FALSE)
    Y_test_1 <- array(NA,dim=c(n/k,1,k))
    Y_Hat_LS_EM <- Y_Hat_Ridge_EM <- Y_Hat_Liu_EM <- array(NA,dim=c(n/k,1,k))
    Y_Hat_LS_CEM <- Y_Hat_Ridge_CEM <- Y_Hat_Liu_CEM <- array(NA,dim=c(n/k,1,k)) 
    Y_Hat_LS_SEM <- Y_Hat_Ridge_SEM <- Y_Hat_Liu_SEM <- array(NA,dim=c(n/k,1,k)) 
    
    for (L in 1:k)
    {
      testIndex <- which(folds==L,arr.ind = T)
      testData <- Mydata[testIndex,]
      traindata <- Mydata[-testIndex,]
      X_train <- traindata[,-1] 
      Y_train <- traindata[,1]
      X_test <- testData[,-1]
      Y_test <- testData[,1]
      
      Y_test_1[,,L] <- Y_test
      
      # EM Alogorithm
      
      EM_Result_LS <- em_ml(epsilon=epsilon,x=X_train,y=Y_train,beta.0=beta.0,
                            beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                            sigma.old=sigma.old,maxit_EM=maxit_EM,printout=F)
      
      #K_LS_EM <- EM_Result_LS$k
      P_ml_EM <- EM_Result_LS$pi
      Sigma_LS_EM <- EM_Result_LS$sigma
      B_ml_EM <- EM_Result_LS$beta
      
      Yhat_LS_EM <- Design_Response_prediction(x=X_test,beta=B_ml_EM,pi=P_ml_EM,M=M)
      Y_Hat_LS_EM[,,L] <- Yhat_LS_EM[,1]
      
      ## Ridge Method EM
      
      k_ridge <- K.hkb(p,Sigma_LS_EM,B_ml_EM)
      
      EM_Result_ridge <- em_ridge(epsilon=epsilon,x=X_train,y=Y_train,beta.0=beta.0,
                                  beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                                  sigma.old=sigma.old,maxit_EM=maxit_EM,k_ridge,printout=F)
      
      #K_Ridge_EM <- EM_Result_ridge$k
      P_Ridge_EM <- EM_Result_ridge$pi
      Sigma_Ridge_EM <- EM_Result_ridge$sigma
      B_Ridge_EM <- EM_Result_ridge$beta
      
      Yhat_Ridge_EM <- Design_Response_prediction(x=X_test,beta=B_Ridge_EM,pi=P_Ridge_EM,M=M)
      Y_Hat_Ridge_EM[,,L] <- Yhat_Ridge_EM[,1]
      
      
      ## Liu EM
      
      k_Liu <- K.hkb(p,Sigma_Ridge_EM,B_Ridge_EM )
      
      EM_Result_Liu <- em_Liu(epsilon=epsilon,x=X_train,y=Y_train,beta.0=beta.0,
                              beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                              sigma.old=sigma.old,maxit_EM=maxit_EM,k_Liu,Sigma_Ridge_EM,printout=F)
      
      #K_Liu_EM <- EM_Result_Liu$k
      P_Liu_EM <- EM_Result_Liu$pi
      Sigma_Liu_EM <- EM_Result_Liu$sigma
      B_Liu_EM <- EM_Result_Liu$beta
      
      Yhat_Liu_EM <- Design_Response_prediction(x=X_test,beta=B_Liu_EM,pi=P_Liu_EM,M=M)
      Y_Hat_Liu_EM[,,L] <- Yhat_Liu_EM[,1]
      
      
      ## Methods of CEM
      
      ## ML Method
      
      CEM_Result_LS <- CEM_ml(epsilon=epsilon,x=X_train,y=Y_train,beta.0=beta.0,
                              beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                              sigma.old=sigma.old,maxit_EM=maxit_EM,printout=F)
      
      #K_ml_CEM <- CEM_Result_LS$k
      P_ml_CEM <- CEM_Result_LS$pi
      Sigma_ml_CEM <- CEM_Result_LS$sigma
      B_ml_CEM <- CEM_Result_LS$beta
      
      Yhat_LS_CEM <- Design_Response_prediction(x=X_test,beta=B_ml_CEM,pi=P_ml_CEM,M=M)
      Y_Hat_LS_CEM[,,L] <- Yhat_LS_CEM[,1]
      
      
      ## Ridge Method
      
      k_ridge <- K.hkb(p,Sigma_ml_CEM,B_ml_CEM)
      
      CEM_Result_ridge <- CEM_ridge(epsilon=epsilon,x=X_train,y=Y_train,beta.0=beta.0,
                                    beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                                    sigma.old=sigma.old,maxit_EM=maxit_EM,k_ridge,printout=F)
      
      #K_Ridge_CEM <- CEM_Result_ridge$k
      P_Ridge_CEM <- CEM_Result_ridge$pi
      Sigma_Ridge_CEM <- CEM_Result_ridge$sigma
      B_Ridge_CEM <- CEM_Result_ridge$beta
      
      Yhat_Ridge_CEM <- Design_Response_prediction(x=X_test,beta=B_Ridge_CEM,pi=P_Ridge_CEM,M=M)
      Y_Hat_Ridge_CEM[,,L] <- Yhat_Ridge_CEM[,1]
      
      
      ## Liu Method
      
      k_Liu <- K.hkb(p,Sigma_Ridge_CEM,B_Ridge_CEM )
      
      CEM_Result_Liu <- CEM_Liu(epsilon=epsilon,x=X_train,y=Y_train,beta.0=beta.0,
                                beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                                sigma.old=sigma.old,maxit_EM=maxit_EM,k_Liu,Sigma_Ridge_CEM,printout=F)
      
      #K_Liu_CEM <- CEM_Result_Liu$k
      P_Liu_CEM <- CEM_Result_Liu$pi
      Sigma_Liu_CEM <- CEM_Result_Liu$sigma
      B_Liu_CEM <- CEM_Result_Liu$beta
      
      Yhat_Liu_CEM <- Design_Response_prediction(x=X_test,beta=B_Liu_CEM,pi=P_Liu_CEM,M=M)
      Y_Hat_Liu_CEM[,,L] <- Yhat_Liu_CEM[,1]
      
      
      ## Methods of SEM
      
      ## ML Method
      
      SEM_Result_LS <- SEM_ml(epsilon=epsilon,x=X_train,y=Y_train,beta.0=beta.0,
                              beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                              sigma.old=sigma.old,maxit_EM=maxit_EM,printout=F)
      
      #K_ml_SEM <- SEM_Result_LS$k
      P_ml_SEM <- SEM_Result_LS$pi
      Sigma_ml_SEM <- SEM_Result_LS$sigma
      B_ml_SEM <- SEM_Result_LS$beta
      
      Yhat_LS_SEM <- Design_Response_prediction(x=X_test,beta=B_ml_SEM,pi=P_ml_SEM,M=M)
      Y_Hat_LS_SEM[,,L] <- Yhat_LS_SEM[,1]
      
      
      ## Ridge Method
      
      k_ridge <- K.hkb(p,Sigma_ml_SEM,B_ml_SEM)
      
      SEM_Result_ridge <- SEM_ridge(epsilon=epsilon,x=X_train,y=Y_train,beta.0=beta.0,
                                    beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                                    sigma.old=sigma.old,maxit_EM=maxit_EM,k_ridge,printout=F)
      
      #K_Ridge_SEM <- SEM_Result_ridge$k
      P_Ridge_SEM <- SEM_Result_ridge$pi
      Sigma_Ridge_SEM <- SEM_Result_ridge$sigma
      B_Ridge_SEM <- SEM_Result_ridge$beta
      
      Yhat_Ridge_SEM <- Design_Response_prediction(x=X_test,beta=B_Ridge_SEM,pi=P_Ridge_SEM,M=M)
      Y_Hat_Ridge_SEM[,,L] <- Yhat_Ridge_SEM[,1]
      
      
      ## Liu Method
      
      k_Liu <- K.hkb(p,Sigma_Ridge_SEM,B_Ridge_SEM )
      
      SEM_Result_Liu <- SEM_Liu(epsilon=epsilon,x=X_train,y=Y_train,beta.0=beta.0,
                                beta.old=beta.old,pi.0=pi.0,pi.old=pi.old,sigma.0=sigma.0,
                                sigma.old=sigma.old,maxit_EM=maxit_EM,k_Liu,Sigma_Ridge_SEM,printout=F)
      
      #K_Liu_SEM <- SEM_Result_Liu$k
      P_Liu_SEM <- SEM_Result_Liu$pi
      Sigma_Liu_SEM <- SEM_Result_Liu$sigma
      B_Liu_SEM <- SEM_Result_Liu$beta
      
      Yhat_Liu_SEM <- Design_Response_prediction(x=X_test,beta=B_Liu_SEM,pi=P_Liu_SEM,M=M)
      Y_Hat_Liu_SEM[,,L] <- Yhat_Liu_SEM[,1]
      
      
    }
    
    RMSEP_LS_EM[i] <- Calc_RMSEP(Y_test_1,Y_Hat_LS_EM)
    RMSEP_Ridge_EM[i] <- Calc_RMSEP(Y_test_1,Y_Hat_Ridge_EM)
    RMSEP_Liu_EM[i] <- Calc_RMSEP(Y_test_1,Y_Hat_Liu_EM)
    
    RMSEP_LS_CEM[i] <- Calc_RMSEP(Y_test_1,Y_Hat_LS_CEM)
    RMSEP_Ridge_CEM[i] <- Calc_RMSEP(Y_test_1,Y_Hat_Ridge_CEM)
    RMSEP_Liu_CEM[i] <- Calc_RMSEP(Y_test_1,Y_Hat_Liu_CEM)
    
    RMSEP_LS_SEM[i] <- Calc_RMSEP(Y_test_1,Y_Hat_LS_SEM)
    RMSEP_Ridge_SEM[i] <- Calc_RMSEP(Y_test_1,Y_Hat_Ridge_SEM)
    RMSEP_Liu_SEM[i] <- Calc_RMSEP(Y_test_1,Y_Hat_Liu_SEM)
    
  }
  
  RMSEP_LS_EM <- Calc_MRSEP(RMSEP_LS_EM)
  RMSEP_Ridge_EM <- Calc_MRSEP(RMSEP_Ridge_EM)
  RMSEP_Liu_EM <- Calc_MRSEP(RMSEP_Liu_EM)
  
  RMSEP_LS_CEM <- Calc_MRSEP(RMSEP_LS_CEM)
  RMSEP_Ridge_CEM <- Calc_MRSEP(RMSEP_Ridge_CEM)
  RMSEP_Liu_CEM <- Calc_MRSEP(RMSEP_Liu_CEM)
  
  RMSEP_LS_SEM <- Calc_MRSEP(RMSEP_LS_SEM)
  RMSEP_Ridge_SEM <- Calc_MRSEP(RMSEP_Ridge_SEM)
  RMSEP_Liu_SEM <- Calc_MRSEP(RMSEP_Liu_SEM)
  
  tabMSE <- rbind(RMSEP_LS_EM,RMSEP_Ridge_EM,RMSEP_Liu_EM,
                  RMSEP_LS_CEM,RMSEP_Ridge_CEM,RMSEP_Liu_CEM,
                  RMSEP_LS_SEM,RMSEP_Ridge_SEM,RMSEP_Liu_SEM)
  
  colnames(tabMSE) <- c("Mean","Median","2.5%","97.5%","length")
  
  rownames(tabMSE) <- c("RMSEP_LS_EM","RMSEP_Ridge_EM","RMSEP_Liu_EM",
                        "RMSEP_LS_CEM","RMSEP_Ridge_CEM","RMSEP_Liu_CEM",
                        "RMSEP_LS_SEM","RMSEP_Ridge_SEM","RMSEP_Liu_SEM")
  # df <- data.frame(MRSEP_LS=c(mean(RMSEP_LS_EM),mean(RMSEP_LS_CEM),mean(RMSEP_LS_SEM)),
  #                  MRSEP_Ridge=c(mean(RMSEP_Ridge_EM),mean(RMSEP_Ridge_CEM),mean(RMSEP_Ridge_SEM)),
  #                  MRSEP_Liu=c(mean(RMSEP_Liu_EM),mean(RMSEP_Liu_CEM),mean(RMSEP_Liu_SEM)),row.names = c("EM","CEM","SEM"))
  # 
  # 
  return(tabMSE)
}

MSE_Result

pre <- Prediction(Data,k_fold,n0t,Nsim,p=p0,epsilon=epsilon_em0,M=M0,
                  beta.0=Beta.01,beta.old=Beta.02,pi.0=pi.01,pi.old=pi.02,sigma.0=sigma.01
                  ,sigma.old=sigma.02,maxit_EM=maxit_em0)
MSE_Result
pre 

print("The mean number of iteration required for each method is:")

df_Iterantion

save(MSE_Result,file="mse_size60_HKP.RData")
save(pre,file="test_size60_HKP.RData")




