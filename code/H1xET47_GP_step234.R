# '@ Alain Mbebi
set.seed(123)
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# '@ libraries
library(data.table)
library(MASS)
library(corpcor)
library(Matrix)          
library(lattice)
library(mvtnorm)
library(glasso)          
library(glassoFast)
library(matrixcalc)
library(glmnet)
library(BGLR)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# '@functions
#function to center data matrices using 'colMeans()' note that we center the col because they represent variables
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function to compute the L2-norm of a given vector x
L2norm <- function(x) {
  sqrt(crossprod(x))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for the estimation in mutiple ridge for B
multi.ridge <- function(X,Y, lam){
	q=dim(Y)[2]
  p=dim(X)[2] 
  Bhat = matrix(0, nrow=p, ncol=q)
	for(kk in 1:q)
	{
	  Bhat[,kk]=as.numeric(glmnet(x=X, y=Y[,kk], family="gaussian", alpha=0,lambda=lam, standardize=FALSE)$beta)
	}
  return(Bhat)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for cross validation in multiple ridge
multi.ridge.cv <-function(X,Y, lambdaB, kfold=3, silent=TRUE){
  set.seed(123)
  n.out = nrow(Y)
  K.out = kfold
  d.out = ceiling(n.out/K.out)
  set.seed(123)
  i.mix.out = sample(1:n.out)
  folds.out = vector(mode="list", length=K.out) 
  for (k in 1:(K.out-1)) {
    folds.out[[k]] = i.mix.out[((k-1)*d.out+1):(k*d.out)]
  }
  folds.out[[K.out]] = i.mix.out[((K.out-1)*d.out+1):n.out]
  lambdaB = sort(lambdaB)
  CV_errors.out = rep(NA, K.out)
  CV_errors.val.out = rep(list(CV_errors.out), K.out)
  
  CV_corr.out = matrix(NA, nrow=ncol(Y), ncol=ncol(Y))
  CV_corr.val.out = rep(list(CV_corr.out), K.out)
  
  opt.i=c()
  opt.lamB.out = c()
  
  for (l in 1:K.out) {
    cat("outer Fold",l,"\n")
    i.tr.out = unlist(folds.out[-k])
    i.val.out = folds.out[[k]]
    
    X.tr.out = X[i.tr.out,]    # training predictors
    Y.tr.out = Y[i.tr.out,]    # training responses
    
    X.val.out = X[i.val.out,]  # validation predictors
    Y.val.out = Y[i.val.out,]  # validation responses
    
    #Inner loop for hyper params  
    n.iner = nrow(Y.tr.out)
    K.iner = 3
    d.iner = ceiling(n.iner/K.iner)
    set.seed(123)
    i.mix.iner = sample(1:n.iner)
    folds.iner = vector(mode="list", length=K.iner) 
    for (k in 1:(K.iner-1)) {
      folds.iner[[k]] = i.mix.iner[((k-1)*d.iner+1):(k*d.iner)]
    }
    folds.iner[[K.iner]] = i.mix.iner[((K.iner-1)*d.iner+1):n.iner]  

    CV_errors.iner = rep(NA, length(lambdaB))
    CV_errors.val.iner = rep(list(CV_errors.iner), K.iner)

    for (k in 1:K.iner) {
      cat("inner Fold",k,"\n")
      
      i.tr.iner = unlist(folds.iner[-k])
      i.val.iner = folds.iner[[k]]
      
      X.tr.iner  = X.tr.out[i.tr.iner,]    # training predictors in the inner loop 
      Y.tr.iner  = Y.tr.out[i.tr.iner,]    # training responses in the inner loop 
      X.val.iner = X.tr.out[i.val.iner,]   # validation predictors in the inner loop 
      Y.val.iner = Y.tr.out[i.val.iner,]   # validation responses in the inner loop 
      
        for(i in 1:length(lambdaB)){
          B_hat.tr.iner = multi.ridge(X.tr.iner, Y.tr.iner, lambdaB[i]) 
          CV_errors.val.iner[[k]][i] = mean(mse.matrix(X.val.iner%*%B_hat.tr.iner,Y.val.iner))
        }
      
      idx.min=which.min((sapply(CV_errors.val.iner, min)))
      opt = which.min(CV_errors.val.iner[[idx.min]])      
      opt.i[k] = opt
      opt.lamB.out[k] = lambdaB[opt.i[k]]
      
    }
    
        B_hat.tr.out = multi.ridge(X.tr.out, Y.tr.out, opt.lamB.out[l]) 
        CV_errors.val.out[[l]] = mse.matrix(X.val.out%*%B_hat.tr.out, Y.val.out, na.rm = TRUE)
        CV_corr.val.out[[l]] = cor(X.val.out%*%B_hat.tr.out, Y.val.out)
  }
  
  best.idx.out=which.min((sapply(CV_errors.val.out, min)))
  best.lamB.final = opt.lamB.out[best.idx.out]
  B_hat_final = multi.ridge(X, Y, best.lamB.final) 
  return(list(best.lamB.final=best.lamB.final, B_hat_final=B_hat_final, cv_corr=CV_corr.val.out, cv.err=CV_errors.val.out)) 
   
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for the estimation in mutiple lasso for B
multi.lasso <- function(X,Y, lam){
	q=dim(Y)[2]
  p=dim(X)[2] 
  Bhat = matrix(0, nrow=p, ncol=q)
	for(kk in 1:q)
	{
	  Bhat[,kk]=as.numeric(glmnet(x=X, y=Y[,kk], family="gaussian", alpha=1,lambda=lam, standardize=FALSE)$beta)
	}
  return(Bhat)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for cross validation in multiple lasso
multi.lasso.cv <-function(X,Y, lambdaB, kfold=3, silent=TRUE){
  set.seed(123)
  n.out = nrow(Y)
  K.out = kfold
  d.out = ceiling(n.out/K.out)
  set.seed(123)
  i.mix.out = sample(1:n.out)
  folds.out = vector(mode="list", length=K.out) 
  for (k in 1:(K.out-1)) {
    folds.out[[k]] = i.mix.out[((k-1)*d.out+1):(k*d.out)]
  }
  folds.out[[K.out]] = i.mix.out[((K.out-1)*d.out+1):n.out]
  lambdaB = sort(lambdaB)
  CV_errors.out = rep(NA, K.out)
  CV_errors.val.out = rep(list(CV_errors.out), K.out)
  
  CV_corr.out = matrix(NA, nrow=ncol(Y), ncol=ncol(Y))
  CV_corr.val.out = rep(list(CV_corr.out), K.out)
  
  opt.i=c()
  opt.lamB.out = c()
  
  for (l in 1:K.out) {
    cat("outer Fold",l,"\n")
    i.tr.out = unlist(folds.out[-k])
    i.val.out = folds.out[[k]]
    
    X.tr.out = X[i.tr.out,]   
    Y.tr.out = Y[i.tr.out,]    
    
    X.val.out = X[i.val.out,] 
    Y.val.out = Y[i.val.out,]  
    n.iner = nrow(Y.tr.out)
    K.iner = 3
    d.iner = ceiling(n.iner/K.iner)
    set.seed(123)
    i.mix.iner = sample(1:n.iner)
    folds.iner = vector(mode="list", length=K.iner) 
    for (k in 1:(K.iner-1)) {
      folds.iner[[k]] = i.mix.iner[((k-1)*d.iner+1):(k*d.iner)]
    }
    folds.iner[[K.iner]] = i.mix.iner[((K.iner-1)*d.iner+1):n.iner]  
    CV_errors.iner = rep(NA, length(lambdaB))
    CV_errors.val.iner = rep(list(CV_errors.iner), K.iner)

    for (k in 1:K.iner) {
      cat("inner Fold",k,"\n")
      
      i.tr.iner = unlist(folds.iner[-k])
      i.val.iner = folds.iner[[k]]
      
      X.tr.iner  = X.tr.out[i.tr.iner,]   
      Y.tr.iner  = Y.tr.out[i.tr.iner,]    
      X.val.iner = X.tr.out[i.val.iner,] 
      Y.val.iner = Y.tr.out[i.val.iner,] 

        for(i in 1:length(lambdaB)){
          B_hat.tr.iner = multi.lasso(X.tr.iner, Y.tr.iner, lambdaB[i]) 
          CV_errors.val.iner[[k]][i] = mean(mse.matrix(X.val.iner%*%B_hat.tr.iner,Y.val.iner))
        }
      
      idx.min=which.min((sapply(CV_errors.val.iner, min)))
      opt = which.min(CV_errors.val.iner[[idx.min]])      
      opt.i[k] = opt
      opt.lamB.out[k] = lambdaB[opt.i[k]]
      
    }

        B_hat.tr.out = multi.lasso(X.tr.out, Y.tr.out, opt.lamB.out[l]) 
        CV_errors.val.out[[l]] = mse.matrix(X.val.out%*%B_hat.tr.out, Y.val.out, na.rm = TRUE)
        CV_corr.val.out[[l]] = cor(X.val.out%*%B_hat.tr.out, Y.val.out)
  }
  
  best.idx.out=which.min((sapply(CV_errors.val.out, min)))
  best.lamB.final = opt.lamB.out[best.idx.out]
  B_hat_final = multi.lasso(X, Y, best.lamB.final) 
  return(list(best.lamB.final=best.lamB.final, B_hat_final=B_hat_final, cv_corr=CV_corr.val.out, cv.err=CV_errors.val.out)) 
   
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for the estimation in mutiple elastic net for B
multi.elasnet <- function(X,Y, lam){
	q=dim(Y)[2]
  p=dim(X)[2] 
  Bhat = matrix(0, nrow=p, ncol=q)
	for(kk in 1:q)
	{
	  Bhat[,kk]=as.numeric(glmnet(x=X, y=Y[,kk], family="gaussian", alpha=.5,lambda=lam, standardize=FALSE)$beta)
	}
  return(Bhat)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for cross validation EN
multi.elasnet.cv <-function(X,Y, lambdaB, kfold=3, silent=TRUE){
  set.seed(123)
  n.out = nrow(Y)
  K.out = kfold
  d.out = ceiling(n.out/K.out)
  set.seed(123)
  i.mix.out = sample(1:n.out)
  folds.out = vector(mode="list", length=K.out) 
  for (k in 1:(K.out-1)) {
    folds.out[[k]] = i.mix.out[((k-1)*d.out+1):(k*d.out)]
  }
  folds.out[[K.out]] = i.mix.out[((K.out-1)*d.out+1):n.out]
  lambdaB = sort(lambdaB)
  CV_errors.out = rep(NA, K.out)
  CV_errors.val.out = rep(list(CV_errors.out), K.out)
  
  CV_corr.out = matrix(NA, nrow=ncol(Y), ncol=ncol(Y))
  CV_corr.val.out = rep(list(CV_corr.out), K.out)
  
  opt.i=c()
  opt.lamB.out = c()
  
  for (l in 1:K.out) {
    cat("outer Fold",l,"\n")
    i.tr.out = unlist(folds.out[-k])
    i.val.out = folds.out[[k]]
    
    X.tr.out = X[i.tr.out,]   
    Y.tr.out = Y[i.tr.out,]    
    
    X.val.out = X[i.val.out,] 
    Y.val.out = Y[i.val.out,]  
    
    n.iner = nrow(Y.tr.out)
    K.iner = 3
    d.iner = ceiling(n.iner/K.iner)
    set.seed(123)
    i.mix.iner = sample(1:n.iner)
    folds.iner = vector(mode="list", length=K.iner) 
    for (k in 1:(K.iner-1)) {
      folds.iner[[k]] = i.mix.iner[((k-1)*d.iner+1):(k*d.iner)]
    }
    folds.iner[[K.iner]] = i.mix.iner[((K.iner-1)*d.iner+1):n.iner]  
    CV_errors.iner = rep(NA, length(lambdaB))
    CV_errors.val.iner = rep(list(CV_errors.iner), K.iner)

    for (k in 1:K.iner) {
      cat("inner Fold",k,"\n")
      
      i.tr.iner = unlist(folds.iner[-k])
      i.val.iner = folds.iner[[k]]
      
      X.tr.iner  = X.tr.out[i.tr.iner,]   
      Y.tr.iner  = Y.tr.out[i.tr.iner,]    
      X.val.iner = X.tr.out[i.val.iner,] 
      Y.val.iner = Y.tr.out[i.val.iner,]  

        for(i in 1:length(lambdaB)){
          B_hat.tr.iner = multi.elasnet(X.tr.iner, Y.tr.iner, lambdaB[i]) 
          CV_errors.val.iner[[k]][i] = mean(mse.matrix(X.val.iner%*%B_hat.tr.iner,Y.val.iner))
        }
      
      idx.min=which.min((sapply(CV_errors.val.iner, min)))
      opt = which.min(CV_errors.val.iner[[idx.min]])      
      opt.i[k] = opt
      opt.lamB.out[k] = lambdaB[opt.i[k]]
      
    }
        B_hat.tr.out = multi.elasnet(X.tr.out, Y.tr.out, opt.lamB.out[l]) 
        CV_errors.val.out[[l]] = mse.matrix(X.val.out%*%B_hat.tr.out, Y.val.out, na.rm = TRUE)
        CV_corr.val.out[[l]] = cor(X.val.out%*%B_hat.tr.out, Y.val.out)
  }
  
  best.idx.out=which.min((sapply(CV_errors.val.out, min)))
  best.lamB.final = opt.lamB.out[best.idx.out]
  B_hat_final = multi.elasnet(X, Y, best.lamB.final) 
  return(list(best.lamB.final=best.lamB.final, B_hat_final=B_hat_final, cv_corr=CV_corr.val.out, cv.err=CV_errors.val.out)) 
   
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for cross validation in Multitrait Bayes B using the R package BGLR
multi.BayesB.cv <-function(X,Y, kfold=3){
  set.seed(123)
  n.out = nrow(Y)
  K.out = kfold
  d.out = ceiling(n.out/K.out)
  set.seed(123)
  i.mix.out = sample(1:n.out)
  folds.out = vector(mode="list", length=K.out) 
  for (k in 1:(K.out-1)) {
    folds.out[[k]] = i.mix.out[((k-1)*d.out+1):(k*d.out)]
  }
  folds.out[[K.out]] = i.mix.out[((K.out-1)*d.out+1):n.out]

  SVD=list()
  U=list()
  D=list()
  V=list()
  BETA_MBayesB=list()
  Bhat_MBayesB=list()
  
  CV_errors.out = rep(NA, K.out)
  CV_errors.val.out = rep(list(CV_errors.out), K.out)
  
  CV_corr.out = matrix(NA, nrow=ncol(Y), ncol=ncol(Y))
  CV_corr.val.out = rep(list(CV_corr.out), K.out)
  
  for (l in 1:K.out) {
    cat("outer Fold",l,"\n")
    i.tr.out = unlist(folds.out[-k])
    i.val.out = folds.out[[k]]
    
    X.tr.out = X[i.tr.out,]    
    Y.tr.out = Y[i.tr.out,]    
    
    X.val.out = X[i.val.out,]  
    Y.val.out = Y[i.val.out,]  
    
    # Singular value decompositions Y=UDV' Bhat_MBayesB
    SVD[[l]]=svd(Y.tr.out)
    U[[l]]=SVD[[l]]$u
    D[[l]]=diag(SVD[[l]]$d)
    V[[l]]=SVD[[l]]$v
  
  
    BETA_MBayesB[[l]]=matrix(nrow=ncol( X.tr.out),ncol=ncol(Y.tr.out))
  
   ETA=list(list(X=X.tr.out, model='BayesB'))
  for(i in 1:ncol(Y.tr.out)){
    fm=BGLR(y=U[[l]][,i],ETA=ETA,verbose=F) 
    BETA_MBayesB[[l]][,i]=fm$ETA[[1]]$b
  }
  
  # Rotating coefficients to put them in marker space
  Bhat_MBayesB[[l]]=BETA_MBayesB[[l]]%*%D[[l]]%*%t(SVD[[l]]$v)
  
 
        CV_errors.val.out[[l]] = mse.matrix(X.val.out%*%Bhat_MBayesB[[l]], Y.val.out, na.rm = TRUE)
        CV_corr.val.out[[l]] = cor(X.val.out%*%Bhat_MBayesB[[l]], Y.val.out)
  }

  # return best lambdaB and other meaningful values
  return(list(B_hat_final=Bhat_MBayesB, cv_corr=CV_corr.val.out, cv.err=CV_errors.val.out)) 
   
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for cross validation in GBLUP, (marker based) using the R package rrBLUP
gblup.cv <-function(X,Y, kfold=3){
  set.seed(123)
  n.out = nrow(Y)
  K.out = kfold
  d.out = ceiling(n.out/K.out)
  set.seed(123)
  i.mix.out = sample(1:n.out)
  folds.out = vector(mode="list", length=K.out) 
  for (k in 1:(K.out-1)) {
    folds.out[[k]] = i.mix.out[((k-1)*d.out+1):(k*d.out)]
  }
  folds.out[[K.out]] = i.mix.out[((K.out-1)*d.out+1):n.out]

  Bhat_gblup=list()
  
  CV_errors.out = rep(NA, K.out)
  CV_errors.val.out = rep(list(CV_errors.out), K.out)
  
  CV_corr.out = matrix(NA, nrow=ncol(Y), ncol=ncol(Y))
  CV_corr.val.out = rep(list(CV_corr.out), K.out)
  
  for (l in 1:K.out) {
    cat("outer Fold",l,"\n")
    i.tr.out = unlist(folds.out[-k])
    i.val.out = folds.out[[k]]
    
    X.tr.out = X[i.tr.out,]    
    Y.tr.out = Y[i.tr.out,]   
    
    X.val.out = X[i.val.out,] 
    Y.val.out = Y[i.val.out,] 
    #define the genotype matrix G
    #G=X.tr.out
    X.tr=X.tr.out
    ETA=list(list(X.tr=X.tr, model="BRR"))
    Bhat_gblup[[l]]=matrix(nrow=ncol( X.tr.out), ncol=ncol(Y.tr.out))
   
    for (i in 1:ncol(Y.tr.out)){
    fmR<-BGLR(y=Y.tr.out[,i],ETA=ETA,nIter=500,burnIn=200,thin=10)     
    Bhat_gblup[[l]][,i]<- fmR$ETA[[1]]$b     
    }
   
    CV_errors.val.out[[l]] = mse.matrix(X.val.out%*%Bhat_gblup[[l]], Y.val.out, na.rm = TRUE)
    CV_corr.val.out[[l]] = cor(X.val.out%*%Bhat_gblup[[l]], Y.val.out)
  }
  return(list(B_hat_final=Bhat_gblup, cv_corr=CV_corr.val.out, cv.err=CV_errors.val.out)) 
   
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for cross validation in GBLUP, (marker based) using the R package rrBLUP
BayesL.cv <-function(X,Y, kfold=3){
  set.seed(123)
  n.out = nrow(Y)
  K.out = kfold
  d.out = ceiling(n.out/K.out)
  set.seed(123)
  i.mix.out = sample(1:n.out)
  folds.out = vector(mode="list", length=K.out) 
  for (k in 1:(K.out-1)) {
    folds.out[[k]] = i.mix.out[((k-1)*d.out+1):(k*d.out)]
  }
  folds.out[[K.out]] = i.mix.out[((K.out-1)*d.out+1):n.out]

  Bhat_BayesL=list()
  
  CV_errors.out = rep(NA, K.out)
  CV_errors.val.out = rep(list(CV_errors.out), K.out)
  
  CV_corr.out = matrix(NA, nrow=ncol(Y), ncol=ncol(Y))
  CV_corr.val.out = rep(list(CV_corr.out), K.out)
  
  for (l in 1:K.out) {
    cat("outer Fold",l,"\n")
    i.tr.out = unlist(folds.out[-k])
    i.val.out = folds.out[[k]]
    
    X.tr.out = X[i.tr.out,]   
    Y.tr.out = Y[i.tr.out,]  
    
    X.val.out = X[i.val.out,]
    Y.val.out = Y[i.val.out,]  
    X.tr=X.tr.out
    ETA=list(list(X.tr=X.tr, model="BL"))
    Bhat_BayesL[[l]]=matrix(nrow=ncol( X.tr.out), ncol=ncol(Y.tr.out))
   
    for (i in 1:ncol(Y.tr.out)){
    fmR<-BGLR(y=Y.tr.out[,i],ETA=ETA,nIter=500,burnIn=200,thin=10)     
    Bhat_BayesL[[l]][,i]<- fmR$ETA[[1]]$b     
    }
   
    CV_errors.val.out[[l]] = mse.matrix(X.val.out%*%Bhat_BayesL[[l]], Y.val.out, na.rm = TRUE)
    CV_corr.val.out[[l]] = cor(X.val.out%*%Bhat_BayesL[[l]], Y.val.out)
  }

  return(list(B_hat_final=Bhat_BayesL, cv_corr=CV_corr.val.out, cv.err=CV_errors.val.out)) 
   
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for the L21-joint estimation for B and Omega 
JointEstim<-function(lambdaO, lambdaB, X, Y){ 
  # t=0 initialization
  diag_C_0=c(rep(1, nrow(t(X))))                    # Initialize C as the identity matrix
  invdiag_C_0=1/diag_C_0  
  B_ridge= ridgeEst(lambdaB, X, Y)
  Omega_ridge=t(t(Y)-t(B_ridge)%*%t(X))%*%(t(Y)-t(B_ridge)%*%t(X))+lambdaB*diag(ncol(t(Y)))
  
  Omega_0=solve(Omega_ridge)                               #initialize the inv Cov as omega ridge 
  B_0= round((2/(nrow(t(Y))*lambdaB))*invdiag_C_0*t(X)%*%Omega_0%*%(Y-(2/(nrow(t(Y))*lambdaB))*solve((diag(ncol(t(Y)))+(2/(nrow(t(Y))*lambdaB))*X%*%(invdiag_C_0*t(X)%*%Omega_0)))%*%X%*%(invdiag_C_0*t(X)%*%Omega_0)%*%Y), 8)
  #----------------------
  #update t=1 
  l21B_1=c(rep(0, (nrow(t(X)))))
  for(i in 1:(nrow(t(X)))){
    l21B_1[i]=2*norm((B_0[i,]),type = "2")
  }
  invdiag_C_1=l21B_1
  S_1=(1/nrow(t(Y)))*t(t(Y)-t(B_0)%*%t(X))%*%(t(Y)-t(B_0)%*%t(X))                                      
  Omega_1=glasso(S_1, rho=lambdaO, thr=1.0e-4, maxit=1e4, penalize.diagonal=TRUE, approx=FALSE)$wi    
  B_1= round((2/(nrow(t(Y))*lambdaB))*invdiag_C_1*t(X)%*%Omega_1%*%(Y-(2/(nrow(t(Y))*lambdaB))*solve((diag(ncol(t(Y)))+(2/(nrow(t(Y))*lambdaB))*X%*%(invdiag_C_1*t(X)%*%Omega_1)))%*%X%*%(invdiag_C_1*t(X)%*%Omega_1)%*%Y), 8)
  #-----------------------                   
  #update t=2 
  l21B_2=c(rep(0, (nrow(t(X)))))
  for(i in 1:(nrow(t(X)))){
    l21B_2[i]=2*norm((B_1[i,]),type = "2")
  }
  invdiag_C_2=l21B_2
  S_2=(1/nrow(t(Y)))*t(t(Y)-t(B_1)%*%t(X))%*%(t(Y)-t(B_1)%*%t(X))                                      
  Omega_2=glasso(S_2, rho=lambdaO, thr=1.0e-4, maxit=1e4, penalize.diagonal=TRUE, approx=FALSE)$wi       
  B_2= round((2/(nrow(t(Y))*lambdaB))*invdiag_C_2*t(X)%*%Omega_2%*%(Y-(2/(nrow(t(Y))*lambdaB))*solve((diag(ncol(t(Y)))+(2/(nrow(t(Y))*lambdaB))*X%*%(invdiag_C_2*t(X)%*%Omega_2)))%*%X%*%(invdiag_C_2*t(X)%*%Omega_2)%*%Y), 8)
  #----------------------- 
  #start iteration 
  tcontj=1
  
  while((sum(l21B_2)<=sum(l21B_1))==TRUE  && tcontj<=itermax){  
    l21B_1=l21B_2
    invdiag_C_1=invdiag_C_2
    S_1=S_2
    Omega_1=Omega_2
    B_1=B_2                                             
    
    #----------------------
    l21B_2=c(rep(0, (nrow(t(X)))))
    for(i in 1:(nrow(t(X)))){
      l21B_2[i]=2*norm((B_1[i,]),type = "2")
    }
    invdiag_C_2=l21B_2
    S_2=(1/nrow(t(Y)))*t(t(Y)-t(B_1)%*%t(X))%*%(t(Y)-t(B_1)%*%t(X))                                      
    Omega_2=glasso(S_2, rho=lambdaO, thr=1.0e-4, maxit=1e4, penalize.diagonal=TRUE, approx=FALSE)$wi 
    B_2= round((2/(nrow(t(Y))*lambdaB))*invdiag_C_2*t(X)%*%Omega_2%*%(Y-(2/(nrow(t(Y))*lambdaB))*solve((diag(ncol(t(Y)))+(2/(nrow(t(Y))*lambdaB))*X%*%(invdiag_C_2*t(X)%*%Omega_2)))%*%X%*%(invdiag_C_2*t(X)%*%Omega_2)%*%Y), 8)
    tcontj <- sum(tcontj, 1)                      
    #---------------------- 
    print(tcontj)
    
  }
  return(list(B_1, Omega_1, tcontj))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#function for cross validation in the joint estimation
CV_joint<-function(lambdaO, lambdaB, X, Y, kfold=3) {
  set.seed(123)
  n.out = nrow(Y)
  K.out = kfold
  d.out = ceiling(n.out/K.out)
  set.seed(123)
  i.mix.out = sample(1:n.out)
  folds.out = vector(mode="list", length=K.out) 
  for (k in 1:(K.out-1)) {
    folds.out[[k]] = i.mix.out[((k-1)*d.out+1):(k*d.out)]
  }
  folds.out[[K.out]] = i.mix.out[((K.out-1)*d.out+1):n.out]
  lambdaO = sort(lambdaO)
  lambdaB = sort(lambdaB)

  CV_errors.out = rep(NA, K.out)
  CV_errors.val.out = rep(list(CV_errors.out), K.out)
  
  CV_corr.out = matrix(NA, nrow=ncol(Y), ncol=ncol(Y))
  CV_corr.val.out = rep(list(CV_corr.out), K.out)
  
  opt.i=c()
  opt.j=c()
  opt.lamO.out = c()
  opt.lamB.out = c()
  
  for (l in 1:K.out) {
    cat("outer Fold",l,"\n")
    i.tr.out = unlist(folds.out[-k])
    i.val.out = folds.out[[k]]
    
    X.tr.out = X[i.tr.out,]  
    Y.tr.out = Y[i.tr.out,]   
    
    X.val.out = X[i.val.out,] 
    Y.val.out = Y[i.val.out,] 
    
    #Inner loop for hyper params  
    n.iner = nrow(Y.tr.out)
    K.iner = 3
    d.iner = ceiling(n.iner/K.iner)
    set.seed(123)
    i.mix.iner = sample(1:n.iner)
    folds.iner = vector(mode="list", length=K.iner) 
    for (k in 1:(K.iner-1)) {
      folds.iner[[k]] = i.mix.iner[((k-1)*d.iner+1):(k*d.iner)]
    }
    folds.iner[[K.iner]] = i.mix.iner[((K.iner-1)*d.iner+1):n.iner]  
    
    CV_errors.iner = matrix(NA, nrow=length(lambdaO), ncol=length(lambdaB))
    CV_errors.val.iner = rep(list(CV_errors.iner), K.iner)

    for (k in 1:K.iner) {
      cat("inner Fold",k,"\n")
      
      i.tr.iner = unlist(folds.iner[-k])
      i.val.iner = folds.iner[[k]]
      
      X.tr.iner  = X.tr.out[i.tr.iner,]    
      Y.tr.iner  = Y.tr.out[i.tr.iner,]    
      X.val.iner = X.tr.out[i.val.iner,]   
      Y.val.iner = Y.tr.out[i.val.iner,]  
      
      for(i in 1:length(lambdaO)){
        for(j in 1:length(lambdaB)){
          Estimes.tr.iner = JointEstim(lambdaO[i], lambdaB[j], X.tr.iner, Y.tr.iner) 
          B_hat.tr.iner = Estimes.tr.iner[[1]]
          Omega_hat_tr.iner =  Estimes.tr.iner[[2]]
          CV_errors.val.iner[[k]][i,j] = mean(mse.matrix(X.val.iner%*%B_hat.tr.iner,Y.val.iner))
        }
      }
      
      
      idx.min=which.min((sapply(CV_errors.val.iner, min)))
      opt = which.min(CV_errors.val.iner[[idx.min]]) %% (dim(CV_errors.val.iner[[idx.min]])[1])
      opt = (opt != 0)*opt + (opt == 0)*(dim(CV_errors.val.iner[[idx.min]])[1])
      
      opt.i[k] = opt
      opt.j[k] = which.min(CV_errors.val.iner[[idx.min]][opt,])
      opt.lamO.out[k] = lambdaO[opt.i[k]]
      opt.lamB.out[k] = lambdaB[opt.j[k]]
      
    }
    
        Estimes.tr.out = JointEstim(opt.lamO.out[l], opt.lamB.out[l], X.tr.out, Y.tr.out) 
        B_hat.tr.out = Estimes.tr.out[[1]]
        Omega_hat_tr.out =  Estimes.tr.out[[2]]
        CV_errors.val.out[[l]] = mse.matrix(X.val.out%*%B_hat.tr.out, Y.val.out, na.rm = TRUE)
  CV_corr.val.out[[l]] = cor(X.val.out%*%B_hat.tr.out, Y.val.out)
  }
  
  best.idx.out=which.min((sapply(CV_errors.val.out, min)))
  best.lamO.final = opt.lamO.out[best.idx.out]
  best.lamB.final = opt.lamB.out[best.idx.out]
  Estimes.all = JointEstim(best.lamO.final, best.lamB.final, X, Y) 
  B_hat_final = Estimes.all[[1]]
  Omega_hat_final =  Estimes.all[[2]]
  return(list(best.lamO.final=best.lamO.final, best.lamB.final=best.lamB.final,
              B_hat_final=B_hat_final,Omega_hat_final=Omega_hat_final, cv_corr=CV_corr.val.out, cv.err=CV_errors.val.out)) 
  
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

itermax=2
eps=10^-5  		  # treshold for convergence 

Y<-as.matrix(read.csv("H1xET47_coffeetraits_Step234.csv", head=TRUE, row.names = 1, sep=",")) 
X<-t(as.matrix(read.csv("FinalSNP_H1xET47_impute_mean.csv", head=TRUE, row.names = 1, sep=",")))

Y = scale(Y, scale = TRUE)
X = scale(X, scale = FALSE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 1 ridge regression with CV using glmnet package
reps=20 #number of replicates
CV_errors_gp_H1xET47_step234_RR = list()
CV_corr_gp_H1xET47_step234_RR = list()
RR_gp_H1xET47_step234=list()
for(rr in 1:reps){
  set.seed(rr+300)
  Y=Y[sample(nrow(Y),size=nrow(Y),replace=FALSE),]
  X=X[sample(nrow(X),size=nrow(X),replace=FALSE),]
RR_gp_H1xET47_step234[[rr]] = multi.ridge.cv(X,Y, lambdaB = c(seq(1e-3,0.1,length=10),seq(0.12,2.5,length=10)), kfold = 3, silent = TRUE)
CV_errors_gp_H1xET47_step234_RR[[rr]] = RR_gp_H1xET47_step234[[rr]]$cv.err
CV_corr_gp_H1xET47_step234_RR[[rr]] = RR_gp_H1xET47_step234[[rr]]$cv_corr
}

save(CV_errors_gp_H1xET47_step234_RR, file="CV_errors_gp_H1xET47_step234_RR.RData")
save(CV_corr_gp_H1xET47_step234_RR, file="CV_corr_gp_H1xET47_step234_RR.RData")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 2 multi Lasso Estimate with CV using glmnet package
reps=20#number of replicates
CV_errors_gp_H1xET47_step234_Mlasso = list()
CV_corr_gp_H1xET47_step234_Mlasso = list()
Mlasso_gp_H1xET47_step234=list()
for(rr in 1:reps){
  set.seed(rr+300)
  Y=Y[sample(nrow(Y),size=nrow(Y),replace=FALSE),]
  X=X[sample(nrow(X),size=nrow(X),replace=FALSE),]
Mlasso_gp_H1xET47_step234[[rr]] = multi.lasso.cv(X,Y, lambdaB = c(seq(1e-3,0.1,length=10),seq(0.12,2.5,length=10)), kfold = 3, silent = TRUE)
CV_errors_gp_H1xET47_step234_Mlasso[[rr]] = Mlasso_gp_H1xET47_step234[[rr]]$cv.err
CV_corr_gp_H1xET47_step234_Mlasso[[rr]] = Mlasso_gp_H1xET47_step234[[rr]]$cv_corr
}

save(CV_errors_gp_H1xET47_step234_Mlasso, file="CV_errors_gp_H1xET47_step234_Mlasso.RData")
save(CV_corr_gp_H1xET47_step234_Mlasso, file="CV_corr_gp_H1xET47_step234_Mlasso.RData")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 3 Elastic net Estimate with CV using glmnet package
reps=20
CV_errors_gp_H1xET47_step234_EN = list()
CV_corr_gp_H1xET47_step234_EN = list()
EN_gp_H1xET47_step234=list()
for(rr in 1:reps){
  set.seed(rr+300)
  Y=Y[sample(nrow(Y),size=nrow(Y),replace=FALSE),]
  X=X[sample(nrow(X),size=nrow(X),replace=FALSE),]
 EN_gp_H1xET47_step234[[rr]]  = multi.elasnet.cv(X,Y, lambdaB = seq(1e-4,.2,length=10) ,kfold = 3, silent = TRUE)
 CV_errors_gp_H1xET47_step234_EN[[rr]] = EN_gp_H1xET47_step234[[rr]]$cv.err
 CV_corr_gp_H1xET47_step234_EN[[rr]] = EN_gp_H1xET47_step234[[rr]]$cv_corr
}

save(CV_errors_gp_H1xET47_step234_EN, file="CV_errors_gp_H1xET47_step234_EN.RData")
save(CV_corr_gp_H1xET47_step234_EN, file="CV_corr_gp_H1xET47_step234_EN.RData")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 4 Multitrait Bayes B with BGLR package
reps=20
CV_errors_gp_H1xET47_step234_MBayesB = list()
CV_corr_gp_H1xET47_step234_MBayesB = list()
MBayesB_gp_H1xET47_step234=list()
for(rr in 1:reps){
  set.seed(rr+300)
  Y=Y[sample(nrow(Y),size=nrow(Y),replace=FALSE),]
  X=X[sample(nrow(X),size=nrow(X),replace=FALSE),]
  MBayesB_gp_H1xET47_step234[[rr]] = multi.BayesB.cv(X,Y, kfold=3)
  CV_errors_gp_H1xET47_step234_MBayesB[[rr]] = MBayesB_gp_H1xET47_step234[[rr]]$cv.err
  CV_corr_gp_H1xET47_step234_MBayesB[[rr]] = MBayesB_gp_H1xET47_step234[[rr]]$cv_corr
}

save(CV_errors_gp_H1xET47_step234_MBayesB, file="CV_errors_gp_H1xET47_step234_MBayesB.RData")
save(CV_corr_gp_H1xET47_step234_MBayesB, file="CV_corr_gp_H1xET47_step234_MBayesB.RData")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 5 GBLUP with BGLR package
reps=20

CV_errors_gp_H1xET47_step234_gblup = list()
CV_corr_gp_H1xET47_step234_gblup = list()
gblup_gp_H1xET47_step234=list()
for(rr in 1:reps){
  set.seed(rr+300)
  Y=Y[sample(nrow(Y),size=nrow(Y),replace=FALSE),]
  X=X[sample(nrow(X),size=nrow(X),replace=FALSE),]
  gblup_gp_H1xET47_step234[[rr]] = gblup.cv(X,Y, kfold=3)
  CV_errors_gp_H1xET47_step234_gblup[[rr]] = gblup_gp_H1xET47_step234[[rr]]$cv.err
  CV_corr_gp_H1xET47_step234_gblup[[rr]] = gblup_gp_H1xET47_step234[[rr]]$cv_corr
}

save(CV_errors_gp_H1xET47_step234_gblup, file="CV_errors_gp_H1xET47_step234_gblup.RData")
save(CV_corr_gp_H1xET47_step234_gblup, file="CV_corr_gp_H1xET47_step234_gblup.RData")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 6 BayesL with BGLR package
reps=20
CV_errors_gp_H1xET47_step234_BayesL=list()
CV_corr_gp_H1xET47_step234_BayesL=list()
BayesL_gp_H1xET47_step234=list()
for(rr in 1:reps){
  set.seed(rr+300)
  Y=Y[sample(nrow(Y),size=nrow(Y),replace=FALSE),]
  X=X[sample(nrow(X),size=nrow(X),replace=FALSE),]
  BayesL_gp_H1xET47_step234[[rr]] = BayesL.cv(X,Y, kfold=3)
  CV_errors_gp_H1xET47_step234_BayesL[[rr]] = BayesL_gp_H1xET47_step234[[rr]]$cv.err
  CV_corr_gp_H1xET47_step234_BayesL[[rr]] = BayesL_gp_H1xET47_step234[[rr]]$cv_corr
}

save(CV_errors_gp_H1xET47_step234_BayesL, file="CV_errors_gp_H1xET47_step234_BayesL.RData")
save(CV_corr_gp_H1xET47_step234_BayesL, file="CV_corr_gp_H1xET47_step234_BayesL.RData")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Compare 7 jointl21 with features selection and CV
reps=20
CV_errors_gp_H1xET47_step234_L21joint = list()
CV_corr_gp_H1xET47_step234_L21joint = list()
L21joint_gp_H1xET47_step234=list()
for(rr in 1:reps){
  set.seed(rr+300)
  Y=Y[sample(nrow(Y),size=nrow(Y),replace=FALSE),]
  X=X[sample(nrow(X),size=nrow(X),replace=FALSE),]
  L21joint_gp_H1xET47_step234[[rr]]  = CV_joint(lambdaO=2^seq(-12,-10,1), lambdaB=2^seq(3,5,1), X, Y, kfold=3) 
  CV_errors_gp_H1xET47_step234_L21joint[[rr]] = L21joint_gp_H1xET47_step234[[rr]]$cv.err
  CV_corr_gp_H1xET47_step234_L21joint[[rr]] = L21joint_gp_H1xET47_step234[[rr]]$cv_corr
}

save(CV_errors_gp_H1xET47_step234_L21joint, file="CV_errors_gp_H1xET47_step234_L21joint.RData")
save(CV_corr_gp_H1xET47_step234_L21joint, file="CV_corr_gp_H1xET47_step234_L21joint.RData")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


