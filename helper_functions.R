### Helper functions for the simulation studies performed in "On some limitations of probabilistic models for dimension-reduction:
# illustration in the case of probabilistic formulations of Partial Least Squares" by Etievant and Viallon ###

library(MASS)
library(pracma)
library(matrixcalc)
library(Matrix)
library(permute)
library(methods)
library(robustbase)
library(xtable)
library(parallel)
library(ggplot2)

#####
# Identity - function to create a square matrix with 1s on the diagonal and 0s elsewhere:
#
# Input:
# N - number of rows and columns of the square matrix.
#
# Output:
# M - the square matrix of size N with 1s one the diagonal and 0s elsewhere.
###

Identity = function(N){
  M             = matrix(0, nrow = N, ncol = N)
  for(i in 1:N){
    M[i,i]      = 1
  }
  return(M)
}


#####
# Identity_rec - function to create a matrix with 1s on the diagonal and 0s elsewhere:
#
# Inputs:
# N1 - number of rows of the matrix.
# N2 - number of columns of the matrix.
#
# Output:
# M - the matrix of size N1 x N2 with 1s one the diagonal and 0s elsewhere.
###

Identity_rec = function(N1,N2){
  M             = matrix(0, nrow = N1, ncol = N2)
  for(i in 1:(min(N1,N2))){
    M[i,i]      = 1
  }
  return(M)
}


#####
# Trace - function to compute the Trace of a square matrix, that is to say the sum of its diagonal elements.
#
# Input:
# M - square matrix.
#
# Output:
# tra - the value of the sum of the diagonal elements of M.
###

Trace = function(M){
  if(nrow(M) != ncol(M)){
    print("The matrix is not a square matrix")
    break
  }
  tra = sum(diag(M))
  return(tra)
}


#####
# Q.function - function to compute the expectation of the complete data log likelihood, conditional 
# on the observed data and the values of the parameters.
#
# Inputs:
# sigma_X2, sigma_M2, sigma_U2, Sigma_T, D, A, D - parameters of the model.
# X, M - observed data.
# E_TT and E_T - conditional expectation of T'T and T, with T the first set of latent variables, conditional on the observed data and the values of the parameters.
# E_UU and E_U - conditional expectation of U'U and U, with U the second set of latent variables, conditional on the observed data and the values of the parameters.
# E_UT - conditional expectation of U'T, conditional on the observed data and the values of the parameters.
# Lamba_A - matrix of Lagrangian multipliers used for the semi-orthogonality constraint on A.
# Lamba_B - matrix of Lagrangian multipliers used for the semi-orthogonality constraint on B.
#
# Output :
# q - the value of the expectation of the complete data log likelihood, conditional on the observed data and the values of the parameters.
###

Q.function = function(sigma_X2, sigma_M2, sigma_U2, Sigma_T, A, B, D, E_TT, E_T, E_UU, E_U, E_UT, Lambda_A, Lambda_B, X, M){
  
  n = nrow(X)
  pX = nrow(A)
  pM = nrow(B)
  r = ncol(A)
  
  I_r = Identity(r)
  
  q = - n * pX / 2 * log(2 * pi * sigma_X2) - 1 / (2 * sigma_X2) * Trace(t(X) %*% X - 2 * t(X) %*% E_T %*% t(A) + A %*% E_TT %*% t(A)) - n * pM / 2 * log(2 * pi * sigma_M2) - 1 / (2 * sigma_M2) * Trace(t(M) %*% M - 2 * t(M) %*% E_U %*% t(B) + B %*% E_UU %*% t(B))  - n * r / 2 * log(2 * pi) - n / 2 * log(prod(diag(Sigma_T))) - 1 / 2 * Trace(E_TT %*% solve(Sigma_T)) - n * r / 2 * log(2 * pi * sigma_U2) - 1 / (2 * sigma_U2) * Trace(E_UU - 2 * E_UT %*% D + D %*% E_TT %*% D)
  q = q - 1 / (2 * sigma_X2) * Trace((t(A) %*% A - I_r) %*% Lambda_A) - 1 / (2 * sigma_M2) * Trace((t(B) %*% B - I_r) %*% Lambda_B)
  return(q)
}


#####
# Sigma.function - function to compute the matrix of variance-covariance of X and M according to the values of the parameters.
#
# Inputs :
# sigma_X2, sigma_M2, sigma_U2, Sigma_T, D, A, D - parameters of the model.
#
# Output :
# Sigma - the matrix of variance-covariance of X and M.
###

Sigma.function = function(sigma_X2, sigma_M2, sigma_U2, Sigma_T, A, B, D){
  
  pX = nrow(A)
  pM = nrow(B)
  r = ncol(A)
  
  I_r             = Identity(r)
  I_pX            = Identity(pX)
  I_pM            = Identity(pM)
  
  Sigma       = matrix(NA, nrow = pX + pM, ncol = pX + pM)
  VT          = Sigma_T
  VU          = D^2 %*% VT + sigma_U2 * I_r
  cov_TU      = VT %*% D
  cov_XT      = A %*% VT
  cov_MT      = B %*% cov_TU
  cov_XU      = A %*% cov_TU
  cov_MU      = B %*% VU
  cov_XM      = A %*% cov_TU %*% t(B)
  VX          = A %*% VT %*% t(A) + sigma_X2 * I_pX
  VM          = B %*% VU %*% t(B) + sigma_M2 * I_pM
  Sigma[1:pX, 1:pX]                                 = VX
  Sigma[1:pX, (pX + 1):(pX + pM)]                   = cov_XM
  Sigma[(pX + 1):(pX + pM), (pX + 1):(pX + pM)]     = VM
  Sigma[(pX + 1):(pX + pM), 1:pX]                   = t(cov_XM)
  return(Sigma)
}


#####
# E_step - function to compute the E-step in the EM algorithm.
#
# Inputs:
# sigma_X2_0, sigma_M2_0, sigma_U2_0, Sigma_T_0, D_0, A_0, B_0 - current parameters values.
# X and M - observed data.
#
# Outputs:
# E_TT and E_T - conditional expectation of T'T and T, conditional on the observed data and current values of the parameters.
# E_UU and E_U - conditional expectation of U'U and U, conditional on the observed data and current values of the parameters.
# E_UT - conditional expectation of U'T, conditional on the observed data and current values of the parameters.
# Sigma_inv - inverse of the matrix of variance-covariance of X and M according to the current values of the parameters.
# Lamba_A - matrix of Lagrangian multipliers to be used for the semi-orthogonality constraint on A.
# Lamba_B - matrix of Lagrangian multipliers to be used for the semi-orthogonality constraint on B.
# Mat_A - quantity depending on Lambda_A and E_TT that can be used directly for the A update in the M-step.
# Mat_B - quantity depending on Lambda_B and E_UU that can be used directly for the B update in the M-step.
###

E_step = function(sigma_X2_0, sigma_M2_0, sigma_U2_0, Sigma_T_0, A_0, B_0, D_0, X, M){
  
  n = nrow(X)
  r = ncol(A_0)
  
  I_r           = Identity(r)
  
  Sigma         = Sigma.function(sigma_X2_0, sigma_M2_0, sigma_U2_0, Sigma_T_0, A_0, B_0, D_0)
  Sigma_inv     = solve(Sigma)
  
  VT            = Sigma_T_0
  VU            = D_0^2 %*% VT + sigma_U2_0 * I_r
  cov_TU        = VT %*% D_0
  cov_XT        = A_0 %*% VT
  cov_MT        = B_0 %*% cov_TU
  cov_XU        = A_0 %*% cov_TU
  cov_MU        = B_0 %*% VU
  
  E_T           = as.matrix(cbind(X, M) %*% Sigma_inv %*% rbind(cov_XT, cov_MT))
  vT            = VT - cbind(t(cov_XT), t(cov_MT)) %*% Sigma_inv %*% rbind(cov_XT, cov_MT)
  E_TT          = t(E_T) %*% E_T + n * vT
  
  E_U           = as.matrix(cbind(X, M) %*% Sigma_inv %*% rbind(cov_XU, cov_MU))
  vU            = VU - cbind(t(cov_XU), t(cov_MU)) %*% Sigma_inv %*% rbind(cov_XU, cov_MU)
  E_UU          = t(E_U) %*% E_U + n * vU
  
  c_UT          = cov_TU - cbind(t(cov_XU), t(cov_MU)) %*% Sigma_inv %*% rbind(cov_XT, cov_MT)
  E_UT          = t(E_U) %*% E_T + n * c_UT
  
  EigA          = eigen(t(E_T) %*% X %*% t(X) %*% E_T)
  Mat_A         = EigA$vectors %*% diag(sqrt(EigA$values)) %*% solve(EigA$vectors)
  Lambda_A      = Mat_A - E_TT
  
  EigB          = eigen(t(E_U) %*% M %*% t(M) %*% E_U)
  Mat_B         = EigB$vectors %*% diag(sqrt(EigB$values)) %*% solve(EigB$vectors) 
  Lambda_B      = Mat_B - E_UU
  
  return(list(E_TT = E_TT, E_T = E_T, E_UU = E_UU, E_U = E_U, E_UT = E_UT, Sigma_inv = Sigma_inv, Mat_A = Mat_A, Lambda_A = Lambda_A, Mat_B = Mat_B, Lambda_B = Lambda_B))
}


#####
# M_step - function to compute the M-step in the EM algorithm.
#
# Inputs:
# sigma_X2_0, sigma_M2_0,sigma_U2_0, Sigma_T_0, D_0, A_0, B_0 - current parameters values. 
# X and M - observed data.
# E_TT and E_T - conditional expectation of T'T and T, conditional on the observed data and current values of the parameters.
# E_UU and E_U - conditional expectation of U'U and U, conditional on the observed data and current values of the parameters.
# E_UT - conditional expectation of U'T, conditional on the observed data and current values of the parameters.
# Lamba_A - matrix of Lagrangian multipliers to use for the semi-orthogonality constraint on A.
# Lamba_B - matrix of Lagrangian multipliers to use for the semi-orthogonality constraint on B.
# Mat_A - quantity depending on Lambda_A and E_TT that can be used directly for the A update in the M-step.
# Mat_B - quantity depending on Lambda_B and E_UU that can be used directly for the B update in the M-step.
#
# Outputs:
# sigma_X2_new, sigma_M2_new,sigma_U2_new, Sigma_T_new, D_new, A_new, B_new- updated parameters values .
###

M_step = function(sigma_X2_0, sigma_M2_0, sigma_U2_0, Sigma_T_0, A_0, B_0, D_0, E_TT, E_T, E_UU, E_U, E_UT, Mat_A, Mat_B, Lambda_A, Lambda_B, X, M){
  
  r                       = ncol(A_0)
  n                       = nrow(X)
  pX                      = nrow(A_0)
  pM                      = nrow(B_0)
  
  I_r                    = Identity(r)
  
  A_new = t(X) %*% E_T %*% solve(Mat_A)
  
  sigma_X2_new = (Trace(t(X) %*% X - 2 *(t(X) %*% E_T %*% t(A_new)) + A_new %*% E_TT %*% t(A_new)) + Trace((t(A_new) %*% A_new - I_r) %*% Lambda_A)) / (n * pX)
  
  B_new = t(M) %*% E_U %*% solve(Mat_B)
  
  sigma_M2_new = (Trace(t(M) %*% M - 2 * t(M) %*% E_U %*% t(B_new) + B_new %*% E_UU %*% t(B_new)) + Trace((t(B_new) %*% B_new - I_r) %*% Lambda_B) ) / (n * pM)
  
  Sigma_T_new = hadamard.prod(E_TT, I_r) / n
  
  D_new = diag(diag(E_UT)) %*% solve(diag(diag(E_TT)))
  
  sigma_U2_new = Trace(E_UU - 2 * E_UT %*% D_new + D_new %*% E_TT %*% D_new) / (n * r)
  
  return(list(A_new = A_new, B_new = B_new, D_new = D_new, Sigma_T_new = Sigma_T_new, sigma_X2_new = sigma_X2_new, sigma_M2_new = sigma_M2_new,  sigma_U2_new = sigma_U2_new))
}


#####
# EM_algo - function to run the EM algorithm.
#
# Inputs:
# sigma_X2_0, sigma_M2_0, sigma_U2_0, Sigma_T_0, D_0, A_0, B_0 - initial parameters values. 
# X and Y - observed data.
# (opt) delta.stop - stopping criteria regarding the parameters evolution. Default is 0.05.
# (opt) Lklh.stop - stopping criteria regarding the penalized observed data log-likelihood evolution. Default is 10^(-6).
# (opt) niter.max - maximum number of iterations to be run. Default is 10^4.
#
# Outputs:
# theta, that contains sigma_X2_new, sigma_M2_new, sigma_U2_new, Sigma_T_new, D_new, A_new, B_new - estimated parameters values.
# log_lklh.hat - observed penalized observed data log-likelihood.
# Niter.hat - number of iterations done up to convergence of the algorithm (< niter.max if the stopping criterion are met).
###

EM_algo = function(sigma_X2_0, sigma_M2_0, sigma_U2_0, Sigma_T_0, A_0, B_0, D_0, X, M, delta.stop = NULL, Lklh.stop = NULL, niter.max = NULL){
  
  data            = cbind(X, M)
  n               = nrow(X)
  r               = ncol(A_0)
  
  I_r             = Identity(r)
  
  delta           = Inf
  Lklh            = Inf
  VRAIS           = NULL
  
  if((is.null(delta.stop))||(delta.stop <= 0)){
    delta.stop = 0.05
    print("delta.stop has not been filled out. Default is 0.05. ")
  }
  if((is.null(Lklh.stop))||(Lklh.stop <= 0)){
    Lklh.stop = 10^(-6)
    print("Lklh.stop has not been filled out. Default is 10^(-6).")
  }
  if((is.null(niter.max))||(niter.max <= 0)){
    niter.max = 10^4
    print("niter.max has not been filled out. Default is 10^4 iterations.")
  }
  
  while((delta > delta.stop)||(Lklh > Lklh.stop)){
    
    post = E_step(sigma_X2_0, sigma_M2_0, sigma_U2_0, Sigma_T_0, A_0, B_0, D_0, X, M)
    
    log_lik       = n / 2 * log(det(post$Sigma_inv)) - 1 / 2 * Trace(t(data) %*% data %*% post$Sigma_inv) - 1 / (2 * sigma_X2_0) * Trace((t(A_0) %*% A_0 - I_r) %*% post$Lambda_A) - 1 / (2 * sigma_M2_0) * Trace((t(B_0) %*% B_0 - I_r) %*% post$Lambda_B)
    VRAIS         = cbind(VRAIS, log_lik)
    
    theta = M_step(sigma_X2_0, sigma_M2_0, sigma_U2_0, Sigma_T_0, A_0, B_0, D_0, post$E_TT, post$E_T, post$E_UU, post$E_U, post$E_UT, post$Mat_A, post$Mat_B, post$Lambda_A, post$Lambda_B, X, M)
    
    j = length(VRAIS)
    if(j > 1){
      Lklh = abs((VRAIS[j-1] - VRAIS[j]) / VRAIS[j])
      if(VRAIS[j-1] > VRAIS[j]){
        print("decreasing log-likelihood")
      }
    }
    
    if(j > niter.max){
      break
    }
    
    delta = max(max(abs(A_0 - theta$A_new) / abs(theta$A_new)), max(0,(abs(B_0 - theta$B_new) / abs(theta$B_new))[which(abs(B_0 - theta$B_new) / abs(theta$B_new) != Inf)]), max(abs(diag(D_0 - theta$D_new)) / abs(diag(theta$D_new))), max(abs(diag(Sigma_T_0 - theta$Sigma_T_new)) / abs(diag(theta$Sigma_T_new))), abs(sigma_X2_0 - theta$sigma_X2_new) / abs(theta$sigma_X2_new), abs(sigma_M2_0 - theta$sigma_M2_new) / abs(theta$sigma_M2_new), abs(sigma_U2_0 - theta$sigma_U2_new) / abs(theta$sigma_U2_new))
    
    A_0 = theta$A_new
    B_0 = theta$B_new
    D_0 = theta$D_new
    Sigma_T_0 = theta$Sigma_T_new
    sigma_X2_0 = theta$sigma_X2_new
    sigma_M2_0 = theta$sigma_M2_new
    sigma_U2_0 = theta$sigma_U2_new
  }
  
  theta = order.function(theta$A_new, theta$B_new, theta$D_new, theta$Sigma_T_new, theta$sigma_X2_new, theta$sigma_M2_new, theta$sigma_U2_new)
  
  return(list(theta = theta, log_lklh.hat = VRAIS[length(VRAIS)], Niter.hat = length(VRAIS)))
}


#####
# order.function - function to re-order the estimated parameters values, in accordance with the assumptions made for the
# identifiability of the model.
#
# Inputs:
# sigma_X2_new, sigma_M2_new, sigma_U2_new, Sigma_T_new, D_new, A_new, B_new - estimated parameters values.
#
# Outputs:
# sigma_X2_new, sigma_M2_new, sigma_U2_new, Sigma_T_new, D_new, A_new, B_new - re-ordered estimated parameters values.
###

order.function = function(A_new, B_new, D_new, Sigma_T_new, sigma_X2_new, sigma_M2_new, sigma_U2_new){
  for(i in 1:r){
    if( D_new[i,i] < 0){
      D_new[i,i]    = -  D_new[i,i]
      B_new[,i]     = -  B_new[,i]
    }
  }
  
  or = order(diag( D_new %*%  Sigma_T_new), decreasing = TRUE)
  
  SIGMAT0 =  Sigma_T_new
  for(i in 1:length(or)){
    Sigma_T_new[i,i] = SIGMAT0[or[i],or[i]]
  }
  
  D0 =  D_new
  for(i in 1:length(or)){
    D_new[i,i] = D0[or[i],or[i]]
  }
  
  A0 =  A_new
  for(i in 1:length(or)){
    A_new[,i] = A0[, or[i]]
  }
  
  B0 =  B_new
  for(i in 1:length(or)){
    B_new[,i] = B0[, or[i]]
  }
  
  return(list(A_new = A_new, B_new = B_new, D_new = D_new, Sigma_T_new = Sigma_T_new, sigma_X2_new = sigma_X2_new, sigma_M2_new = sigma_M2_new, sigma_U2_new = sigma_U2_new))
}


#####
# parameters_initialization - function to randomly initialize the parameters.
#
# Inputs:
# (opt) dir_data - directory where to save the initial parameters values.
# nsim, n - informations on the simulation, to name the saved initial parameters values.
# pM - number of variables in M. Also used to name the saved initial parameters values.
# pX - number of variables in X. Also used to name the saved initial parameters values.
# r - number of latents variables T and U, must be < pX and < pM. Also used to name the saved initial parameters values.
# ratio - informations on the simulation, to name the saved initial parameters values.
#
# Outputs:
# Parameters.init that contains sigma_X2_init, sigma_M2_init, sigma_U2_init, Sigma_T_init, D_init, A_init, B_init - initial parameters values that have been randomly generated.
# (opt) file Parameters.init - if the optional input "dir_data" has been filled.
###

parameters_initialization = function(dir_data = NULL, nsim, n, pM, pX, r, ratio){
  
  if(r >= pX){
    print("r must be strictly inferior to pX")
    break
  }
  if(r >= pM){
    print("r must be strictly inferior to pM")
    break
  }
  
  sigma_X2_init   = abs(rnorm(1, sd = 1))
  sigma_M2_init   = abs(rnorm(1, sd = 1))
  sigma_U2_init   = abs(rnorm(1, sd = 1))
  
  I_r             = Identity(r)
  I_pX_r          = Identity_rec(pX,r)
  I_pM_r          = Identity_rec(pM,r)
  
  Sigma_T_0       = matrix(0, nrow = r, ncol = r)
  diago           = abs(rnorm(r, mean = 1, sd = 0.5))
  for(i in 1:r){
    Sigma_T_0[i,i]= diago[i]
  }
  D_0             = matrix(0, nrow = r, ncol = r)
  diago           = abs(rnorm(r, mean = 1, sd = 0.5))
  for(i in 1:r){
    D_0[i,i]      = diago[i]
  }
  
  A_0             = randortho(n = pX, "orthonormal")
  A_init          = A_0[,1:r]
  B_0             = randortho(n = pM, "orthonormal")
  B_init          = B_0[,1:r]
  
  or = order(diag(Sigma_T_0 %*% D_0), decreasing = TRUE)
  DD = D_0
  for(i in 1:length(or)){
    D_0[i,i] = DD[or[i],or[i]]
  }
  SIGMAT = Sigma_T_0
  for(i in 1:length(or)){
    Sigma_T_0[i,i] = SIGMAT[or[i],or[i]]
  }
  Sigma_T_init    = Sigma_T_0
  D_init          = D_0
  
  Parameters.init = list(A_init = A_init, B_init = B_init, D_init = D_init, Sigma_T_init = Sigma_T_init, sigma_X2_init = sigma_X2_init, sigma_M2_init = sigma_M2_init, sigma_U2_init = sigma_U2_init)
  
  if(!is.null(dir_data)){
    
    myfile  = paste0(dir_data, "Parameters_init-n_", n, "-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, "-nsim_", nsim,".RData")
    save(Parameters.init, file = myfile)
  }
  
  return(Parameters.init)
}


#####
# best.permutation.matrix - function to find the best permutation matrix to reorder the columns of an estimated 
# semi-orthogonal matrix. It is done by finding the permutation matrix so that the scalar products of the columns 
# of the true and estimated matrices is the closest to the identity matrix.
#
# Inputs:
# M - true semi-orthogonal matrix.
# M.hat - estimated semi-orthogonal matrix.
#
# Outputs:
# permut.matrix - permutation matrix so that M.hat %*% t(permut.matrix) has its columns ordered in the same way as M.
###

best.permutation.matrix = function(M, M.hat){
  if(nrow(M) != nrow(M.hat)){
    print("The two matrices have to be of the same size")
    break
  }
  if(ncol(M) != ncol(M.hat)){
    print("The two matrices have to be of the same size")
    break
  }
  r = ncol(M)
  
  permut = rbind(c(1:r), allPerms(1:r))
  permut.matrix = NULL
  crit = 0
  for(i in 1:nrow(permut)){
    mat = as.matrix(as(as.integer(permut[i,]), "pMatrix"))
    crit.bis = Trace(abs(t(M) %*% (M.hat) %*% t(mat)))
    if(crit < crit.bis){
      crit = crit.bis
      permut.matrix = mat
    }
  }
  
  return(permut.matrix)
}


#####
# generate_data - function to generate data under the PPLS model proposed by el Bouhaddani et al. (2018), which is:
#          T ~ N(0,Sigma_T)
#          U = T D + e_U,      e_U ~ N(0, sigma_U2 I_r)
#          X = T A' + e_X,     e_X ~ N(0, sigma_X2 I_pX)
#          M = U B' + e_M,     e_M ~ N(0, sigma_M2 I_pM)
# with 
# A and B semi-orthogonal matrices.
# Sigma_T and D diagonal matrices with strictly positive diagonal elements.
#
# Inputs:
# (opt) dir_data - directory where to save the initial parameters values.
# nsim - information on the simulation, to name the saved the generated data.
# n - sample size. Also used to name the saved the generated data.
# pM - number of variables in M. Also used to name the saved the generated data.
# pX - number of variables in X. Also used to name the saved the generated data.
# r - number of latents variables T and U, must be < pX and < pM. Also used to name the saved the generated data.
# ratio - Signal to Noise Ratio. Also used to name the saved initial parameters values.
# (opt) theta.true - parameters values to be used for data generation. Default values are generated randomly otherwise.
#
# Outputs:
# Data - contains the generated data X, M, T and U, as well as parameters values used for the generation. Also contains n, pX, pM, r, ratio and nsim.
### 

generate_data = function(dir_data = NULL, nsim, n, pX, pM, r, ratio, theta.true = NULL){
  
  if(r >= pX){
    print("r must be strictly inferior to pX")
    break
  }
  if(r >= pM){
    print("r must be strictly inferior to pM")
    break
  }
  
  #############################
  ##### Parameters values #####
  A             = theta.true$A
  B             = theta.true$B
  diag_D        = diag(theta.true$D)
  diag_Sigma_T  = diag(theta.true$Sigma_T)
  Sigma_T       = theta.true$Sigma_T
  D             = theta.true$D
  sigma_X2      = theta.true$sigma_X2
  sigma_M2      = theta.true$sigma_M2
  sigma_U2      = theta.true$sigma_U2
  
  I_r             = Identity(r)
  I_pX            = Identity(pX)
  I_pM            = Identity(pM)
  I_pX_r          = Identity_rec(pX,r)
  I_pM_r          = Identity_rec(pM,r)
  
  if((is.null(diag_Sigma_T))|(length(diag_Sigma_T) != r)){
    Sigma_T       = matrix(0, nrow = r, ncol = r)
    diago         = abs(rnorm(r, mean = 3, sd = 6))
    for(i in 1:r){
      Sigma_T[i,i]= diago[i]
    }
  }else{
    Sigma_T       = diag(diag_Sigma_T)
  }
  if((is.null(diag_D))|(length(diag_D) != r)){
    D             = matrix(0, nrow = r, ncol = r)
    diago         = abs(rnorm(r, mean = 3, sd = 6))
    for(i in 1:r){
      D[i,i]      = diago[i]
    }
  }else{
    D             = diag(diag_D)
  }
  if((is.null(A))||(nrow(A) != pX)||(ncol(A) != r)||(sum(abs(t(A)%*%A - Identity(r))) > 10^(-5))){
    A               = randortho(n = pX, "orthonormal")
    A               = A[,1:r]
  }
  if((is.null(B))||(nrow(B) != pM)||(ncol(B) != r)||(sum(abs(t(B)%*%B - Identity(r))) > 10^(-5))){
    B               = randortho(n = pM, "orthonormal")
    B               = B[,1:r]
  }
  if((is.null(sigma_U2))||(length(sigma_U2) != 1)||(sigma_U2 <= 0)){
    sigma_U2      = round(Trace(D^2 %*% Sigma_T)) / (ratio * r)
  }
  if((is.null(sigma_X2))||(length(sigma_X2) != 1)||(sigma_X2 <= 0)){
    sigma_X2      = round(Trace(A %*% Sigma_T %*% t(A))) / (ratio * pX)
  }
  if((is.null(sigma_M2))||(length(sigma_M2) != 1)||(sigma_M2 <= 0)){
    sigma_M2      = round(Trace(B %*% (D^2 %*% Sigma_T + sigma_U2 * I_r) %*% t(B))) / (ratio * pM)
  }
  
  #####################################
  ######## Data set simulation ########
  T               = mvrnorm(n = n, mu = rep(0, r), Sigma_T, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  
  eU              = mvrnorm(n = n, mu = rep(0, r), sigma_U2 * I_r , tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  U               = T %*% D + eU
  
  eX              = mvrnorm(n = n, mu = rep(0, pX), sigma_X2 * I_pX , tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  
  eM              = mvrnorm(n = n, mu = rep(0, pM), sigma_M2 * I_pM , tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  
  X               = T %*% t(A) + eX
  M               = U %*% t(B) + eM
  
  data            = cbind(X, M)
  
  Data = list(n = n, pX = pX, pM = pM, r = r, ratio = ratio, X = X, M = M, T = T, U = U, A = A, B = B, D = D, Sigma_T = Sigma_T, sigma_X2 = sigma_X2, sigma_M2 = sigma_M2, sigma_U2 = sigma_U2, nsim = nsim)
  
  if(!is.null(dir_data)){
    myfile  = paste0(dir_data,"sample_n_" , n, "-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, "-nsim_", nsim)
    myfile  = paste0(myfile, ".RData")
    save(Data, file = myfile)
  }else{return(Data = Data)}
}


#####
# generate_data_misspecified - function to generate data under the following misspecified PPLS model
#          T ~ N(0,Sigma_T)
#          U = T D + e_U,      e_U ~ N(0, sigma_U2 I_r)
#          X = T A' + e_X,     e_X ~ N(0, Psi_X)
#          M = U B' + e_M,     e_M ~ N(0, Psi_M)
# with 
# Psi_X and Psi_M arbitrary positive semi-definite matrices.
# A and B semi-orthogonal matrices.
# Sigma_T and D diagonal matrices with strictly positive diagonal elements.
#
# Inputs:
# (opt) dir_data - directory where to save the initial parameters values.
# nsim - information on the simulation, to name the saved the generated data.
# n - sample size. Also used to name the saved the generated data.
# pM - number of variables in M. Also used to name the saved the generated data.
# pX - number of variables in X. Also used to name the saved the generated data.
# r - number of latents variables T and U, must be < pX and < pM. Also used to name the saved the generated data.
# ratio - Signal to Noise Ratio. Also used to name the saved initial parameters values.
# (opt) theta.true - parameters values to be used for data generation. Default values are generated randomly otherwise.
#
# Outputs:
# Data - contains the generated data X, M, T and U, as well as parameters values used for the generation. Also contains n, pX, pM, r, ratio and nsim.
### 

generate_data_misspecified = function(dir_data = NULL, nsim, n, pX, pM, r, ratio, theta.true = NULL){
  
  if(r >= pX){
    print("r must be strictly inferior to pX")
    break
  }
  if(r >= pM){
    print("r must be strictly inferior to pM")
    break
  }
  
  #############################
  ##### Parameters values #####
  A             = theta.true$A
  B             = theta.true$B
  diag_D        = diag(theta.true$D)
  diag_Sigma_T  = diag(theta.true$Sigma_T)
  Sigma_T       = theta.true$Sigma_T
  D             = theta.true$D
  Sigma_X      = theta.true$Sigma_X
  Sigma_M      = theta.true$Sigma_M
  sigma_U2      = theta.true$sigma_U2
  
  I_r             = Identity(r)
  I_pX            = Identity(pX)
  I_pM            = Identity(pM)
  I_pX_r          = Identity_rec(pX,r)
  I_pM_r          = Identity_rec(pM,r)
  
  
  if((is.null(diag_Sigma_T))|(length(diag_Sigma_T) != r)){
    Sigma_T       = matrix(0, nrow = r, ncol = r)
    diago         = abs(rnorm(r, mean = 3, sd = 6))
    for(i in 1:r){
      Sigma_T[i,i]= diago[i]
    }
  }else{
    Sigma_T       = diag(diag_Sigma_T)
  }
  if((is.null(diag_D))|(length(diag_D) != r)){
    D             = matrix(0, nrow = r, ncol = r)
    diago         = abs(rnorm(r, mean = 3, sd = 6))
    for(i in 1:r){
      D[i,i]      = diago[i]
    }
  }else{
    D             = diag(diag_D)
  }
  if((is.null(A))||(nrow(A) != pX)||(ncol(A) != r)||(sum(abs(t(A)%*%A - Identity(r))) > 10^(-5))){
    A               = randortho(n = pX, "orthonormal")
    A               = A[,1:r]
  }
  if((is.null(B))||(nrow(B) != pM)||(ncol(B) != r)||(sum(abs(t(B)%*%B - Identity(r))) > 10^(-5))){
    B               = randortho(n = pM, "orthonormal")
    B               = B[,1:r]
  }
  if((is.null(sigma_U2))||(length(sigma_U2) != 1)||(sigma_U2 <= 0)){
    sigma_U2      = round(Trace(D^2 %*% Sigma_T)) / (ratio * r)
  }
  if((is.null(Sigma_X))||(nrow(Sigma_X) != pX)){
    Mat = matrix(rnorm(pX*pX), nrow = pX, ncol = pX)
    Sigma_X = Mat %*% t(Mat)
    Sigma_X = Sigma_X * Trace(Sigma_T) / Trace(Sigma_X) / ratio
  }
  if((is.null(Sigma_M))||(nrow(Sigma_M) != pM)){
    Mat = matrix(rnorm(pM*pM), nrow = pM, ncol = pM)
    Sigma_M = Mat %*% t(Mat)
    Sigma_M = Sigma_M * Trace((D^2 %*% Sigma_T + sigma_U2 * I_r)) / Trace(Sigma_M) / ratio
  }
  
  
  #####################################
  ######## Data set simulation ########
  T               = mvrnorm(n = n, mu = rep(0, r), Sigma_T, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  
  eU              = mvrnorm(n = n, mu = rep(0, r), sigma_U2 * I_r , tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  U               = T %*% D + eU
  
  eX              = mvrnorm(n = n, mu = rep(0, pX), Sigma_X, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  
  eM              = mvrnorm(n = n, mu = rep(0, pM), Sigma_M , tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  
  X               = T %*% t(A) + eX
  M               = U %*% t(B) + eM
  
  data            = cbind(X, M)
  
  Data = list(n = n, pX = pX, pM = pM, r = r, ratio = ratio, X = X, M = M, T = T, U = U, A = A, B = B, D = D, Sigma_T = Sigma_T, Sigma_X = Sigma_X, Sigma_M = Sigma_M, sigma_U2 = sigma_U2, nsim = nsim)
  
  if(!is.null(dir_data)){
    myfile  = paste0(dir_data,"sample_misspecified_n_" , n, "-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, "-nsim_", nsim)
    myfile  = paste0(myfile, ".RData")
    save(Data, file = myfile)
  }else{return(Data = Data)}
  
}


#####
# Comparision_weights - function to compare the true weight matrices A and B with 
#  the estimates returned by the EM alogorithm devised under the PPLS model proposed by el Bouhaddani et al. (2018).
#  the estimates obtained via PLS-SVD on (X,M),
#  the estimates obtained via PLS-W2A on (X,M),
#  the estimates obtained via two distinct PCAs on X and M
# and to compare the estimates returned by the EM alogorithm (devised under the PPLS model proposed by el Bouhaddani et al. (2018))
# with
#  the estimates obtained via PLS-SVD on (X,M),
#  the estimates obtained via PLS-W2A on (X,M),
#  the estimates obtained via two distinct PCAs on X and M,
# under two simulation studies, in our case when X and M have been generated under
#  the PPLS model proposed by el Bouhaddani et al. (2018):
#          T ~ N(0,Sigma_T)
#          U = T D + e_U,      e_U ~ N(0, sigma_U2 I_r)
#          X = T A' + e_X,     e_X ~ N(0, sigma_X2 I_pX)
#          M = U B' + e_M,     e_M ~ N(0, sigma_M2 I_pM)
# with 
# A and B semi-orthogonal matrices.
# Sigma_T and D diagonal matrices with strictly positive diagonal elements.
#
#  the following misspecified PPLS model:
#          T ~ N(0,Sigma_T)
#          U = T D + e_U,      e_U ~ N(0, sigma_U2 I_r)
#          X = T A' + e_X,     e_X ~ N(0, Psi_X)
#          M = U B' + e_M,     e_M ~ N(0, Psi_M)
# with 
# Psi_X and Psi_M arbitrary positive semi-definite matrices.
# A and B semi-orthogonal matrices.
# Sigma_T and D diagonal matrices with strictly positive diagonal elements.
#
# Inputs:
# dir_data_sim1 - directory where true parameters values, generated data and parameters estimates have been saved for the first simulation study (in our case, when the PPLS model proposed by el Bouhaddani et al. (2018) is correctly specified)
# dir_data_sim2 - directory where true parameters values, generated data and parameters estimates have been saved for the second simulation study (in our case, when the PPLS model proposed by el Bouhaddani et al. (2018) is misspecified)
# Nsim - number of replicates chosen for the two simulation studies.
# pX - number of variables in X chosen for the two simulation studies.
# pM - number of variables in M chosen for the two simulation studies.
# r - number of latents variables T and U chosen for the two simulation studies.
# N - sample size(s) chosen for the two simulation studies. Can be a single value or a vector of values. 
# ratio - Signal to Noise Ratio chosen for the two simulation studies.
#
# Outputs:
# res2 - dataframe which contains the results of all the comparisons, for the two simulation studies.
### 

Comparision_weights = function(dir_data_sim1, dir_data_sim2, Nsim, pX, pM, r, N, ratio){
  I_pX = Identity(pX)
  I_pM = Identity(pM)
  
  dir_data = dir_data_sim2
  setwd(dir_data)
  
  RES = NULL
  COL_A_PPLS2 = NULL
  COL_B_PPLS2 = NULL
  COL_A_PCA2 = NULL
  COL_B_PCA2 = NULL
  COL_A_PLS2 = NULL
  COL_B_PLS2 = NULL
  COL_A_PLSW2A2 = NULL
  COL_B_PLSW2A2 = NULL
  RES2 = NULL
  
  for(i in 1:r){
    for(n in N){
      
      COL_A_true_PPLS = NULL
      COL_B_true_PPLS = NULL
      COL_A_true_PCA = NULL
      COL_B_true_PCA = NULL
      COL_A_true_PLS = NULL
      COL_B_true_PLS = NULL
      COL_A_PPLS_PCA = NULL
      COL_B_PPLS_PCA = NULL
      COL_A_PPLS_PLS = NULL
      COL_B_PPLS_PLS = NULL
      COL_A_true_PLSW2A = NULL
      COL_B_true_PLSW2A = NULL
      COL_A_PPLS_PLSW2A = NULL
      COL_B_PPLS_PLSW2A = NULL
      
      load(paste0(dir_data, "Parameters_misspecified-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, ".RData"))
      
      for(nsim in 1:Nsim){
        
        load(paste0(dir_data,"PPLS_misspecified_estimation_n_" , n, "-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, "-nsim_", nsim, ".RData"))
        load(paste0(dir_data,"sample_misspecified_n_" , n, "-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, "-nsim_", nsim, ".RData"))
        
        X = Data$X
        M = Data$M
        EX = eigen(var(X))
        EM = eigen(var(M))
        S = svd(cov(X,M))
        SXM = svd(cov(X,M))
        
        T1 = X %*% SXM$u[,1]
        U1 = M %*% SXM$v[,1]
        
        A1 = SXM$u[,1]
        B1 = SXM$v[,1]
        
        regX1 = lm(X~T1)
        regM1 = lm(M~U1)
        
        X1 = regX1$residuals
        M1 = regM1$residuals
        
        SXM1 = svd(cov(X1,M1))
        
        T2 = X1 %*% SXM1$u[,1]
        U2 = M1 %*% SXM1$v[,1]
        
        A2 = (I_pX - A1 %*% solve(t(A1) %*% t(X) %*% X %*% A1) %*% t(A1) %*% t(X) %*% X) %*% SXM1$u[,1]
        B2 = (I_pM - B1 %*% solve(t(B1) %*% t(M) %*% M %*% B1) %*% t(B1) %*% t(M) %*% M) %*% SXM1$v[,1]
        
        regX2 = lm(X1~T2)
        regM2 = lm(M1~U2)
        
        X2 = regX2$residuals
        M2 = regM2$residuals
        
        SXM2 = svd(cov(X2,M2))
        
        T3 = X2 %*% SXM2$u[,1]
        U3 = M2 %*% SXM2$v[,1]
        
        A3 = (I_pX - A1 %*% solve(t(A1) %*% t(X) %*% X %*% A1) %*% t(A1) %*% t(X) %*% X) %*% (I_pX - A2 %*% solve(t(A2) %*% t(X1) %*% X1%*% A2) %*% t(A2) %*% t(X1) %*% X1) %*% SXM2$u[,1]
        B3 = (I_pM - B1 %*% solve(t(B1) %*% t(M) %*% M %*% B1) %*% t(B1) %*% t(M) %*% M) %*% (I_pM - B2 %*% solve(t(B2) %*% t(M1) %*% M1%*% B2) %*% t(B2) %*% t(M1) %*% M1) %*% SXM2$v[,1]
        
        A_PLSW2A = cbind(A1, A2, A3)
        B_PLSW2A = cbind(B1, B2, B3)
        
        A_PLSW2A = A_PLSW2A / matrix(sqrt(colSums((A_PLSW2A)^2)), nrow = pX, ncol = r, byrow = TRUE)
        B_PLSW2A = B_PLSW2A / matrix(sqrt(colSums((B_PLSW2A)^2)), nrow = pM, ncol = r, byrow = TRUE)
        
        PA_true_PPLS = best.permutation.matrix(Data$A, theta.hat_PPLS$A_new)
        PB_true_PPLS = best.permutation.matrix(Data$B, theta.hat_PPLS$B_new)
        
        PA_true_PCA = best.permutation.matrix(Data$A, EX$vectors[,1:r])
        PB_true_PCA = best.permutation.matrix(Data$B, EM$vectors[,1:r])
        
        PA_true_PLS = best.permutation.matrix(Data$A, S$u[,1:r])
        PB_true_PLS = best.permutation.matrix(Data$B, S$v[,1:r])
        
        PA_true_PLSW2A = best.permutation.matrix(Data$A, A_PLSW2A)
        PB_true_PLSW2A = best.permutation.matrix(Data$B, B_PLSW2A)
        
        name = paste("colA_true_PPLS", i, sep = "")
        assign(name, abs(diag(t(Data$A) %*% theta.hat_PPLS$A_new %*% t(PA_true_PPLS))[i]))
        COL_A_true_PPLS = rbind(COL_A_true_PPLS, get(name))
        
        name = paste("colB_true_PPLS", i, sep = "")
        assign(name, abs(diag(t(Data$B) %*% theta.hat_PPLS$B_new %*% t(PB_true_PPLS))[i]))
        COL_B_true_PPLS = rbind(COL_B_true_PPLS, get(name))
        
        name = paste("colA_true_PCA", i, sep = "")
        assign(name, abs(diag(t(Data$A) %*% EX$vectors[,1:r] %*% t(PA_true_PCA))[i]))
        COL_A_true_PCA = rbind(COL_A_true_PCA, get(name))
        
        name = paste("colB_true_PCA", i, sep = "")
        assign(name, abs(diag(t(Data$B) %*% EM$vectors[,1:r] %*% t(PB_true_PCA))[i]))
        COL_B_true_PCA = rbind(COL_B_true_PCA, get(name))
        
        name = paste("colA_true_PLS", i, sep = "")
        assign(name, abs(diag(t(Data$A) %*% S$u[,1:r] %*% t(PA_true_PLS))[i]))
        COL_A_true_PLS = rbind(COL_A_true_PLS, get(name))
        
        name = paste("colB_true_PLS", i, sep = "")
        assign(name, abs(diag(t(Data$B) %*% S$v[,1:r] %*% t(PB_true_PLS))[i]))
        COL_B_true_PLS = rbind(COL_B_true_PLS, get(name))
        
        name = paste("colA_true_PLSW2A", i, sep = "")
        assign(name, abs(diag(t(Data$A) %*% A_PLSW2A %*% t(PA_true_PLSW2A))[i]))
        COL_A_true_PLSW2A = rbind(COL_A_true_PLSW2A, get(name))
        
        name = paste("colB_true_PLSW2A", i, sep = "")
        assign(name, abs(diag(t(Data$B) %*% B_PLSW2A %*% t(PB_true_PLSW2A))[i]))
        COL_B_true_PLSW2A = rbind(COL_B_true_PLSW2A, get(name))
        
        name = paste("colA_PPLS_PCA", i, sep = "")
        assign(name, abs(diag(PA_true_PCA %*% t(EX$vectors[,1:r]) %*% theta.hat_PPLS$A_new %*% t(PA_true_PPLS))[i]))
        COL_A_PPLS_PCA = rbind(COL_A_PPLS_PCA, get(name))
        
        name = paste("colB_PPLS_PCA", i, sep = "")
        assign(name, abs(diag(PB_true_PCA %*% t(EM$vectors[,1:r]) %*% theta.hat_PPLS$B_new %*% t(PB_true_PPLS))[i]))
        COL_B_PPLS_PCA = rbind(COL_B_PPLS_PCA, get(name))
        
        name = paste("colA_PPLS_PLS", i, sep = "")
        assign(name, abs(diag(PA_true_PLS %*% t(S$u[,1:r]) %*% theta.hat_PPLS$A_new %*% t(PA_true_PPLS))[i]))
        COL_A_PPLS_PLS = rbind(COL_A_PPLS_PLS, get(name))
        
        name = paste("colB_PPLS_PLS", i, sep = "")
        assign(name, abs(diag(PB_true_PLS %*% t(S$v[,1:r]) %*% theta.hat_PPLS$B_new %*% t(PB_true_PPLS))[i]))
        COL_B_PPLS_PLS = rbind(COL_B_PPLS_PLS, get(name))   
        
        name = paste("colA_PPLS_PLSW2A", i, sep = "")
        assign(name, abs(diag(PA_true_PLSW2A %*% t(A_PLSW2A) %*% theta.hat_PPLS$A_new %*% t(PA_true_PPLS))[i]))
        COL_A_PPLS_PLSW2A = rbind(COL_A_PPLS_PLSW2A, get(name))
        
        name = paste("colB_PPLS_PLSW2A", i, sep = "")
        assign(name, abs(diag(PB_true_PLSW2A %*% t(B_PLSW2A) %*% theta.hat_PPLS$B_new %*% t(PB_true_PPLS))[i]))
        COL_B_PPLS_PLSW2A = rbind(COL_B_PPLS_PLSW2A, get(name))   
        
        col = cbind(COL_A_true_PPLS, COL_B_true_PPLS, COL_A_true_PCA, COL_B_true_PCA, COL_A_true_PLS, COL_B_true_PLS, COL_A_true_PLSW2A, COL_B_true_PLSW2A, COL_A_PPLS_PCA, COL_B_PPLS_PCA, COL_A_PPLS_PLS, COL_B_PPLS_PLS, COL_A_PPLS_PLSW2A, COL_B_PPLS_PLSW2A)
      }
      
      sd_A_true_PPLS2 = sd(COL_A_true_PPLS)
      sd_B_true_PPLS2 = sd(COL_B_true_PPLS)
      sd_A_true_PCA2 = sd(COL_A_true_PCA)
      sd_B_true_PCA2 = sd(COL_B_true_PCA)
      sd_A_true_PLS2 = sd(COL_A_true_PLS)
      sd_B_true_PLS2 = sd(COL_B_true_PLS)
      sd_A_true_PLSW2A2 = sd(COL_A_true_PLSW2A)
      sd_B_true_PLSW2A2 = sd(COL_B_true_PLSW2A)
      sd_A_PPLS_PCA2 = sd(COL_A_PPLS_PCA)
      sd_B_PPLS_PCA2 = sd(COL_B_PPLS_PCA)
      sd_A_PPLS_PLS2 = sd(COL_A_PPLS_PLS)
      sd_B_PPLS_PLS2 = sd(COL_B_PPLS_PLS)
      sd_A_PPLS_PLSW2A2 = sd(COL_A_PPLS_PLSW2A)
      sd_B_PPLS_PLSW2A2 = sd(COL_B_PPLS_PLSW2A)
      
      median_A_true_PPLS2 = median(COL_A_true_PPLS)
      median_B_true_PPLS2 = median(COL_B_true_PPLS)
      median_A_true_PCA2 = median(COL_A_true_PCA)
      median_B_true_PCA2 = median(COL_B_true_PCA)
      median_A_true_PLS2 = median(COL_A_true_PLS)
      median_B_true_PLS2 = median(COL_B_true_PLS)
      median_A_true_PLSW2A2 = median(COL_A_true_PLSW2A)
      median_B_true_PLSW2A2 = median(COL_B_true_PLSW2A)
      median_A_PPLS_PCA2 = median(COL_A_PPLS_PCA)
      median_B_PPLS_PCA2 = median(COL_B_PPLS_PCA)
      median_A_PPLS_PLS2 = median(COL_A_PPLS_PLS)
      median_B_PPLS_PLS2 = median(COL_B_PPLS_PLS)
      median_A_PPLS_PLSW2A2 = median(COL_A_PPLS_PLSW2A)
      median_B_PPLS_PLSW2A2 = median(COL_B_PPLS_PLSW2A)
      
      COL_A_PPLS2 = rbind(COL_A_PPLS2, cbind(n = n, column = i, "misspecified", "true", "COL_A_PPLS", "A", "PPLS", COL_A_true_PPLS, sd_A_true_PPLS2, median_A_true_PPLS2))
      COL_B_PPLS2 = rbind(COL_B_PPLS2, cbind(n = n, column = i, "misspecified", "true", "COL_B_PPLS", "B", "PPLS", COL_B_true_PPLS, sd_B_true_PPLS2, median_B_true_PPLS2))
      COL_A_PCA2 = rbind(COL_A_PCA2, cbind(n = n, column = i, "misspecified","true",  "COL_A_PCA", "A", "PCA", COL_A_true_PCA, sd_A_true_PCA2, median_A_true_PCA2))
      COL_B_PCA2 = rbind(COL_B_PCA2, cbind(n = n, column = i, "misspecified", "true", "COL_B_PCA", "B", "PCA", COL_B_true_PCA, sd_B_true_PCA2, median_B_true_PCA2))
      COL_A_PLS2 = rbind(COL_A_PLS2, cbind(n = n, column = i, "misspecified", "true", "COL_A_PLS", "A", "PLS-SVD", COL_A_true_PLS, sd_A_true_PLS2, median_A_true_PLS2))
      COL_B_PLS2 = rbind(COL_B_PLS2, cbind(n = n, column = i, "misspecified", "true", "COL_B_PLS", "B", "PLS-SVD", COL_B_true_PLS, sd_B_true_PLS2, median_B_true_PLS2))
      COL_A_PLSW2A2 = rbind(COL_A_PLSW2A2, cbind(n = n, column = i, "misspecified", "true", "COL_A_PLSW2A", "A", "PLS-W2A", COL_A_true_PLSW2A, sd_A_true_PLSW2A2, median_A_true_PLSW2A2))
      COL_B_PLSW2A2 = rbind(COL_B_PLSW2A2, cbind(n = n, column = i, "misspecified", "true", "COL_B_PLSW2A", "B", "PLS-W2A", COL_B_true_PLSW2A, sd_B_true_PLSW2A2, median_B_true_PLSW2A2))
      COL_A_PCA2 = rbind(COL_A_PCA2, cbind(n = n, column = i, "misspecified", "PPLS", "COL_A_PCA", "A", "PCA", COL_A_PPLS_PCA, sd_A_PPLS_PCA2, median_A_PPLS_PCA2))
      COL_B_PCA2 = rbind(COL_B_PCA2, cbind(n = n, column = i, "misspecified", "PPLS", "COL_B_PCA", "B", "PCA", COL_B_PPLS_PCA, sd_B_PPLS_PCA2, median_B_PPLS_PCA2))
      COL_A_PLS2 = rbind(COL_A_PLS2, cbind(n = n, column = i, "misspecified", "PPLS", "COL_A_PLS", "A", "PLS-SVD", COL_A_PPLS_PLS, sd_A_PPLS_PLS2, median_A_PPLS_PLS2))
      COL_B_PLS2 = rbind(COL_B_PLS2, cbind(n = n, column = i, "misspecified", "PPLS", "COL_B_PLS", "B", "PLS-SVD", COL_B_PPLS_PLS, sd_B_PPLS_PLS2, median_B_PPLS_PLS2))
      COL_A_PLSW2A2 = rbind(COL_A_PLSW2A2, cbind(n = n, column = i, "misspecified", "PPLS", "COL_A_PLSW2A", "A", "PLS-W2A", COL_A_PPLS_PLSW2A, sd_A_PPLS_PLSW2A2, median_A_PPLS_PLSW2A2))
      COL_B_PLSW2A2 = rbind(COL_B_PLSW2A2, cbind(n = n, column = i, "misspecified", "PPLS", "COL_B_PLSW2A", "B", "PLS-W2A", COL_B_PPLS_PLSW2A, sd_B_PPLS_PLSW2A2, median_B_PPLS_PLSW2A2))
      
      RES = rbind(RES, cbind(n = n, column = i, ratio = ratio, "misspecified", col))
      
    }
  }
  
  dir_data        = dir_data_sim1
  setwd(dir_data)
  
  for(i in 1:r){
    for(n in N){
      
      COL_A_true_PPLS = NULL
      COL_B_true_PPLS = NULL
      COL_A_true_PCA = NULL
      COL_B_true_PCA = NULL
      COL_A_true_PLS = NULL
      COL_B_true_PLS = NULL
      COL_A_PPLS_PCA = NULL
      COL_B_PPLS_PCA = NULL
      COL_A_PPLS_PLS = NULL
      COL_B_PPLS_PLS = NULL
      COL_A_true_PLSW2A = NULL
      COL_B_true_PLSW2A = NULL
      COL_A_PPLS_PLSW2A = NULL
      COL_B_PPLS_PLSW2A = NULL
      
      load(paste0(dir_data, "Parameters-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, ".RData"))
      
      for(nsim in 1:Nsim){
        
        load(paste0(dir_data,"PPLS_estimation_n_" , n, "-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, "-nsim_", nsim, ".RData"))
        load(paste0(dir_data,"sample_n_" , n, "-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, "-nsim_", nsim, ".RData"))
        
        X = Data$X
        M = Data$M
        
        EX = eigen(var(X))
        EM = eigen(var(M))
        S = svd(cov(X,M))
        
        SXM = svd(cov(X,M))
        
        T1 = X %*% SXM$u[,1]
        U1 = M %*% SXM$v[,1]
        
        A1 = SXM$u[,1]
        B1 = SXM$v[,1]
        
        regX1 = lm(X~T1)
        regM1 = lm(M~U1)
        
        X1 = regX1$residuals
        M1 = regM1$residuals
        
        SXM1 = svd(cov(X1,M1))
        
        T2 = X1 %*% SXM1$u[,1]
        U2 = M1 %*% SXM1$v[,1]
        
        A2 = (I_pX - A1 %*% solve(t(A1) %*% t(X) %*% X %*% A1) %*% t(A1) %*% t(X) %*% X) %*% SXM1$u[,1]
        B2 = (I_pM - B1 %*% solve(t(B1) %*% t(M) %*% M %*% B1) %*% t(B1) %*% t(M) %*% M) %*% SXM1$v[,1]
        
        regX2 = lm(X1~T2)
        regM2 = lm(M1~U2)
        
        X2 = regX2$residuals
        M2 = regM2$residuals
        
        SXM2 = svd(cov(X2,M2))
        
        T3 = X2 %*% SXM2$u[,1]
        U3 = M2 %*% SXM2$v[,1]
        
        A3 = (I_pX - A1 %*% solve(t(A1) %*% t(X) %*% X %*% A1) %*% t(A1) %*% t(X) %*% X) %*% (I_pX - A2 %*% solve(t(A2) %*% t(X1) %*% X1%*% A2) %*% t(A2) %*% t(X1) %*% X1) %*% SXM2$u[,1]
        B3 = (I_pM - B1 %*% solve(t(B1) %*% t(M) %*% M %*% B1) %*% t(B1) %*% t(M) %*% M) %*% (I_pM - B2 %*% solve(t(B2) %*% t(M1) %*% M1%*% B2) %*% t(B2) %*% t(M1) %*% M1) %*% SXM2$v[,1]
        
        A_PLSW2A = cbind(A1, A2, A3)
        B_PLSW2A = cbind(B1, B2, B3)
        
        A_PLSW2A = A_PLSW2A / matrix(sqrt(colSums((A_PLSW2A)^2)), nrow = pX, ncol = r, byrow = TRUE)
        B_PLSW2A = B_PLSW2A / matrix(sqrt(colSums((B_PLSW2A)^2)), nrow = pM, ncol = r, byrow = TRUE)
        
        PA_true_PPLS = best.permutation.matrix(Data$A, theta.hat_PPLS$A_new)
        PB_true_PPLS = best.permutation.matrix(Data$B, theta.hat_PPLS$B_new)
        
        PA_true_PCA = best.permutation.matrix(Data$A, EX$vectors[,1:r])
        PB_true_PCA = best.permutation.matrix(Data$B, EM$vectors[,1:r])
        
        PA_true_PLS = best.permutation.matrix(Data$A, S$u[,1:r])
        PB_true_PLS = best.permutation.matrix(Data$B, S$v[,1:r])
        
        PA_true_PLSW2A = best.permutation.matrix(Data$A, A_PLSW2A)
        PB_true_PLSW2A = best.permutation.matrix(Data$B, B_PLSW2A)
        
        name = paste("colA_true_PPLS", i, sep = "")
        assign(name, abs(diag(t(Data$A) %*% theta.hat_PPLS$A_new %*% t(PA_true_PPLS))[i]))
        COL_A_true_PPLS = rbind(COL_A_true_PPLS, get(name))
        
        name = paste("colB_true_PPLS", i, sep = "")
        assign(name, abs(diag(t(Data$B) %*% theta.hat_PPLS$B_new %*% t(PB_true_PPLS))[i]))
        COL_B_true_PPLS = rbind(COL_B_true_PPLS, get(name))
        
        name = paste("colA_true_PCA", i, sep = "")
        assign(name, abs(diag(t(Data$A) %*% EX$vectors[,1:r] %*% t(PA_true_PCA))[i]))
        COL_A_true_PCA = rbind(COL_A_true_PCA, get(name))
        
        name = paste("colB_true_PCA", i, sep = "")
        assign(name, abs(diag(t(Data$B) %*% EM$vectors[,1:r] %*% t(PB_true_PCA))[i]))
        COL_B_true_PCA = rbind(COL_B_true_PCA, get(name))
        
        name = paste("colA_true_PLS", i, sep = "")
        assign(name, abs(diag(t(Data$A) %*% S$u[,1:r] %*% t(PA_true_PLS))[i]))
        COL_A_true_PLS = rbind(COL_A_true_PLS, get(name))
        
        name = paste("colB_true_PLS", i, sep = "")
        assign(name, abs(diag(t(Data$B) %*% S$v[,1:r] %*% t(PB_true_PLS))[i]))
        COL_B_true_PLS = rbind(COL_B_true_PLS, get(name))
        
        name = paste("colA_true_PLSW2A", i, sep = "")
        assign(name, abs(diag(t(Data$A) %*% A_PLSW2A %*% t(PA_true_PLSW2A))[i]))
        COL_A_true_PLSW2A = rbind(COL_A_true_PLSW2A, get(name))
        
        name = paste("colB_true_PLSW2A", i, sep = "")
        assign(name, abs(diag(t(Data$B) %*% B_PLSW2A %*% t(PB_true_PLSW2A))[i]))
        COL_B_true_PLSW2A = rbind(COL_B_true_PLSW2A, get(name))
        
        name = paste("colA_PPLS_PCA", i, sep = "")
        assign(name, abs(diag(PA_true_PCA %*% t(EX$vectors[,1:r]) %*% theta.hat_PPLS$A_new %*% t(PA_true_PPLS))[i]))
        COL_A_PPLS_PCA = rbind(COL_A_PPLS_PCA, get(name))
        
        name = paste("colB_PPLS_PCA", i, sep = "")
        assign(name, abs(diag(PB_true_PCA %*% t(EM$vectors[,1:r]) %*% theta.hat_PPLS$B_new %*% t(PB_true_PPLS))[i]))
        COL_B_PPLS_PCA = rbind(COL_B_PPLS_PCA, get(name))
        
        name = paste("colA_PPLS_PLS", i, sep = "")
        assign(name, abs(diag(PA_true_PLS %*% t(S$u[,1:r]) %*% theta.hat_PPLS$A_new %*% t(PA_true_PPLS))[i]))
        COL_A_PPLS_PLS = rbind(COL_A_PPLS_PLS, get(name))
        
        name = paste("colB_PPLS_PLS", i, sep = "")
        assign(name, abs(diag(PB_true_PLS %*% t(S$v[,1:r]) %*% theta.hat_PPLS$B_new %*% t(PB_true_PPLS))[i]))
        COL_B_PPLS_PLS = rbind(COL_B_PPLS_PLS, get(name))   
        
        name = paste("colA_PPLS_PLSW2A", i, sep = "")
        assign(name, abs(diag(PA_true_PLSW2A %*% t(A_PLSW2A) %*% theta.hat_PPLS$A_new %*% t(PA_true_PPLS))[i]))
        COL_A_PPLS_PLSW2A = rbind(COL_A_PPLS_PLSW2A, get(name))
        
        name = paste("colB_PPLS_PLSW2A", i, sep = "")
        assign(name, abs(diag(PB_true_PLSW2A %*% t(B_PLSW2A) %*% theta.hat_PPLS$B_new %*% t(PB_true_PPLS))[i]))
        COL_B_PPLS_PLSW2A = rbind(COL_B_PPLS_PLSW2A, get(name))   
        
        col = cbind(COL_A_true_PPLS, COL_B_true_PPLS, COL_A_true_PCA, COL_B_true_PCA, COL_A_true_PLS, COL_B_true_PLS, COL_A_true_PLSW2A, COL_B_true_PLSW2A, COL_A_PPLS_PCA, COL_B_PPLS_PCA, COL_A_PPLS_PLS, COL_B_PPLS_PLS, COL_A_PPLS_PLSW2A, COL_B_PPLS_PLSW2A)
        
      }
      
      sd_A_true_PPLS2 = sd(COL_A_true_PPLS)
      sd_B_true_PPLS2 = sd(COL_B_true_PPLS)
      sd_A_true_PCA2 = sd(COL_A_true_PCA)
      sd_B_true_PCA2 = sd(COL_B_true_PCA)
      sd_A_true_PLS2 = sd(COL_A_true_PLS)
      sd_B_true_PLS2 = sd(COL_B_true_PLS)
      sd_A_true_PLSW2A2 = sd(COL_A_true_PLSW2A)
      sd_B_true_PLSW2A2 = sd(COL_B_true_PLSW2A)
      sd_A_PPLS_PCA2 = sd(COL_A_PPLS_PCA)
      sd_B_PPLS_PCA2 = sd(COL_B_PPLS_PCA)
      sd_A_PPLS_PLS2 = sd(COL_A_PPLS_PLS)
      sd_B_PPLS_PLS2 = sd(COL_B_PPLS_PLS)
      sd_A_PPLS_PLSW2A2 = sd(COL_A_PPLS_PLSW2A)
      sd_B_PPLS_PLSW2A2 = sd(COL_B_PPLS_PLSW2A)
      
      median_A_true_PPLS2 = median(COL_A_true_PPLS)
      median_B_true_PPLS2 = median(COL_B_true_PPLS)
      median_A_true_PCA2 = median(COL_A_true_PCA)
      median_B_true_PCA2 = median(COL_B_true_PCA)
      median_A_true_PLS2 = median(COL_A_true_PLS)
      median_B_true_PLS2 = median(COL_B_true_PLS)
      median_A_true_PLSW2A2 = median(COL_A_true_PLSW2A)
      median_B_true_PLSW2A2 = median(COL_B_true_PLSW2A)
      median_A_PPLS_PCA2 = median(COL_A_PPLS_PCA)
      median_B_PPLS_PCA2 = median(COL_B_PPLS_PCA)
      median_A_PPLS_PLS2 = median(COL_A_PPLS_PLS)
      median_B_PPLS_PLS2 = median(COL_B_PPLS_PLS)
      median_A_PPLS_PLSW2A2 = median(COL_A_PPLS_PLSW2A)
      median_B_PPLS_PLSW2A2 = median(COL_B_PPLS_PLSW2A)
      
      COL_A_PPLS2 = rbind(COL_A_PPLS2, cbind(n = n, column = i, "correctly specified", "true", "COL_A_PPLS", "A", "PPLS", COL_A_true_PPLS, sd_A_true_PPLS2, median_A_true_PPLS2))
      COL_B_PPLS2 = rbind(COL_B_PPLS2, cbind(n = n, column = i, "correctly specified", "true", "COL_B_PPLS", "B", "PPLS", COL_B_true_PPLS, sd_B_true_PPLS2, median_B_true_PPLS2))
      COL_A_PCA2 = rbind(COL_A_PCA2, cbind(n = n, column = i, "correctly specified","true",  "COL_A_PCA", "A", "PCA", COL_A_true_PCA, sd_A_true_PCA2, median_A_true_PCA2))
      COL_B_PCA2 = rbind(COL_B_PCA2, cbind(n = n, column = i, "correctly specified", "true", "COL_B_PCA", "B", "PCA", COL_B_true_PCA, sd_B_true_PCA2, median_B_true_PCA2))
      COL_A_PLS2 = rbind(COL_A_PLS2, cbind(n = n, column = i, "correctly specified", "true", "COL_A_PLS", "A", "PLS-SVD", COL_A_true_PLS, sd_A_true_PLS2, median_B_true_PCA2))
      COL_B_PLS2 = rbind(COL_B_PLS2, cbind(n = n, column = i, "correctly specified", "true", "COL_B_PLS", "B", "PLS-SVD", COL_B_true_PLS, sd_B_true_PLS2, median_B_true_PLS2))
      COL_A_PLSW2A2 = rbind(COL_A_PLSW2A2, cbind(n = n, column = i, "correctly specified", "true", "COL_A_PLSW2A", "A", "PLS-W2A", COL_A_true_PLSW2A, sd_A_true_PLSW2A2, median_A_true_PLSW2A2))
      COL_B_PLSW2A2 = rbind(COL_B_PLSW2A2, cbind(n = n, column = i, "correctly specified", "true", "COL_B_PLSW2A", "B", "PLS-W2A", COL_B_true_PLSW2A, sd_B_true_PLSW2A2, median_B_true_PLSW2A2))
      COL_A_PCA2 = rbind(COL_A_PCA2, cbind(n = n, column = i, "correctly specified", "PPLS", "COL_A_PCA", "A", "PCA", COL_A_PPLS_PCA, sd_A_PPLS_PCA2, median_A_PPLS_PCA2))
      COL_B_PCA2 = rbind(COL_B_PCA2, cbind(n = n, column = i, "correctly specified", "PPLS", "COL_B_PCA", "B", "PCA", COL_B_PPLS_PCA, sd_B_PPLS_PCA2, median_B_PPLS_PCA2))
      COL_A_PLS2 = rbind(COL_A_PLS2, cbind(n = n, column = i, "correctly specified", "PPLS", "COL_A_PLS", "A", "PLS-SVD", COL_A_PPLS_PLS, sd_A_PPLS_PLS2, median_A_PPLS_PLS2))
      COL_B_PLS2 = rbind(COL_B_PLS2, cbind(n = n, column = i, "correctly specified", "PPLS", "COL_B_PLS", "B", "PLS-SVD", COL_B_PPLS_PLS, sd_B_PPLS_PLS2, median_B_PPLS_PLS2))
      COL_A_PLSW2A2 = rbind(COL_A_PLSW2A2, cbind(n = n, column = i, "correctly specified", "PPLS", "COL_A_PLSW2A", "A", "PLS-W2A", COL_A_PPLS_PLSW2A, sd_A_PPLS_PLSW2A2, median_A_PPLS_PLSW2A2))
      COL_B_PLSW2A2 = rbind(COL_B_PLSW2A2, cbind(n = n, column = i, "correctly specified", "PPLS", "COL_B_PLSW2A", "B", "PLS-W2A", COL_B_PPLS_PLSW2A, sd_B_PPLS_PLSW2A2, median_B_PPLS_PLSW2A2))
      
      RES = rbind(RES, cbind(n = n, column = i, ratio = ratio, "correctly specified", col))
      
    }
  }
  
  RES2 = rbind(COL_A_PPLS2, COL_B_PPLS2, COL_A_PCA2, COL_B_PCA2, COL_A_PLS2, COL_B_PLS2, COL_A_PLSW2A2, COL_B_PLSW2A2)
  res = as.data.frame(RES)
  colnames(res) = c("n", "column", "ratio", "data generated from", "COL_A_true_PPLS", "COL_B_true_PPLS", "COL_A_true_PCA", "COL_B_true_PCA", "COL_A_true_PLS", "COL_B_true_PLS", "COL_A_true_PLSW2A", "COL_B_true_PLSW2A", "COL_A_PPLS_PCA", "COL_B_PPLS_PCA", "COL_A_PPLS_PLS", "COL_B_PPLS_PLS", "COL_A_PPLS_PLSW2A", "COL_B_PPLS_PLSW2A")
  res2 = as.data.frame(RES2)
  colnames(res2) = c("n", "Column", "data generated from", "Comparison with", "Type_tot", "Loading", "Method", "Inner Product", "Standard Error", "Median Inner Product")
  res2$Method = as.factor(res2$Method)
  res2$Type_tot = as.factor(res2$Type_tot)
  res2$Loading = as.factor(res2$Loading)
  res2$`Inner Product` = as.numeric(as.character(res2$`Inner Product`))
  res2$`Median Inner Product` = as.numeric(as.character(res2$`Median Inner Product`))
  res2$`Standard Error` = as.numeric(as.character(res2$`Standard Error`))
  res2$n = as.factor(res2$n)
  res2$Column = as.factor(res2$Column)
  res2$`data generated from` = as.factor(res2$`data generated from`)
  res2$`Comparison with` = as.factor(res2$`Comparison with`)
  res2$`Comparison with` = factor(res2$`Comparison with`, levels = c("true", "PPLS"))
  
  myfile  = paste0(dir_data, "Dataframe_weights_comparison-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, ".RData")
  save(res2, file = myfile) # save the dataframe with the results of all the comparisons, for the two simulations.
}



#####
# run_Simul - function to 
## Generate Data (X,M) under the PPLS model proposed by el Bouhaddani et al. (2018), which is:
#          T ~ N(0,Sigma_T)
#          U = T D + e_U,      e_U ~ N(0, sigma_U2 I_r)
#          X = T A' + e_X,     e_X ~ N(0, sigma_X2 I_pX)
#          M = U B' + e_M,     e_M ~ N(0, sigma_M2 I_pM)
# with 
# A and B semi-orthogonal matrices, sigma_U2 = 5.33, sigma_X2 = 0.4, sigma_M2 = 4
# Sigma_T = diag(exp(−(i−1)/5))_{i = 1, ... r} and D = diag(1.5exp(3(i−1)/10))_{i = 1, ... r}.
## Estimate parameters values by the EM algorithm devised under the PPLS model proposed by el Bouhaddani et al. (2018).
## Generate Data (X,M) under the following misspecified PPLS model:
#          T ~ N(0,Sigma_T)
#          U = T D + e_U,      e_U ~ N(0, sigma_U2 I_r)
#          X = T A' + e_X,     e_X ~ N(0, Psi_X)
#          M = U B' + e_M,     e_M ~ N(0, Psi_M)
# with 
## Psi_X and Psi_M arbitrary positive semi-definite matrices.
## Estimate parameters values by the EM algorithm devised under the PPLS model proposed by el Bouhaddani et al. (2018).
## Save generated data sets, true parameters values and estimated values.
## In both cases, compare the true weight matrices A and B with 
#  the estimates returned by the EM alogorithm devised under the PPLS model proposed by el Bouhaddani et al. (2018).
#  the estimates obtained via PLS-SVD on (X,M),
#  the estimates obtained via PLS-W2A on (X,M),
#  the estimates obtained via two distinct PCAs on X and M,
# and compare the estimates returned by the EM alogorithm (devised under the PPLS model proposed by el Bouhaddani et al. (2018))
# with
#  the estimates obtained via PLS-SVD on (X,M),
#  the estimates obtained via PLS-W2A on (X,M),
#  the estimates obtained via two distinct PCAs on X and M,
#
# Inputs:
# dir_data_sim1 - directory where true parameters values, generated data and parameters estimates have been saved for the first simulation study.
# dir_data_sim2 - directory where true parameters values, generated data and parameters estimates have been saved for the second simulation study.
# Nsim - number of replicates chosen for the two simulation studies.
# pX - number of variables in X chosen for the two simulation studies.
# pM - number of variables in M chosen for the two simulation studies.
# r - number of latents variables T and U chosen for the two simulation studies.
# N - sample size(s) chosen for the two simulation studies. Can be a single value or a vector of values. 
# RATIO - Signal to Noise Ratio chosen for the two simulation studies.
# Ncores - number of cores to use to run the data generation and parameters estimation simultaneously.
#
# Outputs:
# Parameters-pX-pM-r-ratio.RData - file containing the true parameters values for the first simulation study.
# Parameters_misspecified-pX-pM-r-ratio.RData - file containing the true parameters values for the second simulation study.
# sample_n-pX-pM-r-ratio-nsim.RData - files containing the generated data for each replicate, for the first simulation study.
# sample_misspecified_n-pX-pM-r-ratio-nsim.RData - files containing the generated data for each replicate, for the second simulation study.
# PPLS_estimation_n-pX-pM-r-ratio-nsim.RData - parameters estimated for each replicate, for the first simulation study.
# PPLS_misspecified_estimation_n-pX-pM-r-ratio-nsim.RData - parameters estimated for each replicate, for the second simulation study.
# Dataframe_weights_comparison-pX-pM-r-ratio.RData - file containing the comparisons decribed above.
###

run_Simul = function(dir_data_sim1, dir_data_sim2, Nsim, pX, pM, r, RATIO, N, Ncores){
  
  if(r >= pX){
    print("r must be strictly inferior to pX")
    break
  }
  if(r >= pM){
    print("r must be strictly inferior to pM")
    break
  }
  
  dir_data = dir_data_sim1
  setwd(dir_data)
  I_r             = Identity(r)   # Identity matrix of size r
  I_pX            = Identity(pX)
  I_pM            = Identity(pM)
  
  # Parameters values which will be used to generate the data
  Sigma_T       = matrix(0, nrow = r, ncol = r)
  for(i in 1:r){
    Sigma_T[i,i]= exp(-(i-1)/5)
  }
  D             = matrix(0, nrow = r, ncol = r)
  for(i in 1:r){
    D[i,i]      = exp(log(1.5)-3*(i-1)/10)
  }
  A               = randortho(n = pX, "orthonormal")
  A               = as.matrix(A[,1:r])
  B               = randortho(n = pM, "orthonormal")
  B               = as.matrix(B[,1:r])
  
  # First simulation study
  for(ratio in RATIO){
    
    sigma_U2      = round(Trace(D^2 %*% Sigma_T)) / (r * ratio)
    sigma_X2      = round(Trace(A %*% Sigma_T %*% t(A))) / (pX * ratio)
    sigma_M2      = round(Trace(B %*% (D^2 %*% Sigma_T + sigma_U2 * I_r) %*% t(B))) / (pM * ratio)
    
    theta.true = list(A = A, B = B, D = D, Sigma_T = Sigma_T, sigma_X2 = sigma_X2, sigma_M2 = sigma_M2, sigma_U2 = sigma_U2)
    
    myfile  = paste0(dir_data, "Parameters-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, ".RData")
    save(theta.true, file = myfile)     # save the true parameters values.
    
    for(n in N){
      Onerun = function(nsim = nsim){
        
        generate_data(dir_data = dir_data, nsim = nsim, n = n, pX = pX, pM = pM, r = r, ratio = ratio, theta.true = theta.true)
        load(paste0(dir_data,"sample_n_" , n, "-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, "-nsim_", nsim, ".RData")) # save generated data
        
        X = Data$X
        M = Data$M
        
        # initialization of the parameters values for the EM algorithm
        theta.init = parameters_initialization(dir_data = NULL, nsim = nsim, n, pM, pX, r, ratio = ratio)
        sigma_X2_0     = theta.init$sigma_X2_init
        sigma_M2_0     = theta.init$sigma_M2_init
        sigma_U2_0     = theta.init$sigma_U2_init
        Sigma_T_0     = as.matrix(theta.init$Sigma_T_init)
        D_0            = as.matrix(theta.init$D_init)
        A_0            = as.matrix(theta.init$A_init)
        B_0            = as.matrix(theta.init$B_init)
        
        delta.stop = NULL 
        Lklh.stop = NULL
        niter.max = NULL
        theta.hat = EM_algo(sigma_X2_0, sigma_M2_0, sigma_U2_0, Sigma_T_0, A_0, B_0, D_0, X, M, delta.stop = NULL, Lklh.stop = NULL, niter.max = NULL)
        theta.hat_PPLS = theta.hat$theta                  # estimates returned by the EM alogorithm
        myfile  = paste0(dir_data,"PPLS_estimation_n_" , n, "-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, "-nsim_", nsim, ".RData")
        save(theta.hat_PPLS, file = myfile)               # save the parameters estimates
      }
      result = mclapply(1:Nsim, Onerun, mc.cores = Ncores)  # parallelization
    }
  }
  
  # Second simulation study
  for(ratio in RATIO){  
    
    dir_data = dir_data_sim1
    setwd(dir_data)
    load(paste0(dir_data, "Parameters-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, ".RData"))
    
    dir_data        = dir_data_sim2
    setwd(dir_data)
    
    A = theta.true$A      # Use the same parameters values A, B, D, Sigma_T and sigma_U2 as in the first simulation study.
    B = theta.true$B
    D = theta.true$D
    Sigma_T = theta.true$Sigma_T
    
    sigma_U2      = round(Trace(D^2 %*% Sigma_T)) / (r * ratio)
    
    # Acceptance-rejection procedure to choose the values of Psi_X and Psi_M below
    nsim = 1
    n = 5000
    crit = 0
    critA = 0
    critB = 0
    while(crit != 1){
      if(critA == 0){
        Mat = matrix(rnorm(pX*pX), nrow = pX, ncol = pX)
        Psi_X = Mat %*% t(Mat)
        Psi_X = Psi_X * Trace(Sigma_T) / Trace(Psi_X) / ratio
      }
      
      if(critB == 0){
        Mat = matrix(rnorm(pM*pM), nrow = pM, ncol = pM)
        Psi_M = Mat %*% t(Mat)
        Psi_M = Psi_M * Trace((D^2 %*% Sigma_T + sigma_U2 * I_r)) / Trace(Psi_M) / ratio
      }
      
      theta.true = list(A = A, B = B, D = D, Sigma_T = Sigma_T, Sigma_X = Psi_X, Sigma_M = Psi_M, sigma_U2 = sigma_U2)
      
      generate_data_misspecified(dir_data = dir_data, nsim = nsim, n = n, pX = pX, pM = pM, r = r, ratio = ratio, theta.true = theta.true)
      load(paste0(dir_data,"sample_misspecified_n_" , n, "-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, "-nsim_", nsim, ".RData"))
      
      X = Data$X
      M = Data$M
      
      EX = eigen(var(X))
      EM = eigen(var(M))
      
      APCA = EX$vectors[,1:r]
      BPCA = EM$vectors[,1:r]
      
      PA_true_PCA = best.permutation.matrix(A, APCA)
      PB_true_PCA = best.permutation.matrix(B, BPCA)
      
      if((sum(abs(diag(t(A) %*% APCA %*% t(PA_true_PCA))) < 0.5) == 3)){
        critA = 0.5
      }
      
      if((sum(abs(diag(t(B) %*% BPCA %*% t(PB_true_PCA))) < 0.5)) == 3){
        critB = 0.5
      }
      
      crit = critA + critB
    }
    
    myfile  = paste0(dir_data, "Parameters_misspecified-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, ".RData")
    save(theta.true, file = myfile)       # save the true parameters values.
    
    for(n in N){
      Onerun = function(nsim = nsim){
        
        # data generation under the misspecified model (presented above)
        generate_data_misspecified(dir_data = dir_data, nsim = nsim, n = n, pX = pX, pM = pM, r = r, ratio = ratio, theta.true = theta.true)
        load(paste0(dir_data,"sample_misspecified_n_" , n, "-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, "-nsim_", nsim, ".RData")) # save generated data
        
        X = Data$X
        M = Data$M
        
        # initialization of the parameters values for the EM algorithm (according to the PPLS model proposed by el Bouhaddani et al. (2018))
        theta.init = parameters_initialization(dir_data = NULL, nsim = nsim, n, pM, pX, r, ratio = ratio)
        sigma_X2_0     = theta.init$sigma_X2_init
        sigma_M2_0     = theta.init$sigma_M2_init
        sigma_U2_0     = theta.init$sigma_U2_init
        Sigma_T_0     = as.matrix(theta.init$Sigma_T_init)
        D_0            = as.matrix(theta.init$D_init)
        A_0            = as.matrix(theta.init$A_init)
        B_0            = as.matrix(theta.init$B_init)
        
        delta.stop = NULL 
        Lklh.stop = NULL
        niter.max = NULL
        theta.hat = EM_algo(sigma_X2_0, sigma_M2_0, sigma_U2_0, Sigma_T_0, A_0, B_0, D_0, X, M, delta.stop = NULL, Lklh.stop = NULL, niter.max = NULL)
        theta.hat_PPLS = theta.hat$theta                  # estimates returned by the EM alogorithm devised under the PPLS model proposed by el Bouhaddani et al. (2018)
        myfile  = paste0(dir_data,"PPLS_misspecified_estimation_n_" , n, "-pX_", pX, "-pM_", pM, "-r_", r, "-ratio_", ratio, "-nsim_", nsim, ".RData")
        save(theta.hat_PPLS, file = myfile)               # save the parameters estimates
      }
      result = mclapply(1:Nsim, Onerun, mc.cores = Ncores)  # parallelization
    }
  }
  
  Comparision_weights(dir_data_sim1, dir_data_sim2, Nsim, pX, pM, r, N, RATIO)
}




