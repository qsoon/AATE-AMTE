library(MASS)
library(R.matlab)
library(rmatio)
# library(GeodRegr)
library(optmatch)
library(senstrat)
library(abind)
library(ggplot2)
library(matrixStats)

source("efficiency.R")
source("GeodRegr.R")
source("normal.R")
source("others.R")


########## private functions: don't need to use these

spdip <- function(pinv, v1, v2) {
  # pinv: inverse of a point p in SPD space (m, m)
  # v1, v2: vectors in T_pM (m, m)
  # result: inner product <v_1, v_2>_p in T_pM
  result <- sum(diag(pinv %*% v1 %*% pinv %*% v2))
  return(result)
}

spdexp2 <- function(phalf, phalfinv, v) {
  # phalf: square root of a point p in SPD space (m, m)
  # phalfinv: inverse of phalf (m, m)
  # v: vector in T_pM (m, m)
  # result : exp_p(v)
  prod <- phalfinv %*% v %*% phalfinv
  prodeigen <- eigen(prod)
  prodexp <- prodeigen$vectors %*% diag(exp(prodeigen$values)) %*% t(prodeigen$vectors)
  result <- phalf %*% prodexp %*% phalf
  return(result)
}

spdlog2 <- function(phalf, phalfinv, x) {
  # phalf: square root of a point p in SPD space (m, m)
  # phalfinv: inverse of phalf (m, m)
  # x: point in SPD space (m, m)
  # result: log_p(x) in T_pM (m, m)
  prod <- phalfinv %*% x %*% phalfinv
  prodeigen <- eigen(prod)
  prodlog <- prodeigen$vectors %*% diag(log(prodeigen$values)) %*% t(prodeigen$vectors)
  result <- phalf %*% prodlog %*% phalf
  return(result)
}

spddistance2 <- function(phalf, pinv, phalfinv, x) {
  # phalf: square root of a point p in SPD space (m, m)
  # pinv: inverse of a point p in SPD space (m, m)
  # phalfinv: inverse of phalf (m, m)
  # x: point in SPD space (m, m)
  # result: Riemannian distance between p and x
  v <- spdlog2(phalf, phalfinv, x)
  result <- spdip(pinv, v, v)^0.5
  return(result)
}

spdL2loss <- function(phalf, pinv, phalfinv, X, W) {
  # phalf: square root of a point p in SPD space (m, m)
  # pinv: inverse of a point p in SPD space (m, m)
  # phalfinv: inverse of phalf (m, m)
  # X: collection of N mxm SPD matrices (m, m, N)
  # W: weights for points in X; should sum to 1 (N)
  # result: squared loss function for p, X
  result <- 0
  for (i in 1:dim(X)[3]) {
    result <- result+W[i]*spddistance2(phalf, pinv, phalfinv, X[, , i])^2
  }
  return(result)
}

spdL1loss <- function(phalf, pinv, phalfinv, X, W) {
  # phalf: square root of a point p in SPD space (m, m)
  # pinv: inverse of a point p in SPD space (m, m)
  # phalfinv: inverse of phalf (m, m)
  # X: collection of N mxm SPD matrices (m, m, N)
  # W: weights for points in X; should sum to 1 (N)
  # result: absolute loss function for p, X
  result <- 0
  for (i in 1:dim(X)[3]) {
    result <- result+W[i]*spddistance2(phalf, pinv, phalfinv, X[, , i])
  }
  return(result)
}

spdL2grad <- function(phalf, pinv, phalfinv, X, W) {
  # phalf: square root of a point p in SPD space (m, m)
  # pinv: inverse of a point p in SPD space (m, m)
  # phalfinv: inverse of phalf (m, m)
  # X: collection of N mxm SPD matrices (m, m, N)
  # W: weights for points in X; should sum to 1 (N)
  # result: gradient in T_pM (m, m)
  m <- dim(X)[1]
  result <- matrix(numeric(m*m), nrow = m)
  for (i in 1:dim(X)[3]) {
    result <- result-W[i]*spdlog2(phalf, phalfinv, X[, , i])
  }
  return(result)
}

spdL1grad <- function(phalf, pinv, phalfinv, X, W) {
  # phalf: square root of a point p in SPD space (m, m)
  # pinv: inverse of a point p in SPD space (m, m)
  # phalfinv: inverse of phalf (m, m)
  # X: collection of N mxm SPD matrices (m, m, N)
  # W: weights for points in X; should sum to 1 (N)
  # result: gradient in T_pM (m, m)
  m <- dim(X)[1]
  result <- matrix(numeric(m*m), nrow = m)
  for (i in 1:dim(X)[3]) {
    vect <- spdlog2(phalf, phalfinv, X[, , i])
    norm <- spdip(pinv, vect, vect)^0.5
    if (norm > 1e-10) {
      result <- result-W[i]*vect/norm
    }
  }
  return(result)
}

########## public functions: use these
  
spdL2center <- function(X, W, tol=1e-10) {
  # X: collection of N mxm SPD matrices (m, m, N)
  # W: weights for points in X; should sum to 1 (N)
  # result: Frechet mean (m, m)
  current_p <- apply(X, c(1,2), mean)
  current_peigen <- eigen(current_p)
  current_phalf <-current_peigen$vectors %*% diag(current_peigen$values^0.5) %*% t(current_peigen$vectors)
  current_pinv <- solve(current_p)
  current_phalfinv <- solve(current_phalf)
  current_loss <- spdL2loss(current_phalf, current_pinv, current_phalfinv, X, W)
  step <- spdL2grad(current_phalf, current_pinv, current_phalfinv, X, W)
  lr <- 1
  count <- 0
  acount <- 0
  while ((spdip(current_pinv, step, step)>tol) & count<100) {
    new_p <- spdexp2(current_phalf, current_phalfinv, -lr*step)
    new_peigen <- eigen(new_p)
    new_phalf <-new_peigen$vectors %*% diag(new_peigen$values^0.5) %*% t(new_peigen$vectors)
    new_pinv <- solve(new_p)
    new_phalfinv <- solve(new_phalf)
    new_loss <- spdL2loss(new_phalf, new_pinv, new_phalfinv, X, W)
    if (new_loss < current_loss) {
      current_p <- new_p
      current_phalf <- new_phalf
      current_pinv <- new_pinv
      current_phalfinv <- new_phalfinv
      current_loss <- new_loss
      step <- spdL2grad(current_phalf, current_pinv, current_phalfinv, X, W)
      lr <- 1.1*lr
      count <- count+1
    } else {
      lr <- lr/2
      acount <- acount+1
    }
  }
  result <- current_p
  #print(count)
  #print(acount)
  return(result)
}

spdL1center <- function(X, W, tol=1e-10) {
  # X: collection of N mxm SPD matrices (m, m, N)
  # W: weights for points in X; should sum to 1 (N)
  # result: Frechet median (m, m)
  current_p <- apply(X, c(1,2), mean)
  current_peigen <- eigen(current_p)
  current_phalf <-current_peigen$vectors %*% diag(current_peigen$values^0.5) %*% t(current_peigen$vectors)
  current_pinv <- solve(current_p)
  current_phalfinv <- solve(current_phalf)
  current_loss <- spdL1loss(current_phalf, current_pinv, current_phalfinv, X, W)
  step <- spdL1grad(current_phalf, current_pinv, current_phalfinv, X, W)
  lr <- 1
  count <- 0
  acount <- 0
  while ((spdip(current_pinv, step, step)>tol) & count<100) {
    new_p <- spdexp2(current_phalf, current_phalfinv, -lr*step)
    new_peigen <- eigen(new_p)
    new_phalf <-new_peigen$vectors %*% diag(new_peigen$values^0.5) %*% t(new_peigen$vectors)
    new_pinv <- solve(new_p)
    new_phalfinv <- solve(new_phalf)
    new_loss <- spdL1loss(new_phalf, new_pinv, new_phalfinv, X, W)
    if (new_loss < current_loss) {
      current_p <- new_p
      current_phalf <- new_phalf
      current_pinv <- new_pinv
      current_phalfinv <- new_phalfinv
      current_loss <- new_loss
      step <- spdL1grad(current_phalf, current_pinv, current_phalfinv, X, W)
      lr <- 1.1*lr
      count <- count+1
    } else {
      lr <- lr/2
      acount <- acount+1
    }
  }
  result <- current_p
  #print(count)
  #print(acount)
  return(result)
}

spdexp <- function(p, v) {
  # p: a point p in SPD space (m, m)
  # v: vector in T_pM (m, m)
  # result : exp_p(v)
  peigen <- eigen(p)
  phalf <- peigen$vectors %*% diag(peigen$values^0.5) %*% t(peigen$vectors)
  phalfinv <- solve(phalf)
  prod <- phalfinv %*% v %*% phalfinv
  prodeigen <- eigen(prod)
  prodexp <- prodeigen$vectors %*% diag(exp(prodeigen$values)) %*% t(prodeigen$vectors)
  result <- phalf %*% prodexp %*% phalf
  return(result)
}

spdlog <- function(p, x) {
  # p: a point p in SPD space (m, m)
  # x: point in SPD space (m, m)
  # result: log_p(x) in T_pM (m, m)
  peigen <- eigen(p)
  phalf <- peigen$vectors %*% diag(peigen$values^0.5) %*% t(peigen$vectors)
  phalfinv <- solve(phalf)
  prod <- phalfinv %*% x %*% phalfinv
  prodeigen <- eigen(prod)
  prodlog <- prodeigen$vectors %*% diag(log(prodeigen$values)) %*% t(prodeigen$vectors)
  result <- phalf %*% prodlog %*% phalf
  return(result)
}

spddistance <- function(p, x) {
  # p: a point p in SPD space (m, m)
  # x: point in SPD space (m, m)
  # result: Riemannian distance between p and x
  peigen <- eigen(p)
  phalf <- peigen$vectors %*% diag(peigen$values^0.5) %*% t(peigen$vectors)
  pinv <- solve(p)
  phalfinv <- solve(phalf)
  v <- spdlog2(phalf, phalfinv, x)
  result <- spdip(pinv, v, v)^0.5
  return(result)
}

# Note that the following is an orthonormal basis in the tangent space at the identity
# [1,0,0]  (1/sqrt(2)) * [0,1,0]  (1/sqrt(2)) * [0,0,1]  [0,0,0]  (1/sqrt(2)) * [0,0,0]  [0,0,0]
# [0,0,0]                [1,0,0]                [0,0,0]  [0,1,0]                [0,0,1]  [0,0,0]
# [0,0,0],               [0,0,0],               [1,0,0], [0,0,0],               [0,1,0], [0,0,1]

vectorform <- function(A) {
  # A: symmetric matrix (3, 3)
  # result : A written as a vector in terms of the above orthonormal basis (6)
  a <- A[1,1]
  b <- A[1,2]*sqrt(2)
  c <- A[1,3]*sqrt(2)
  d <- A[2,2]
  e <- A[2,3]*sqrt(2)
  f <- A[3,3]
  return(c(a,b,c,d,e,f))
}

fa <- function(p) {
  # p: a point p in 3x3 SPD space (3, 3)
  # result: fractional anisotropy of p
  peigenvalues <- eigen(p)$values
  result <- var(peigenvalues)*3/sum(peigenvalues^2)
  return(result)
}

md <- function(p) {
  # p: a point p in 3x3 SPD space (3, 3)
  # result: mean diffusitivity of p
  result <- mean(eigen(p)$values)
  return(result)
}


calculateWeight <- function(z, strata, beta.T, beta.C){
  N <- length(z)
  L <- length(unique(strata))
  m.T <- as.numeric(table(strata[z==1]))
  m.C <- as.numeric(table(strata[z==0]))
  lambda.hat <- as.numeric(table(strata)) / N
  # beta.T <- sum(z) / N
  # beta.C <- 1-beta.T
  wt <- beta.T*((lambda.hat/m.T)[strata])*z + beta.C*((lambda.hat/m.C)[strata])*(1-z)
  return(wt)
}

calculateWeight.T <- function(z, strata){
  N <- length(z)
  L <- length(unique(strata))
  m.T <- as.numeric(table(strata[z==1]))
  lambda.hat <- as.numeric(table(strata)) / N
  wt <- ((lambda.hat/m.T)[strata])*z 
  return(wt)
}

calculateWeight.C <- function(z, strata){
  N <-length(z)
  L <- length(unique(strata))
  m.C <- as.numeric(table(strata[z==0]))
  lambda.hat <- as.numeric(table(strata)) / N
  wt <- ((lambda.hat/m.C)[strata])*(1-z)
  return(wt)
}

calculateWeight.T.sre <- function(z, strata){
  N <- length(z)
  L <- length(unique(strata))
  m.T <- as.numeric(table(strata[z==1]))
  lambda.hat <- 0.5
  wt <- ((lambda.hat/m.T)[strata])*z 
  return(wt)
}

calculateWeight.C.sre <- function(z, strata){
  N <-length(z)
  L <- length(unique(strata))
  m.C <- as.numeric(table(strata[z==0]))
  lambda.hat <- 0.5
  wt <- ((lambda.hat/m.C)[strata])*(1-z)
  return(wt)
}

getpval_spd <- function(z, strata, y.list, T.obs, B=500, estimator, iter_total){
  N <- length(z)
  L <- length(unique(strata))
  strata.list <- list()
  z.list <- list()
  for(i in 1:L){
    strata.list[[i]] <- which(strata==i)
    z.list[[i]] <- z[which(strata==i)]
  }
  
  T.shuffle <- matrix(0, nrow=B, ncol=5)
  p.val <- numeric(5)
  if(estimator=="l1"){
    for(b in 1:B){
      set.seed(iter_total*10000 + b*10)
      print(paste("total_iter:", iter_total, "iter:", b))
      
      z.shuffle <- vector(length=N)
      z.shuffle.list <- lapply(z.list, function(x) x[sample(1:length(x))])
      z.shuffle[unlist(strata.list)] <- unlist(z.shuffle.list)
      
      wt.t.shuffle <- calculateWeight.T.sre(z.shuffle, strata)
      wt.c.shuffle <- calculateWeight.C.sre(z.shuffle, strata)
      if((sum(is.na(wt.t.shuffle))!=0) | (!isTRUE(all.equal(sum(wt.t.shuffle), 1)))){
        print("check!")
        next
      }
      if((sum(is.na(wt.c.shuffle))!=0) | (!isTRUE(all.equal(sum(wt.c.shuffle), 1)))){
        print("check!")
        next
      }
      
      # method 1 (proposed)
      T.shuffle[b,1] <- spddistance(spdL1center(y.list[[1]], wt.t.shuffle), spdL1center(y.list[[1]], wt.c.shuffle))
      if(T.shuffle[b,1] >= T.obs[1]){
        p.val[1] <- p.val[1] + 1
      }
      # method 2 (through vectorform)
      T.shuffle[b,2] <- geo_dist(manifold="euclidean",intrinsic_location(manifold="euclidean", y.list[[2]], w = wt.t.shuffle, estimator = "l1", p_tol=1e-6, V_tol=1e-6),
                                 intrinsic_location(manifold="euclidean", y.list[[2]], w = wt.c.shuffle, estimator = "l1", p_tol=1e-6, V_tol=1e-6))
      if(T.shuffle[b,2] >= T.obs[2]){
        p.val[2] <- p.val[2] + 1
      }
      # method 3 (through spdlog + vectorform)
      T.shuffle[b,3] <- geo_dist(manifold="euclidean",intrinsic_location(manifold="euclidean", y.list[[3]], w = wt.t.shuffle, estimator = "l1", p_tol=1e-6, V_tol=1e-6),
                                 intrinsic_location(manifold="euclidean", y.list[[3]], w = wt.c.shuffle, estimator = "l1", p_tol=1e-6, V_tol=1e-6))
      if(T.shuffle[b,3] >= T.obs[3]){
        p.val[3] <- p.val[3] + 1
      }
      # method 4 (through fa)
      T.shuffle[b,4] <- abs(weightedMedian(y.list[[4]], wt.t.shuffle) - weightedMedian(y.list[[4]], wt.c.shuffle))
      if(T.shuffle[b,4] >= T.obs[4]){
        p.val[4] <- p.val[4] + 1
      }
      # method 5 (through md)
      T.shuffle[b,5] <- abs(weightedMedian(y.list[[5]], wt.t.shuffle) - weightedMedian(y.list[[5]], wt.c.shuffle))
      if(T.shuffle[b,5] >= T.obs[5]){
        p.val[5] <- p.val[5] + 1
      }
    }
  } else if(estimator=="l2"){
    for(b in 1:B){
      set.seed(iter_total*10000 + b*10)
      print(paste("total_iter:", iter_total, "iter:", b))
      z.shuffle <- vector(length=N)
      z.shuffle.list <- lapply(z.list, function(x) x[sample(1:length(x))])
      z.shuffle[unlist(strata.list)] <- unlist(z.shuffle.list)
      
      wt.t.shuffle <- calculateWeight.T.sre(z.shuffle, strata)
      wt.c.shuffle <- calculateWeight.C.sre(z.shuffle, strata)
      if((sum(is.na(wt.t.shuffle))!=0) | (!isTRUE(all.equal(sum(wt.t.shuffle), 1)))){
        print("check!")
        next
      }
      if((sum(is.na(wt.c.shuffle))!=0) | (!isTRUE(all.equal(sum(wt.c.shuffle), 1)))){
        print("check!")
        next
      }
      
      # method 1 (proposed)
      T.shuffle[b,1] <- spddistance(spdL2center(y.list[[1]], wt.t.shuffle), spdL2center(y.list[[1]], wt.c.shuffle))
      if(T.shuffle[b,1] >= T.obs[1]){
        p.val[1] <- p.val[1] + 1
      }
      
      # method 2 (through vectorform)
      T.shuffle[b,2] <- geo_dist(manifold="euclidean",intrinsic_location(manifold="euclidean", y.list[[2]], w = wt.t.shuffle, estimator = "l1", p_tol=1e-6, V_tol=1e-6),
                                 intrinsic_location(manifold="euclidean", y.list[[2]], w = wt.c.shuffle, estimator = "l2", p_tol=1e-6, V_tol=1e-6))
      if(T.shuffle[b,2] >= T.obs[2]){
        p.val[2] <- p.val[2] + 1
      }
      # method 3 (through spdlog + vectorform)
      T.shuffle[b,3] <- geo_dist(manifold="euclidean",intrinsic_location(manifold="euclidean", y.list[[3]], w = wt.t.shuffle, estimator = "l1", p_tol=1e-6, V_tol=1e-6),
                                 intrinsic_location(manifold="euclidean", y.list[[3]], w = wt.c.shuffle, estimator = "l2", p_tol=1e-6, V_tol=1e-6))
      if(T.shuffle[b,3] >= T.obs[3]){
        p.val[3] <- p.val[3] + 1
      }
      # method 4 (through fa)
      T.shuffle[b,4] <- abs(weightedMean(y.list[[4]], wt.t.shuffle) - weightedMean(y.list[[4]], wt.c.shuffle))
      if(T.shuffle[b,4] >= T.obs[4]){
        p.val[4] <- p.val[4] + 1
      }
      # method 5 (through md)
      T.shuffle[b,5] <- abs(weightedMean(y.list[[5]], wt.t.shuffle) - weightedMean(y.list[[5]], wt.c.shuffle))
      if(T.shuffle[b,5] >= T.obs[5]){
        p.val[5] <- p.val[5] + 1
      }
    }
  }
  
  p.val <- p.val / B
  
  res <- list()
  res$rejection_region_proposed <- quantile(T.shuffle[,1], probs = 0.95)
  res$pval <- p.val
  return(res)  
}

# data generation


##########################
## SRE
##########################
m <- 3 # 3 x 3 matrices
N <- 80 # number of patients
v1 <- matrix(c(0,1,0,1,0,0,0,0,0), nrow=3, ncol=3)
v2 <- matrix(c(0,0,1,0,0,0,1,0,0), nrow=3, ncol=3)
mu_T <- spdexp(diag(3),0.5*matrix(c(1,0,0,0,0,0,0,0,0), nrow=3, ncol=3)) # treatment mean
mu_C <- spdexp(diag(3),0.5*matrix(c(-1,0,0,0,0,0,0,0,0), nrow=3, ncol=3)) # control mean
spddistance(mu_T, mu_C)
J <- 500 # total number of simulations

res.mat.l1.spd <- matrix(0, nrow=500, ncol=11)
res.mat.l2.spd <- matrix(0, nrow=500, ncol=11)

for (j in 1:J) {
  set.seed(10000*j)
  covariates <- matrix(runif(2*N),nrow=N, ncol=2) - 0.5
  r_Ts <- array(numeric(m*m*N), dim=c(m, m, N)) # (m, m, N) array where ith layer is the control result for patient i
  r_Cs <- array(numeric(m*m*N), dim=c(m, m, N)) # (m, m, N) array where ith layer is the treatment result for patient i
  for (i in 1:N) {
    noise <- matrix(rnorm(m^2, 0, 1), nrow=m, ncol=m)
    noise <- (noise+t(noise))/2 # to ensure noise is symmetric
    r_Ts[,,i] <- spdexp(mu_T, covariates[i,1]*v1+covariates[i,2]*v2+noise)
    r_Cs[,,i] <- spdexp(mu_C, covariates[i,1]*v1+covariates[i,2]*v2+noise)
  }
  
  strata.sre <- ifelse(covariates[,1] >= 0, 1, 2) # make 2 strata based on covariate x1
  strata.sre1.ind <- which(strata.sre==1)
  len.s1 <- length(strata.sre1.ind)
  strata.sre2.ind <- which(strata.sre==2)
  len.s2 <- length(strata.sre2.ind)
  
  z.sre <- rep(0,N) # independent of (r_t, r_c) given S, here, 50:50 in each stratum
  z.sre[sort(sample(strata.sre1.ind, floor((len.s1+1)/2), replace=FALSE))] <- 1 
  z.sre[sort(sample(strata.sre2.ind, floor((len.s2+1)/2), replace=FALSE))] <- 1
  
  wt.t.sre <- calculateWeight.T.sre(z.sre, strata.sre)
  wt.c.sre <- calculateWeight.C.sre(z.sre, strata.sre)
  if((sum(is.na(wt.t.sre))!=0) | (!isTRUE(all.equal(sum(wt.t.sre), 1)))){
    print("check!")
    next
  }
  if((sum(is.na(wt.c.sre))!=0) | (!isTRUE(all.equal(sum(wt.c.sre), 1)))){
    print("check!")
    next
  }
  
  y.sre <- array(numeric(m*m*N), dim=c(m, m, N))
  y.sre[,,which(z.sre==1)] <- r_Ts[,,which(z.sre==1)]
  y.sre[,,which(z.sre==0)] <- r_Cs[,,which(z.sre==0)]
  
  # method 1 (proposed)
  T.l1.proposed <- spddistance(spdL1center(y.sre, wt.t.sre), spdL1center(y.sre, wt.c.sre))
  T.l2.proposed <- spddistance(spdL2center(y.sre, wt.t.sre), spdL2center(y.sre, wt.c.sre))
  
  y.sre.vectorform <- matrix(0, nrow=6, ncol=N)
  y.sre.spdlog <- matrix(0, nrow=6, ncol=N)
  y.sre.fa <- numeric(N)
  y.sre.md <- numeric(N)
  for(i in 1:N){
    y.sre.vectorform[,i] <- vectorform(y.sre[,,i])
    y.sre.spdlog[,i] <- vectorform(spdlog(diag(3), y.sre[,,i]))
    y.sre.fa[i] <- fa(y.sre[,,i])
    y.sre.md[i] <- md(y.sre[,,i])
  }
  
  y.sre.list <- list()
  y.sre.list[[1]] <- y.sre
  y.sre.list[[2]] <- y.sre.vectorform
  y.sre.list[[3]] <- y.sre.spdlog
  y.sre.list[[4]] <- y.sre.fa
  y.sre.list[[5]] <- y.sre.md
  
  
  # method 2 (through vectorform)
  T.l1.vectorform <- geo_dist(manifold="euclidean",intrinsic_location(manifold="euclidean", y.sre.vectorform, w = wt.t.sre, estimator = "l1", p_tol=1e-6, V_tol=1e-6),
                              intrinsic_location(manifold="euclidean", y.sre.vectorform, w = wt.c.sre, estimator = "l1", p_tol=1e-6, V_tol=1e-6))
  T.l2.vectorform <- geo_dist(manifold="euclidean",intrinsic_location(manifold="euclidean", y.sre.vectorform, w = wt.t.sre, estimator = "l2", p_tol=1e-6, V_tol=1e-6),
                              intrinsic_location(manifold="euclidean", y.sre.vectorform, w = wt.c.sre, estimator = "l2", p_tol=1e-6, V_tol=1e-6))
  
  # method 3 (through spdlog + vectorform)
  T.l1.spdlog <- geo_dist(manifold="euclidean",intrinsic_location(manifold="euclidean", y.sre.spdlog, w = wt.t.sre, estimator = "l1", p_tol=1e-6, V_tol=1e-6),
                              intrinsic_location(manifold="euclidean", y.sre.spdlog, w = wt.c.sre, estimator = "l1", p_tol=1e-6, V_tol=1e-6))
  T.l2.spdlog <- geo_dist(manifold="euclidean",intrinsic_location(manifold="euclidean", y.sre.spdlog, w = wt.t.sre, estimator = "l2", p_tol=1e-6, V_tol=1e-6),
                              intrinsic_location(manifold="euclidean", y.sre.spdlog, w = wt.c.sre, estimator = "l2", p_tol=1e-6, V_tol=1e-6))
  
  # method 4 (through fa)
  T.l1.fa <- abs(weightedMedian(y.sre.fa, wt.t.sre) - weightedMedian(y.sre.fa, wt.c.sre))
  T.l2.fa <- abs(weightedMean(y.sre.fa, wt.t.sre) - weightedMean(y.sre.fa, wt.c.sre))
  
  # method 5 (through md)
  T.l1.md <- abs(weightedMedian(y.sre.md, wt.t.sre) - weightedMedian(y.sre.md, wt.c.sre))
  T.l2.md <- abs(weightedMean(y.sre.md, wt.t.sre) - weightedMean(y.sre.md, wt.c.sre))
  
  T.l1.collection <- c(T.l1.proposed, T.l1.vectorform, T.l1.spdlog,
                       T.l1.fa, T.l1.md)
  T.l2.collection <- c(T.l2.proposed, T.l2.vectorform, T.l2.spdlog,
                       T.l2.fa, T.l2.md)
  
  res.dist.l1.spd <- getpval_spd(z.sre, strata.sre, y.sre.list, T.l1.collection, B=500, estimator = 'l1', iter_total=j)
  res.dist.l2.spd <- getpval_spd(z.sre, strata.sre, y.sre.list, T.l2.collection, B=500, estimator = 'l2', iter_total=j)
  
  res.mat.l1.spd[j,] <- c(res.dist.l1.spd$rejection_region_proposed, res.dist.l1.spd$pval, T.l1.collection) 
  res.mat.l2.spd[j,] <- c(res.dist.l2.spd$rejection_region_proposed, res.dist.l2.spd$pval, T.l2.collection) 
  
  # do both randomized experiments and matched observational studies
  # assign strata for randomized experiments, and then assign treatment in the exact same way as in our existing experiments
  # test Fisher's sharp null (with 500 permutations) in 5 ways, comparing our method with 4 Euclidean alternatives. Use the same strata/weights each time:
  # 1. using our manifold-based test statistic using the spdL2center(), spdL1center() and spddistance() functions.
  # 2. put all observations through vectorform() to get 6D vectors, then do multivariate Euclidean version.
  # 3. put all observations through spdlog(diag(3), ...) and then vectorform() to get 6D vectors in tangent space at diag(3), then do multivariate Euclidean version.
  # 4. put all observations through fa() to get scalars, then do univariate Euclidean version.
  # 5. put all observations through md() to get scalars, then do univariate Euclidean version.
  # save p-values for each of these 5 tests for each j. also save the 95\% rejection region for our manifold-based test for each j.
}

colnames(res.mat.l1.spd) <- c("rejection_region", "proposed.pval","vectorform.pval",
                              "spdlog.pval", "fa.pval", "md.pval", "proposed.T","vectorform.T",
                              "spdlog.T", "fa.T", "md.T")
colnames(res.mat.l2.spd) <- c("rejection_region", "proposed.pval","vectorform.pval",
                              "spdlog.pval", "fa.pval", "md.pval", "proposed.T","vectorform.T",
                              "spdlog.T", "fa.T", "md.T")



##########################
## Observational study
##########################
library(MASS)
rankmahal <- function(z,X){
  X <- as.matrix(X)
  n <- dim(X)[1]
  rownames(X) <- 1:n
  k <- dim(X)[2]
  m <- sum(z)
  
  for(j in 1:k){
    X[,j] <- rank(X[,j])
  }
  cv <- cov(X)
  vuntied <- var(1:n)
  rat <- sqrt(vuntied/diag(cv))
  cv <- diag(rat)%*%cv%*%diag(rat)
  out <- matrix(NA,m,n-m)
  Xc <- X[z==0,]
  Xt <- X[z==1,]
  rownames(out) <- rownames(X)[z==1]
  colnames(out) <- rownames(X)[z==0]
  
  icov <- ginv(cv)
  if(m==1){
    out[1,] <- mahalanobis(Xc, Xt, icov, inverted=T)
  } else{
    for(i in 1:m){
      out[i,] <- mahalanobis(Xc, Xt[i,], icov, inverted=T)
    }
  }
  
  
  return(out)
}



addcaliper <- function(dmat, z, logitp, calipersd=.2, penalty=1000){
  sd.logitp <- sd(logitp)
  adif <- abs(outer(logitp[z==1],logitp[z==0],"-"))
  adif <- (adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
  dmat <- dmat+adif*penalty
  return(dmat)
}

m <- 3 # 3 x 3 matrices
N <- 100 # number of patients
v1 <- matrix(c(0,1,0,1,0,0,0,0,0), nrow=3, ncol=3)
v2 <- matrix(c(0,0,1,0,0,0,1,0,0), nrow=3, ncol=3)
mu_T <- spdexp(diag(3), 0.5*matrix(c(1,0,0,0,0,0,0,0,0), nrow=3, ncol=3)) # treatment mean
mu_C <- spdexp(diag(3), 0.5*matrix(c(-1,0,0,0,0,0,0,0,0), nrow=3, ncol=3)) # control mean
J <- 500 # total number of simulations

res.mat.l1.spd.os <- matrix(0, nrow=J, ncol=11)
res.mat.l2.spd.os <- matrix(0, nrow=J, ncol=11)

for (j in 1:J) {
  set.seed(10000*j)
  covariates <- matrix(runif(2*N),nrow=N, ncol=2) - 0.5
  r_Ts <- array(numeric(m*m*N), dim=c(m, m, N)) # (m, m, N) array where ith layer is the control result for patient i
  r_Cs <- array(numeric(m*m*N), dim=c(m, m, N)) # (m, m, N) array where ith layer is the treatment result for patient i
  for (i in 1:N) {
    noise <- matrix(rnorm(m^2, 0, 1), nrow=m, ncol=m)
    noise <- (noise+t(noise))/2 # to ensure noise is symmetric
    r_Ts[,,i] <- spdexp(mu_T, covariates[i,1]*v1+covariates[i,2]*v2+noise)
    r_Cs[,,i] <- spdexp(mu_C, covariates[i,1]*v1+covariates[i,2]*v2+noise)
  }
  
  z.os <- rbinom(N, 1, 1 / (1+exp(-(covariates[,1]+covariates[,2])))) # dependent on (r_t, r_c)
  
  y.os <- array(numeric(m*m*N), dim=c(m, m, N))
  y.os[,,which(z.os==1)] <- r_Ts[,,which(z.os==1)]
  y.os[,,which(z.os==0)] <- r_Cs[,,which(z.os==0)]
  
  x_data.os <- data.frame(V1 = covariates[,1], V2 = covariates[,2], Z = z.os)
  
  # matching
  propscore.model.os <- glm(Z ~ V1 + V2,
                            family=binomial, x=TRUE, data=x_data.os)
  
  Xmat.all.os <- (propscore.model.os$x[,-1])
  
  # Rank based Mahalanobis distance
  distmat.all.os <- rankmahal(z.os, Xmat.all.os)
  # Add caliper
  logit.propscore.all.os <- (predict(propscore.model.os))
  distmat2.all.os <- addcaliper(distmat.all.os, z.os, logit.propscore.all.os) 
  
  ## full matching 
  matchvec.full.all.os <- fullmatch(distmat2.all.os)
  
  treated.subject.index.full.all.os <- rep(0,sum(z.os==1))
  # The subject indices in the order of matchvec
  matchedset.index.full.all.os <- substr(matchvec.full.all.os, start=3, stop=10) # group name
  matchedset.index.full.numeric.all.os <- as.numeric(matchedset.index.full.all.os) # group name (numeric)
  subjects.match.order.full.all.os <- as.numeric(names(matchvec.full.all.os)) # subject indices
  
  # Create a numeric variable for which stratum each unit belongs to
  # 0 denotes that the unit was not matched
  # there is not units not matched b/c we use full matching
  stratum.short.full.all.os <- substr(matchvec.full.all.os, start=3, stop=10)
  stratum.full.numeric.all.os <- as.numeric(stratum.short.full.all.os) # i-th element : group name of unit i
  
  # Reassign numbers to each stratum that go from 1,..., no. of straum
  sort.unique.stratum.full.all.os <- sort(unique(stratum.full.numeric.all.os)) 
  stratum.myindex.matchvecorder.full.all.os <- rep(0,length(stratum.full.numeric.all.os))
  for(ii in 1:length(sort.unique.stratum.full.all.os)){
    stratum.myindex.matchvecorder.full.all.os[stratum.full.numeric.all.os==sort.unique.stratum.full.all.os[ii]] <- ii  
  }
  
  stratum.myindex.full.all.os <- rep(0, length(stratum.myindex.matchvecorder.full.all.os)) 
  stratum.myindex.full.all.os[subjects.match.order.full.all.os] <- stratum.myindex.matchvecorder.full.all.os # i-th element : matched group index for unit i 
  
  
  group.use.os <- stratum.myindex.full.all.os
  
  wt.t.os <- calculateWeight.T(z.os, group.use.os)
  wt.c.os <- calculateWeight.C(z.os, group.use.os)
  if((sum(is.na(wt.t.os))!=0) | (!isTRUE(all.equal(sum(wt.t.os), 1)))){
    print("check!")
    next
  }
  if((sum(is.na(wt.c.os))!=0) | (!isTRUE(all.equal(sum(wt.c.os), 1)))){
    print("check!")
    next
  }
  
  # method 1 (proposed)
  T.l1.proposed <- spddistance(spdL1center(y.os, wt.t.os), spdL1center(y.os, wt.c.os))
  T.l2.proposed <- spddistance(spdL2center(y.os, wt.t.os), spdL2center(y.os, wt.c.os))
  
  y.os.vectorform <- matrix(0, nrow=6, ncol=N)
  y.os.spdlog <- matrix(0, nrow=6, ncol=N)
  y.os.fa <- numeric(N)
  y.os.md <- numeric(N)
  for(i in 1:N){
    y.os.vectorform[,i] <- vectorform(y.os[,,i])
    y.os.spdlog[,i] <- vectorform(spdlog(diag(3), y.os[,,i]))
    y.os.fa[i] <- fa(y.os[,,i])
    y.os.md[i] <- md(y.os[,,i])
  }
  
  y.os.list <- list()
  y.os.list[[1]] <- y.os
  y.os.list[[2]] <- y.os.vectorform
  y.os.list[[3]] <- y.os.spdlog
  y.os.list[[4]] <- y.os.fa
  y.os.list[[5]] <- y.os.md
  
  
  # method 2 (through vectorform)
  T.l1.vectorform <- geo_dist(manifold="euclidean",intrinsic_location(manifold="euclidean", y.os.vectorform, w = wt.t.os, estimator = "l1", p_tol=1e-6, V_tol=1e-6),
                              intrinsic_location(manifold="euclidean", y.os.vectorform, w = wt.c.os, estimator = "l1", p_tol=1e-6, V_tol=1e-6))
  T.l2.vectorform <- geo_dist(manifold="euclidean",intrinsic_location(manifold="euclidean", y.os.vectorform, w = wt.t.os, estimator = "l2", p_tol=1e-6, V_tol=1e-6),
                              intrinsic_location(manifold="euclidean", y.os.vectorform, w = wt.c.os, estimator = "l2", p_tol=1e-6, V_tol=1e-6))
  
  # method 3 (through spdlog + vectorform)
  T.l1.spdlog <- geo_dist(manifold="euclidean",intrinsic_location(manifold="euclidean", y.os.spdlog, w = wt.t.os, estimator = "l1", p_tol=1e-6, V_tol=1e-6),
                          intrinsic_location(manifold="euclidean", y.os.spdlog, w = wt.c.os, estimator = "l1", p_tol=1e-6, V_tol=1e-6))
  T.l2.spdlog <- geo_dist(manifold="euclidean",intrinsic_location(manifold="euclidean", y.os.spdlog, w = wt.t.os, estimator = "l2", p_tol=1e-6, V_tol=1e-6),
                          intrinsic_location(manifold="euclidean", y.os.spdlog, w = wt.c.os, estimator = "l2", p_tol=1e-6, V_tol=1e-6))
  
  # method 4 (through fa)
  T.l1.fa <- abs(weightedMedian(y.os.fa, wt.t.os) - weightedMedian(y.os.fa, wt.c.os))
  T.l2.fa <- abs(weightedMean(y.os.fa, wt.t.os) - weightedMean(y.os.fa, wt.c.os))
  
  # method 5 (through md)
  T.l1.md <- abs(weightedMedian(y.os.md, wt.t.os) - weightedMedian(y.os.md, wt.c.os))
  T.l2.md <- abs(weightedMean(y.os.md, wt.t.os) - weightedMean(y.os.md, wt.c.os))
  
  T.l1.collection <- c(T.l1.proposed, T.l1.vectorform, T.l1.spdlog,
                       T.l1.fa, T.l1.md)
  T.l2.collection <- c(T.l2.proposed, T.l2.vectorform, T.l2.spdlog,
                       T.l2.fa, T.l2.md)
  
  res.dist.l1.spd.os <- getpval_spd(z.os, group.use.os, y.os.list, T.l1.collection, B=500, estimator = 'l1', iter_total=j)
  res.dist.l2.spd.os <- getpval_spd(z.os, group.use.os, y.os.list, T.l2.collection, B=500, estimator = 'l2', iter_total=j)
  
  res.mat.l1.spd.os[j,] <- c(res.dist.l1.spd.os$rejection_region_proposed, res.dist.l1.spd.os$pval, T.l1.collection) 
  res.mat.l2.spd.os[j,] <- c(res.dist.l2.spd.os$rejection_region_proposed, res.dist.l2.spd.os$pval, T.l2.collection) 
  
  # do both randomized experiments and matched observational studies
  # assign strata for randomized experiments, and then assign treatment in the exact same way as in our existing experiments
  # test Fisher's sharp null (with 500 permutations) in 5 ways, comparing our method with 4 Euclidean alternatives. Use the same strata/weights each time:
  # 1. using our manifold-based test statistic using the spdL2center(), spdL1center() and spddistance() functions.
  # 2. put all observations through vectorform() to get 6D vectors, then do multivariate Euclidean version.
  # 3. put all observations through spdlog(diag(3), ...) and then vectorform() to get 6D vectors in tangent space at diag(3), then do multivariate Euclidean version.
  # 4. put all observations through fa() to get scalars, then do univariate Euclidean version.
  # 5. put all observations through md() to get scalars, then do univariate Euclidean version.
  # save p-values for each of these 5 tests for each j. also save the 95\% rejection region for our manifold-based test for each j.
}


colnames(res.mat.l1.spd.os) <- c("rejection_region", "proposed.pval","vectorform.pval",
                              "spdlog.pval", "fa.pval", "md.pval", "proposed.T","vectorform.T",
                              "spdlog.T", "fa.T", "md.T")
colnames(res.mat.l2.spd.os) <- c("rejection_region", "proposed.pval","vectorform.pval",
                              "spdlog.pval", "fa.pval", "md.pval", "proposed.T","vectorform.T",
                              "spdlog.T", "fa.T", "md.T")



# summary and visualization 
res.mat.l1.spd.sre <- as.data.frame(res.mat.l1.spd.sre)
res.mat.l2.spd.sre <- as.data.frame(res.mat.l2.spd.sre)
res.mat.l1.spd.os <- as.data.frame(res.mat.l1.spd.os)
res.mat.l2.spd.os <- as.data.frame(res.mat.l2.spd.os)

res_mat <- NULL
for(i in c(2,3,5,6)){
  tmp <- c(sum(res.mat.l2.spd.sre[,i]>0.1) / 500, sum(res.mat.l2.spd.sre[,i]>0.05) / 500,
           sum(res.mat.l2.spd.sre[,i]>0.01) / 500, sum(res.mat.l2.spd.sre[,i]>0.005) / 500,
           sum(res.mat.l2.spd.sre[,i]>0.001) / 500)
  res_mat <- rbind(res_mat, c(tmp, mean(res.mat.l2.spd.sre[,i])))
  
  res_mat <- rbind(res_mat, 
                   c(sqrt(tmp*(1-tmp)/499), sd(res.mat.l2.spd.sre[,i])))
}

for(i in c(2,3,5,6)){
  tmp <- c(sum(res.mat.l1.spd.sre[,i]>0.1) / 500, sum(res.mat.l1.spd.sre[,i]>0.05) / 500,
           sum(res.mat.l1.spd.sre[,i]>0.01) / 500, sum(res.mat.l1.spd.sre[,i]>0.005) / 500,
           sum(res.mat.l1.spd.sre[,i]>0.001) / 500)
  res_mat <- rbind(res_mat, c(tmp, mean(res.mat.l1.spd.sre[,i])))
  
  res_mat <- rbind(res_mat, 
                   c(sqrt(tmp*(1-tmp)/499), sd(res.mat.l1.spd.sre[,i])))
}

for(i in c(2,3,5,6)){
  tmp <- c(sum(res.mat.l2.spd.os[,i]>0.1) / 500, sum(res.mat.l2.spd.os[,i]>0.05) / 500,
           sum(res.mat.l2.spd.os[,i]>0.01) / 500, sum(res.mat.l2.spd.os[,i]>0.005) / 500,
           sum(res.mat.l2.spd.os[,i]>0.001) / 500)
  res_mat <- rbind(res_mat, c(tmp, mean(res.mat.l2.spd.os[,i])))
  
  res_mat <- rbind(res_mat, 
                   c(sqrt(tmp*(1-tmp)/499), sd(res.mat.l2.spd.os[,i])))
}


for(i in c(2,3,5,6)){
  tmp <- c(sum(res.mat.l1.spd.os[,i]>0.1) / 500, sum(res.mat.l1.spd.os[,i]>0.05) / 500,
           sum(res.mat.l1.spd.os[,i]>0.01) / 500, sum(res.mat.l1.spd.os[,i]>0.005) / 500,
           sum(res.mat.l1.spd.os[,i]>0.001) / 500)
  res_mat <- rbind(res_mat, c(tmp, mean(res.mat.l1.spd.os[,i])))
  
  res_mat <- rbind(res_mat, 
                   c(sqrt(tmp*(1-tmp)/499), sd(res.mat.l1.spd.os[,i])))
}


res_mat


res_mat2 <- rbind(c(min(res.mat.l2.spd.sre[,1]), as.numeric(quantile(res.mat.l2.spd.sre[,1], c(0.1*c(1:9)))), max(res.mat.l2.spd.sre[,1])),
                  c(min(res.mat.l1.spd.sre[,1]), as.numeric(quantile(res.mat.l1.spd.sre[,1], c(0.1*c(1:9)))), max(res.mat.l1.spd.sre[,1])),
                  c(min(res.mat.l2.spd.os[,1]), as.numeric(quantile(res.mat.l2.spd.os[,1], c(0.1*c(1:9)))), max(res.mat.l2.spd.os[,1])),
                  c(min(res.mat.l1.spd.os[,1]), as.numeric(quantile(res.mat.l1.spd.os[,1], c(0.1*c(1:9)))), max(res.mat.l1.spd.os[,1])))

res_mat2



## Combine rejection values
rejection_values <- list(
  "L2 SRE" = res.mat.l2.spd.sre[, 1],
  "L1 SRE" = res.mat.l1.spd.sre[, 1],
  "L2 OBS" = res.mat.l2.spd.os[, 1],
  "L1 OBS" = res.mat.l1.spd.os[, 1]
)

## Boxplot
## Colors for boxplots
cols <- c("steelblue3", "skyblue", "tomato", "orange")  # 4가지 완전 구분

par(mar = c(5, 5, 4, 2)) 

## Boxplot
boxplot(
  rejection_values,
  col = cols,
  las = 1,
  ylab = "Rejection value",
  cex.lab = 1.3,      # y label 글씨 키우기
  xaxt = "n",         # x축 label 제거
  main = ""
)

## Legend
legend(
  "topleft",
  legend = c(expression(alpha == 2~", Stratified randomized experiment"),
             expression(alpha == 1~", Stratified randomized experiment"),
             expression(alpha == 2~", Observational study"),
             expression(alpha == 1~", Observational study")),
  fill = cols,
  bty = "n",
  y.intersp = 1.6  # 줄간격 1.4배
)
