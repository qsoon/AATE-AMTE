library(MASS)
library(R.matlab)
library(rmatio)
# library(GeodRegr)
library(optmatch)
library(senstrat)
library(abind)
library(ggplot2)

##################################################
## Simulations: consistency and boostrap C.I check
##################################################

setwd('/home/kyu9510/AATE-AMTE')

source("efficiency.R")
source("GeodRegr.R")
source("normal.R")
source("others.R")


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

##########################
## SRE
##########################

## initializations

manifold <- 'sphere'
# manifold <- 'hyperbolic'

dim <- 2
embed <- dim+1 # vector length
k <- 2 # number of covariates

p <- numeric(embed)
p[1] <- 1 # p = (1,0,0)
v <- matrix(0L, nrow = embed, ncol = k)
v[2, 1] <- pi / 4
v[embed, 2] <- -pi  /6
# sigma <- pi / 8
sigma <- pi / 8


## consistency check
N.list <- c(32,64,128,256,512,1024)
iter.num <- 1000
error.list.l2.sre <- c()
error.list.l1.sre <- c()

set.seed(1)
for(N in N.list){
  e.l2 <- 0
  e.l1 <- 0
  
  j <- 1
  while(j <= iter.num){
    print(paste("SRE","N:", N, "iter num:",j))
    
    x <- matrix(runif(k * N), ncol = k) - 0.5
    r_t <- matrix(0L, nrow = embed, ncol = N)
    r_c <- matrix(0L, nrow = embed, ncol = N)
    sims <- rnormtangents(manifold, N = N, n = dim, sigma_sq = sigma ^ 2)
    shifts <- v %*% t(x)
    mu_t <- exp_map(manifold, p, c(0,1,0))
    mu_c <- exp_map(manifold, p, -c(0,1,0))
    
    for (i in 1:N) {
      new_shift <- par_trans(manifold, p, mu_t, shifts[, i])
      r_t[, i] <- exp_map(manifold, mu_t, new_shift)
      error <- par_trans(manifold, p, r_t[, i], sims[, i])
      r_t[, i] <- exp_map(manifold, r_t[ ,i], error)
      new_shift <- par_trans(manifold, p, mu_c, shifts[, i])
      r_c[, i] <- exp_map(manifold, mu_c, new_shift)
      error <- par_trans(manifold, p, r_c[, i], sims[, i])
      r_c[, i] <- exp_map(manifold, r_c[ ,i], error)
    }
    
    strata.sre <- ifelse(x[,1] >= 0, 1, 2) # make 2 strata based on covariate x1
    strata.sre1.ind <- which(strata.sre==1)
    len.s1 <- length(strata.sre1.ind)
    strata.sre2.ind <- which(strata.sre==2)
    len.s2 <- length(strata.sre2.ind)
    
    z.sre <- rep(0,N) # independent of (r_t, r_c) given S, here, 50:50 in each stratum
    z.sre[sort(sample(strata.sre1.ind, floor((len.s1+1)/2), replace=FALSE))] <- 1 
    z.sre[sort(sample(strata.sre2.ind, floor((len.s2+1)/2), replace=FALSE))] <- 1
    
    wt.t.sre <- calculateWeight.T.sre(z.sre, strata.sre)
    wt.c.sre <- calculateWeight.C.sre(z.sre, strata.sre)
    if((sum(is.na(wt.t.sre))!=0) | (sum(wt.t.sre)!=1)){
      next
    }
    if((sum(is.na(wt.c.sre))!=0) | (sum(wt.c.sre)!=1)){
      next
    }
    
    y.sre <- matrix(0L, nrow = embed, ncol=N)
    y.sre[,which(z.sre==1)] <- r_t[,which(z.sre==1)]
    y.sre[,which(z.sre==0)] <- r_c[,which(z.sre==0)]
    
    estimator <- 'l2'
    T1.l2.sre <- geo_dist(manifold,intrinsic_location(manifold, y.sre, w = wt.t.sre, estimator, p_tol=1e-6, V_tol=1e-6),
                          intrinsic_location(manifold, y.sre, w = wt.c.sre, estimator, p_tol=1e-6, V_tol=1e-6))
    
    # T2.res.l2.sre <- geo_reg(manifold, z.sre, y.sre, 
    #                            w = calculateWeight(z.sre, strata.sre, sum(z.sre)/length(z.sre), 1-sum(z.sre)/length(z.sre)), 
    #                            estimator, max_iter = 200)
    # 
    # T2.l2.sre <- sqrt(abs(sum(Conj(T2.res.l2.sre$V)*T2.res.l2.sre$V)))
    
    
    estimator <- 'l1'
    T1.l1.sre <- geo_dist(manifold,intrinsic_location(manifold, y.sre, w = wt.t.sre, estimator, p_tol=1e-6, V_tol=1e-6),
                          intrinsic_location(manifold, y.sre, w = wt.c.sre, estimator, p_tol=1e-6, V_tol=1e-6))
    
    # T2.res.l1.sre <- geo_reg(manifold, z.sre, y.sre, 
    #                          w = calculateWeight(z.sre, strata.sre, sum(z.sre)/length(z.sre), 1-sum(z.sre)/length(z.sre)), 
    #                          estimator, max_iter = 200)
    # 
    # T2.l1.sre <- sqrt(abs(sum(Conj(T2.res.l1.sre$V)*T2.res.l1.sre$V)))
    
    
    e.l2 <- e.l2 + abs(2-T1.l2.sre)
    e.l1 <- e.l1 + abs(2-T1.l1.sre)
    
    # e.l2 <- e.l2 + abs(2-T2.l2.sre)
    # e.l1 <- e.l1 + abs(2-T2.l1.sre)
    
    j <- j + 1
  }
  error.list.l2.sre <- c(error.list.l2.sre, e.l2/iter.num)
  error.list.l1.sre <- c(error.list.l1.sre, e.l1/iter.num)
  
}

print(error.list.l2.sre)
print(error.list.l1.sre)


## C.I check
N.list <- c(32,64,128,256,512,1024)
iter.num <- 500
M <- 500
ci.list.l2.sre <- c()
ci.list.l1.sre <- c()

set.seed(1)
for(N in N.list){
  ci.cnt.l2 <- 0
  ci.cnt.l1 <- 0
  
  j <- 1
  while(j <= iter.num){
    print(paste("SRE","N:", N, "iter num:",j))
    
    x <- matrix(runif(k * N), ncol = k) - 0.5
    r_t <- matrix(0L, nrow = embed, ncol = N)
    r_c <- matrix(0L, nrow = embed, ncol = N)
    sims <- rnormtangents(manifold, N = N, n = dim, sigma_sq = sigma ^ 2)
    shifts <- v %*% t(x)
    mu_t <- exp_map(manifold, p, c(0,1,0))
    mu_c <- exp_map(manifold, p, -c(0,1,0))
    
    for (i in 1:N) {
      new_shift <- par_trans(manifold, p, mu_t, shifts[, i])
      r_t[, i] <- exp_map(manifold, mu_t, new_shift)
      error <- par_trans(manifold, p, r_t[, i], sims[, i])
      r_t[, i] <- exp_map(manifold, r_t[ ,i], error)
      new_shift <- par_trans(manifold, p, mu_c, shifts[, i])
      r_c[, i] <- exp_map(manifold, mu_c, new_shift)
      error <- par_trans(manifold, p, r_c[, i], sims[, i])
      r_c[, i] <- exp_map(manifold, r_c[ ,i], error)
    }
    
    strata.sre <- ifelse(x[,1] >= 0, 1, 2) # make 2 strata based on covariate x1
    strata.sre1.ind <- which(strata.sre==1)
    len.s1 <- length(strata.sre1.ind)
    strata.sre2.ind <- which(strata.sre==2)
    len.s2 <- length(strata.sre2.ind)
    
    z.sre <- rep(0,N) # independent of (r_t, r_c) given S, here, 50:50 in each stratum
    z.sre[sort(sample(strata.sre1.ind, floor((len.s1+1)/2), replace=FALSE))] <- 1 
    z.sre[sort(sample(strata.sre2.ind, floor((len.s2+1)/2), replace=FALSE))] <- 1
    
    wt.t.sre <- calculateWeight.T.sre(z.sre, strata.sre)
    wt.c.sre <- calculateWeight.C.sre(z.sre, strata.sre)
    if((sum(is.na(wt.t.sre))!=0) | (sum(wt.t.sre)!=1)){
      next
    }
    if((sum(is.na(wt.c.sre))!=0) | (sum(wt.c.sre)!=1)){
      next
    }
    
    y.sre <- matrix(0L, nrow = embed, ncol=N)
    y.sre[,which(z.sre==1)] <- r_t[,which(z.sre==1)]
    y.sre[,which(z.sre==0)] <- r_c[,which(z.sre==0)]
    
    estimator <- 'l2'
    T1.l2.sre <- geo_dist(manifold,intrinsic_location(manifold, y.sre, w = wt.t.sre, estimator, p_tol=1e-6, V_tol=1e-6),
                          intrinsic_location(manifold, y.sre, w = wt.c.sre, estimator, p_tol=1e-6, V_tol=1e-6))
    
    # T2.res.l2.sre <- geo_reg(manifold, z.sre, y.sre, 
    #                            w = calculateWeight(z.sre, strata.sre, sum(z.sre)/length(z.sre), 1-sum(z.sre)/length(z.sre)), 
    #                            estimator, max_iter = 200)
    # 
    # T2.l2.sre <- sqrt(abs(sum(Conj(T2.res.l2.sre$V)*T2.res.l2.sre$V)))
    
    
    estimator <- 'l1'
    T1.l1.sre <- geo_dist(manifold,intrinsic_location(manifold, y.sre, w = wt.t.sre, estimator, p_tol=1e-6, V_tol=1e-6),
                          intrinsic_location(manifold, y.sre, w = wt.c.sre, estimator, p_tol=1e-6, V_tol=1e-6))
    
    # T2.res.l1.sre <- geo_reg(manifold, z.sre, y.sre, 
    #                          w = calculateWeight(z.sre, strata.sre, sum(z.sre)/length(z.sre), 1-sum(z.sre)/length(z.sre)), 
    #                          estimator, max_iter = 200)
    # 
    # T2.l1.sre <- sqrt(abs(sum(Conj(T2.res.l1.sre$V)*T2.res.l1.sre$V)))
    
    
    T1.l2.sre.bstrp.list <- c()
    T1.l1.sre.bstrp.list <- c()
    
    # T2.l2.sre.bstrp.list <- c()
    # T2.l1.sre.bstrp.list <- c()
    
    m <- 1
    while(m <= M){
      bstrp.sample.sre <- sort(sample(1:N, N, replace=TRUE))
      x.bstrp.sre <- x[bstrp.sample.sre,]
      
      strata.sre.bstrp <- strata.sre[bstrp.sample.sre]
      # strata.sre1.bstrp.ind <- which(strata.sre.bstrp==1)
      # len.s1.bstrp <- length(strata.sre1.bstrp.ind)
      # strata.sre2.bstrp.ind <- which(strata.sre.bstrp==2)
      # len.s2.bstrp <- length(strata.sre2.bstrp.ind)
      
      # z.sre.bstrp <- rep(0,N) # independent of (r_t, r_c) given S, here, 50:50 in each stratum
      # z.sre.bstrp[sort(sample(strata.sre1.bstrp.ind, floor((len.s1.bstrp+1)/2), replace=FALSE))] <- 1 
      # z.sre.bstrp[sort(sample(strata.sre2.bstrp.ind, floor((len.s2.bstrp+1)/2), replace=FALSE))] <- 1
      
      # y.sre.bstrp <- matrix(0L, nrow = embed, ncol=N)
      # y.sre.bstrp[,which(z.sre.bstrp==1)] <- r_t[,bstrp.sample.sre][,which(z.sre.bstrp==1)]
      # y.sre.bstrp[,which(z.sre.bstrp==0)] <- r_c[,bstrp.sample.sre][,which(z.sre.bstrp==0)]
      
      z.sre.bstrp <- z.sre[bstrp.sample.sre]
      
      if(length(unique(z.sre.bstrp))==1){
        next
      }
      
      y.sre.bstrp <- y.sre[,bstrp.sample.sre]
      
      wt.t.sre.bstrp <- calculateWeight.T.sre(z.sre.bstrp, strata.sre.bstrp)
      wt.c.sre.bstrp <- calculateWeight.C.sre(z.sre.bstrp, strata.sre.bstrp)
      if((sum(is.na(wt.t.sre.bstrp))!=0) | (sum(wt.t.sre.bstrp)!=1)){
        next
      }
      if((sum(is.na(wt.c.sre.bstrp))!=0) | (sum(wt.c.sre.bstrp)!=1)){
        next
      }
      
      estimator <- 'l2'
      T1.l2.sre.bstrp <- geo_dist(manifold,intrinsic_location(manifold, y.sre.bstrp, w = wt.t.sre.bstrp, estimator, p_tol=1e-6, V_tol=1e-6),
                                  intrinsic_location(manifold, y.sre.bstrp, w = wt.c.sre.bstrp, estimator, p_tol=1e-6, V_tol=1e-6))
      
      T1.l2.sre.bstrp.list <- c(T1.l2.sre.bstrp.list, T1.l2.sre.bstrp)
      
      
      # T2.res.l2.sre.bstrp <- geo_reg(manifold, z.sre.bstrp, y.sre[,bstrp.sample.sre], 
      #                          w = calculateWeight(z.sre.bstrp, strata.sre.bstrp, sum(z.sre.bstrp)/length(z.sre.bstrp), 1-sum(z.sre.bstrp)/length(z.sre.bstrp)), 
      #                          estimator, max_iter = 200)
      # 
      # T2.l2.sre.bstrp <- sqrt(abs(sum(Conj(T2.res.l2.sre.bstrp$V)*T2.res.l2.sre.bstrp$V)))
      # 
      # T2.l2.sre.bstrp.list <- c(T2.l2.sre.bstrp.list, T2.l2.sre.bstrp)
      
      
      estimator <- 'l1'
      T1.l1.sre.bstrp <- geo_dist(manifold,intrinsic_location(manifold, y.sre.bstrp, w = wt.t.sre.bstrp, estimator, p_tol=1e-6, V_tol=1e-6),
                                  intrinsic_location(manifold, y.sre.bstrp, w = wt.c.sre.bstrp, estimator, p_tol=1e-6, V_tol=1e-6))
      
      T1.l1.sre.bstrp.list <- c(T1.l1.sre.bstrp.list, T1.l1.sre.bstrp)
      
      # T2.res.l1.sre.bstrp <- geo_reg(manifold, z.sre.bstrp, y.sre[,bstrp.sample.sre], 
      #                                w = calculateWeight(z.sre.bstrp, strata.sre.bstrp, sum(z.sre.bstrp)/length(z.sre.bstrp), 1-sum(z.sre.bstrp)/length(z.sre.bstrp)), 
      #                                estimator, max_iter = 200)
      # 
      # T2.l1.sre.bstrp <- sqrt(abs(sum(Conj(T2.res.l1.sre.bstrp$V)*T2.res.l1.sre.bstrp$V)))
      # 
      # T2.l1.sre.bstrp.list <- c(T2.l1.sre.bstrp.list, T2.l1.sre.bstrp)
      
      m <- m + 1
    }
    
    bstrp.pivot.ci.l2.sre <- as.vector(c(2*T1.l2.sre - quantile(T1.l2.sre.bstrp.list, 0.975), 2*T1.l2.sre - quantile(T1.l2.sre.bstrp.list, 0.025)))
    bstrp.pivot.ci.l1.sre <- as.vector(c(2*T1.l1.sre - quantile(T1.l1.sre.bstrp.list, 0.975), 2*T1.l1.sre - quantile(T1.l1.sre.bstrp.list, 0.025)))
    
    # bstrp.pivot.ci.l2.sre <- as.vector(c(2*T2.l2.sre - quantile(T2.l2.sre.bstrp.list, 0.975), 2*T2.l2.sre - quantile(T2.l2.sre.bstrp.list, 0.025)))
    # bstrp.pivot.ci.l1.sre <- as.vector(c(2*T2.l1.sre - quantile(T2.l1.sre.bstrp.list, 0.975), 2*T2.l1.sre - quantile(T2.l1.sre.bstrp.list, 0.025)))
    
    
    if(bstrp.pivot.ci.l2.sre[1]<=2 & bstrp.pivot.ci.l2.sre[2]>=2){
      ci.cnt.l2 <- ci.cnt.l2 + 1
    }
    if(bstrp.pivot.ci.l1.sre[1]<=2 & bstrp.pivot.ci.l1.sre[2]>=2){
      ci.cnt.l1 <- ci.cnt.l1 + 1
    }
    
    print(paste("L2",ci.cnt.l2/j, "L1", ci.cnt.l1/j))
    j <- j + 1
  }
  
  ci.list.l2.sre <- c(ci.list.l2.sre, ci.cnt.l2/iter.num)
  ci.list.l1.sre <- c(ci.list.l1.sre, ci.cnt.l1/iter.num)
}
print(ci.list.l2.sre)
print(ci.list.l1.sre)



##########################
## Observational study
##########################

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




## initializations

# manifold <- 'sphere'
manifold <- 'hyperbolic'

dim <- 2
embed <- dim+1 # vector length
k <- 2 # number of covariates

p <- numeric(embed)
p[1] <- 1 # p = (1,0,0)
v <- matrix(0L, nrow = embed, ncol = k)
v[2, 1] <- pi / 4
v[embed, 2] <- -pi  /6
# sigma <- pi / 8
sigma <- pi / 8


## consistency check
N.list <- c(32,64,128,256,512,1024)
iter.num <- 1000
error.list.l2.os <- c()
error.list.l1.os <- c()

set.seed(1)
for(N in N.list){
  e.l2 <- 0
  e.l1 <- 0
  
  j <- 1
  while(j <= iter.num){
    print(paste("Obs.study","N:", N, "iter num:",j))
    
    x <- matrix(runif(k * N), ncol = k) - 0.5
    r_t <- matrix(0L, nrow = embed, ncol = N)
    r_c <- matrix(0L, nrow = embed, ncol = N)
    sims <- rnormtangents(manifold, N = N, n = dim, sigma_sq = sigma ^ 2)
    shifts <- v %*% t(x)
    mu_t <- exp_map(manifold, p, c(0,1,0))
    mu_c <- exp_map(manifold, p, -c(0,1,0))
    
    for (i in 1:N) {
      new_shift <- par_trans(manifold, p, mu_t, shifts[, i])
      r_t[, i] <- exp_map(manifold, mu_t, new_shift)
      error <- par_trans(manifold, p, r_t[, i], sims[, i])
      r_t[, i] <- exp_map(manifold, r_t[ ,i], error)
      new_shift <- par_trans(manifold, p, mu_c, shifts[, i])
      r_c[, i] <- exp_map(manifold, mu_c, new_shift)
      error <- par_trans(manifold, p, r_c[, i], sims[, i])
      r_c[, i] <- exp_map(manifold, r_c[ ,i], error)
    }
    
    
    z.os <- rbinom(N, 1, 1 / (1+exp(-(x[,1]+x[,2])))) # dependent on (r_t, r_c)
    
    y.os <- matrix(0L, nrow = embed, ncol=N)
    y.os[,which(z.os==1)] <- r_t[,which(z.os==1)]
    y.os[,which(z.os==0)] <- r_c[,which(z.os==0)]
    
    x_data.os <- data.frame(V1 = x[,1], V2 = x[,2], Z = z.os)
    
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
    if((sum(is.na(wt.t.os))!=0) | (sum(wt.t.os)!=1)){
      next
    }
    if((sum(is.na(wt.c.os))!=0) | (sum(wt.c.os)!=1)){
      next
    }
    
    estimator <- 'l2'
    T1.l2.os <- geo_dist(manifold,intrinsic_location(manifold, y.os, w = wt.t.os, estimator, p_tol=1e-6, V_tol=1e-6),
                         intrinsic_location(manifold, y.os, w = wt.c.os, estimator, p_tol=1e-6, V_tol=1e-6))
    
    # T2.res.l2.os <- geo_reg(manifold, z.os, y.os, 
    #                            w = calculateWeight(z.os, strata.os, sum(z.os)/length(z.os), 1-sum(z.os)/length(z.os)), 
    #                            estimator, max_iter = 200)
    # 
    # T2.l2.os <- sqrt(abs(sum(Conj(T2.res.l2.os$V)*T2.res.l2.os$V)))
    
    
    estimator <- 'l1'
    T1.l1.os <- geo_dist(manifold,intrinsic_location(manifold, y.os, w = wt.t.os, estimator, p_tol=1e-6, V_tol=1e-6),
                         intrinsic_location(manifold, y.os, w = wt.c.os, estimator, p_tol=1e-6, V_tol=1e-6))
    
    # T2.res.l1.os <- geo_reg(manifold, z.os, y.os, 
    #                          w = calculateWeight(z.os, strata.os, sum(z.os)/length(z.os), 1-sum(z.os)/length(z.os)), 
    #                          estimator, max_iter = 200)
    # 
    # T2.l1.os <- sqrt(abs(sum(Conj(T2.res.l1.os$V)*T2.res.l1.os$V)))
    
    
    e.l2 <- e.l2 + abs(2-T1.l2.os)
    e.l1 <- e.l1 + abs(2-T1.l1.os)
    
    # e.l2 <- e.l2 + abs(2-T2.l2.os)
    # e.l1 <- e.l1 + abs(2-T2.l1.os)
    
    j <- j + 1
  }
  error.list.l2.os <- c(error.list.l2.os, e.l2/iter.num)
  error.list.l1.os <- c(error.list.l1.os, e.l1/iter.num)
  
}

print(error.list.l2.os)
print(error.list.l1.os)


## C.I check
N.list <- c(32,64,128,256,512,1024)
iter.num <- 500
M <- 500
ci.list.l2.os <- c()
ci.list.l1.os <- c()

set.seed(1)
for(N in N.list){
  ci.cnt.l2 <- 0
  ci.cnt.l1 <- 0
  
  j <- 1
  while(j <= iter.num){
    print(paste("Obs.study","N:", N, "iter num:",j))
    
    x <- matrix(runif(k * N), ncol = k) - 0.5
    r_t <- matrix(0L, nrow = embed, ncol = N)
    r_c <- matrix(0L, nrow = embed, ncol = N)
    sims <- rnormtangents(manifold, N = N, n = dim, sigma_sq = sigma ^ 2)
    shifts <- v %*% t(x)
    mu_t <- exp_map(manifold, p, c(0,1,0))
    mu_c <- exp_map(manifold, p, -c(0,1,0))
    
    for (i in 1:N) {
      new_shift <- par_trans(manifold, p, mu_t, shifts[, i])
      r_t[, i] <- exp_map(manifold, mu_t, new_shift)
      error <- par_trans(manifold, p, r_t[, i], sims[, i])
      r_t[, i] <- exp_map(manifold, r_t[ ,i], error)
      new_shift <- par_trans(manifold, p, mu_c, shifts[, i])
      r_c[, i] <- exp_map(manifold, mu_c, new_shift)
      error <- par_trans(manifold, p, r_c[, i], sims[, i])
      r_c[, i] <- exp_map(manifold, r_c[ ,i], error)
    }
    
    
    z.os <- rbinom(N, 1, 1 / (1+exp(-(x[,1]+x[,2])))) # dependent on (r_t, r_c)
    
    y.os <- matrix(0L, nrow = embed, ncol=N)
    y.os[,which(z.os==1)] <- r_t[,which(z.os==1)]
    y.os[,which(z.os==0)] <- r_c[,which(z.os==0)]
    
    x_data.os <- data.frame(V1 = x[,1], V2 = x[,2], Z = z.os)
    
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
    if((sum(is.na(wt.t.os))!=0) | (sum(wt.t.os)!=1)){
      next
    }
    if((sum(is.na(wt.c.os))!=0) | (sum(wt.c.os)!=1)){
      next
    }
    
    estimator <- 'l2'
    T1.l2.os <- geo_dist(manifold,intrinsic_location(manifold, y.os, w = wt.t.os, estimator, p_tol=1e-6, V_tol=1e-6),
                         intrinsic_location(manifold, y.os, w = wt.c.os, estimator, p_tol=1e-6, V_tol=1e-6))
    
    # T2.res.l2.os <- geo_reg(manifold, z.os, y.os, 
    #                            w = calculateWeight(z.os, strata.os, sum(z.os)/length(z.os), 1-sum(z.os)/length(z.os)), 
    #                            estimator, max_iter = 200)
    # 
    # T2.l2.os <- sqrt(abs(sum(Conj(T2.res.l2.os$V)*T2.res.l2.os$V)))
    
    
    estimator <- 'l1'
    T1.l1.os <- geo_dist(manifold,intrinsic_location(manifold, y.os, w = wt.t.os, estimator, p_tol=1e-6, V_tol=1e-6),
                         intrinsic_location(manifold, y.os, w = wt.c.os, estimator, p_tol=1e-6, V_tol=1e-6))
    
    # T2.res.l1.os <- geo_reg(manifold, z.os, y.os, 
    #                          w = calculateWeight(z.os, strata.os, sum(z.os)/length(z.os), 1-sum(z.os)/length(z.os)), 
    #                          estimator, max_iter = 200)
    # 
    # T2.l1.os <- sqrt(abs(sum(Conj(T2.res.l1.os$V)*T2.res.l1.os$V)))
    
    
    T1.l2.os.bstrp.list <- c()
    T1.l1.os.bstrp.list <- c()
    
    # T2.l2.os.bstrp.list <- c()
    # T2.l1.os.bstrp.list <- c()
    
    m <- 1
    while(m <= M){
      bstrp.sample.os <- sort(sample(1:N, N, replace=TRUE))
      x.bstrp.os <- x[bstrp.sample.os,]
      z.os.bstrp <- z.os[bstrp.sample.os]
      if(length(unique(z.os.bstrp))==1){
        next
      }
      y.os.bstrp <- y.os[,bstrp.sample.os]
      
      x_data.os.bstrp <- data.frame(V1 = x.bstrp.os[,1], V2 = x.bstrp.os[,2], Z = z.os.bstrp)
      propscore.model.os.bstrp <- glm(Z ~ V1 + V2,
                                      family=binomial, x=TRUE, data=x_data.os.bstrp)
      
      Xmat.all.os.bstrp <- (propscore.model.os.bstrp$x[,-1])
      
      # Rank based Mahalanobis distance
      distmat.all.os.bstrp <- rankmahal(z.os.bstrp, Xmat.all.os.bstrp)
      # Add caliper
      logit.propscore.all.os.bstrp <- (predict(propscore.model.os.bstrp))
      distmat2.all.os.bstrp <- addcaliper(distmat.all.os.bstrp, z.os.bstrp, logit.propscore.all.os.bstrp) 
      
      ## full matching 
      matchvec.full.all.os.bstrp <- fullmatch(distmat2.all.os.bstrp)
      
      treated.subject.index.full.all.os.bstrp <- rep(0,sum(z.os.bstrp==1))
      # The subject indices in the order of matchvec
      matchedset.index.full.all.os.bstrp <- substr(matchvec.full.all.os.bstrp, start=3, stop=10) # group name
      matchedset.index.full.numeric.all.os.bstrp <- as.numeric(matchedset.index.full.all.os.bstrp) # group name (numeric)
      subjects.match.order.full.all.os.bstrp <- as.numeric(names(matchvec.full.all.os.bstrp)) # subject indices
      
      # Create a numeric variable for which stratum each unit belongs to
      # 0 denotes that the unit was not matched
      # there is not units not matched b/c we use full matching
      stratum.short.full.all.os.bstrp <- substr(matchvec.full.all.os.bstrp, start=3, stop=10)
      stratum.full.numeric.all.os.bstrp <- as.numeric(stratum.short.full.all.os.bstrp) # i-th element : group name of unit i
      
      # Reassign numbers to each stratum that go from 1,..., no. of straum
      sort.unique.stratum.full.all.os.bstrp <- sort(unique(stratum.full.numeric.all.os.bstrp)) 
      stratum.myindex.matchvecorder.full.all.os.bstrp <- rep(0,length(stratum.full.numeric.all.os.bstrp))
      for(ii in 1:length(sort.unique.stratum.full.all.os.bstrp)){
        stratum.myindex.matchvecorder.full.all.os.bstrp[stratum.full.numeric.all.os.bstrp==sort.unique.stratum.full.all.os.bstrp[ii]] <- ii  
      }
      
      stratum.myindex.full.all.os.bstrp <- rep(0, length(stratum.myindex.matchvecorder.full.all.os.bstrp)) 
      stratum.myindex.full.all.os.bstrp[subjects.match.order.full.all.os.bstrp] <- stratum.myindex.matchvecorder.full.all.os.bstrp # i-th element : matched group index for unit i 
      
      group.use.os.bstrp <- stratum.myindex.full.all.os.bstrp
      
      wt.t.os.bstrp <- calculateWeight.T(z.os.bstrp, group.use.os.bstrp)
      wt.c.os.bstrp <- calculateWeight.C(z.os.bstrp, group.use.os.bstrp)
      if((sum(is.na(wt.t.os.bstrp))!=0) | (sum(wt.t.os.bstrp)!=1)){
        next
      }
      if((sum(is.na(wt.c.os.bstrp))!=0) | (sum(wt.c.os.bstrp)!=1)){
        next
      }
      
      estimator <- 'l2'
      T1.l2.os.bstrp <- geo_dist(manifold,intrinsic_location(manifold, y.os.bstrp, w = wt.t.os.bstrp, estimator, p_tol=1e-6, V_tol=1e-6),
                                 intrinsic_location(manifold, y.os.bstrp, w = wt.c.os.bstrp, estimator, p_tol=1e-6, V_tol=1e-6))
      
      T1.l2.os.bstrp.list <- c(T1.l2.os.bstrp.list, T1.l2.os.bstrp)
      
      
      # T2.res.l2.os.bstrp <- geo_reg(manifold, z.os.bstrp, y.os[,bstrp.sample.os], 
      #                          w = calculateWeight(z.os.bstrp, strata.os.bstrp, sum(z.os.bstrp)/length(z.os.bstrp), 1-sum(z.os.bstrp)/length(z.os.bstrp)), 
      #                          estimator, max_iter = 200)
      # 
      # T2.l2.os.bstrp <- sqrt(abs(sum(Conj(T2.res.l2.os.bstrp$V)*T2.res.l2.os.bstrp$V)))
      # 
      # T2.l2.os.bstrp.list <- c(T2.l2.os.bstrp.list, T2.l2.os.bstrp)
      
      
      estimator <- 'l1'
      T1.l1.os.bstrp <- geo_dist(manifold,intrinsic_location(manifold, y.os.bstrp, w = wt.t.os.bstrp, estimator, p_tol=1e-6, V_tol=1e-6),
                                 intrinsic_location(manifold, y.os.bstrp, w = wt.c.os.bstrp, estimator, p_tol=1e-6, V_tol=1e-6))
      
      T1.l1.os.bstrp.list <- c(T1.l1.os.bstrp.list, T1.l1.os.bstrp)
      
      # T2.res.l1.os.bstrp <- geo_reg(manifold, z.os.bstrp, y.os[,bstrp.sample.os], 
      #                                w = calculateWeight(z.os.bstrp, strata.os.bstrp, sum(z.os.bstrp)/length(z.os.bstrp), 1-sum(z.os.bstrp)/length(z.os.bstrp)), 
      #                                estimator, max_iter = 200)
      # 
      # T2.l1.os.bstrp <- sqrt(abs(sum(Conj(T2.res.l1.os.bstrp$V)*T2.res.l1.os.bstrp$V)))
      # 
      # T2.l1.os.bstrp.list <- c(T2.l1.os.bstrp.list, T2.l1.os.bstrp)
      
      m <- m + 1
    }
    
    bstrp.pivot.ci.l2.os <- as.vector(c(2*T1.l2.os - quantile(T1.l2.os.bstrp.list, 0.975), 2*T1.l2.os - quantile(T1.l2.os.bstrp.list, 0.025)))
    bstrp.pivot.ci.l1.os <- as.vector(c(2*T1.l1.os - quantile(T1.l1.os.bstrp.list, 0.975), 2*T1.l1.os - quantile(T1.l1.os.bstrp.list, 0.025)))
    
    # bstrp.pivot.ci.l2.os <- as.vector(c(2*T2.l2.os - quantile(T2.l2.os.bstrp.list, 0.975), 2*T2.l2.os - quantile(T2.l2.os.bstrp.list, 0.025)))
    # bstrp.pivot.ci.l1.os <- as.vector(c(2*T2.l1.os - quantile(T2.l1.os.bstrp.list, 0.975), 2*T2.l1.os - quantile(T2.l1.os.bstrp.list, 0.025)))
    
    
    if(bstrp.pivot.ci.l2.os[1]<=2 & bstrp.pivot.ci.l2.os[2]>=2){
      ci.cnt.l2 <- ci.cnt.l2 + 1
    }
    if(bstrp.pivot.ci.l1.os[1]<=2 & bstrp.pivot.ci.l1.os[2]>=2){
      ci.cnt.l1 <- ci.cnt.l1 + 1
    }
    
    print(paste("L2",ci.cnt.l2/j, "L1", ci.cnt.l1/j))
    j <- j + 1
  }
  
  ci.list.l2.os <- c(ci.list.l2.os, ci.cnt.l2/iter.num)
  ci.list.l1.os <- c(ci.list.l1.os, ci.cnt.l1/iter.num)
}

print(ci.list.l2.os)
print(ci.list.l1.os)

