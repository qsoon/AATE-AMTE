# Inner product. v1 and v2 should either be matrices of the same size or one of them (v1 in the hyperbolic case) should be a matrix and the other a vector of length equal to the number of rows in the matrix.
ip <- function(manifold, v1, v2) {
  if ((manifold == 'euclidean') | (manifold == 'sphere')) {
    result <- colSums(v1 * v2)
  } else if (manifold == 'hyperbolic') {
    v1[1, ] <- -v1[1, ]
    result <- colSums(v1 * v2)
  } else if (manifold == 'kendall') {
    result <- colSums(v1 * Conj(v2))
  }
  return(result)
}

# Magnitude of a vector. v should be a matrix.
mag <- function(manifold, v) {
  return(Re(sqrt(ip(manifold, v, v) + 0i)))
}

# Function for internal use in the geo_reg function. Similar to exp_map, but vectorized and without the errors, for speed and so as to not unexpectedly stop the algorithm because of one small error that can be safely ignored without changing the result. Difference from expo2: p and v should be matrices of the same dimensions.
expo <- function(manifold, p, v) {
  if (manifold == 'euclidean') {
    result <- p + v
  } else if (manifold == 'sphere') {
    theta <- mag(manifold, v)
    e1 <- t(t(p) / mag(manifold, p)) # reprojects p onto the manifold, for precision
    e2 <- t(t(v) / theta)
    result <- t(t(e1) * cos(theta)) + t(t(e2) * sin(theta))
    index <- which(theta == 0) # theta == 0 case must be dealt with separately due to division by theta
    result[, index] <- p[, index]
    result <- t(t(result) / mag(manifold, result)) # reprojects result onto the manifold, for precision
  } else if (manifold == 'hyperbolic') {
    theta <- mag(manifold, v)
    e1 <- t(t(p) / sqrt(-ip(manifold, p, p))) # reprojects p onto the manifold, for precision
    e2 <- t(t(v) / theta)
    result <- t(t(e1) * cosh(theta)) + t(t(e2) * sinh(theta))
    index <- which(theta == 0) # theta == 0 case must be dealt with separately due to division by theta
    result[, index] <- p[, index] # if theta == 0, the result is p
    result <- t(t(result) / sqrt(-ip(manifold, result, result))) # reprojects result onto the manifold, for precision
  } else if (manifold == 'kendall') {
    meanp <- colMeans(p)
    theta <- mag(manifold, v)
    e1 <- t((t(p) - meanp) / mag(manifold, t(t(p) - meanp))) # reprojects p onto the manifold, for precision
    e2 <- t(t(v) / theta)
    result <- t(t(e1) * cos(theta) + t(e2) * sin(theta))
    index <- which(theta == 0) # theta == 0 case must be dealt with separately due to division by theta
    result[, index] <- p[, index]
    result <- t((t(result) - colMeans(result)) / mag(manifold, t(t(result) - colMeans(result)))) # reprojects result onto the manifold, for precision
  }
  return(result)
}

# Function for internal use in the geo_reg function. Similar to exp_map, but vectorized and without the errors, for speed and so as to not unexpectedly stop the algorithm because of one small error that can be safely ignored without changing the result. Difference from expo: p should be a column matrix, v should be a matrix with the same number of rows as p.
expo2 <- function(manifold, p, v) {
  if (manifold == 'euclidean') {
    result <- as.vector(p) + v
  } else if (manifold == 'sphere') {
    theta <- mag(manifold, v)
    e1 <- t(t(p) / mag(manifold, p))
    e2 <- t(t(v) / theta)
    result <- e1 %*% cos(theta) + t(t(e2) * sin(theta))
    index <- which(theta == 0)
    result[, index] <- p
    result <- t(t(result) / mag(manifold, result))
  } else if (manifold == 'hyperbolic') {
    theta <- mag(manifold, v)
    e1 <- t(t(p) / sqrt(-ip(manifold, p, p)))
    e2 <- t(t(v) / theta)
    result <- e1 %*% cosh(theta) + t(t(e2) * sinh(theta))
    index <- which(theta == 0)
    result[, index] <- p
    result <- t(t(result) / sqrt(-ip(manifold, result, result)))
  } else if (manifold == 'kendall') {
    meanp <- colMeans(p)
    theta <- mag(manifold, v)
    e1 <- t((t(p) - meanp) / mag(manifold, t(t(p) - meanp)))
    e2 <- t(t(v) / theta)
    result <- e1 %*% t(cos(theta)) + t(t(e2) * sin(theta))
    index <- which(theta == 0)
    result[, index] <- p
    result <- t((t(result) - colMeans(result)) / mag(manifold, t(t(result) - colMeans(result))))
  }
  return(result)
}

# Function for internal use in the geo_reg function. Similar to log_map, but vectorized and without the errors, for speed and so as to not unexpectedly stop the algorithm because of one small error that can be safely ignored without changing the result. Difference from loga2 and loga3: p1 and p2 should be matrices of the same dimensions.
loga <- function(manifold, p1, p2) {
  if (manifold == 'euclidean') {
    result <- p2 - p1
  } else if (manifold == 'sphere') {
    p1 <- t(t(p1) / mag(manifold, p1)) # reprojects p1 onto the manifold, for precision
    p2 <- t(t(p2) / mag(manifold, p2)) # reprojects p2 onto the manifold, for precision
    a <- pmax(pmin(ip(manifold, p1, p2), 1), -1) # ensures a is in [-1, 1]
    theta <- acos(a)
    tang <- p2 - t(t(p1) * a)
    t <- mag(manifold, tang)
    if (any(t == 0)) { # any(t == 0) case must be dealt with separately due to division by t
      if (any(mag(manifold, p1 - p2) < 1e-6)) { # determining whether any(t == 0) because of p1 = p2 or p1 = -p2
        result <- numeric(dim(p1)[1])
      } else if (any(mag(manifold, p1 + p2) < 1e-6)) { # determining whether any(t == 0) because of p1 = p2 or p1 = -p2
        stop('p2 is the antipode of p1 and is therefore not in the domain of the log map at p1') ## change to continue
      }
    }
    result <- t(t(tang) * (theta / t))
    result[, which(t == 0)] <- 0
  } else if (manifold == 'hyperbolic') {
    p1 <- t(t(p1) / sqrt(-ip(manifold, p1, p1))) # reprojects p1 onto the manifold, for precision
    p2 <- t(t(p2) / sqrt(-ip(manifold, p2, p2))) # reprojects p2 onto the manifold, for precision
    a <- pmin(ip(manifold, p1, p2), -1) # ensures -a is at least 1
    theta <- acosh(-a)
    tang <- p2 + t(t(p1) * a)
    t <- mag(manifold, tang)
    if (any(t == 0)) { # any(t == 0) case must be dealt with separately due to division by t
      result <- numeric(dim(p1)[1])
    }
    result <- t(t(tang) * (theta / t))
    result[, which(t == 0)] <- 0
  } else if (manifold == 'kendall') {
    meanp1 <- colMeans(p1)
    meanp2 <- colMeans(p2)
    p1 <- t((t(p1) - meanp1) / mag(manifold, t(t(p1) - meanp1))) # reprojects p1 onto the manifold, for precision
    p2 <-  t((t(p2) - meanp2) / mag(manifold, t(t(p2) - meanp2))) # reprojects p2 onto the manifold, for precision
    a <- ip(manifold, p1, p2)
    theta <- acos(pmax(pmin(abs(a), 1), -1)) # ensures argument is in [-1, 1]
    tang <- t(t(p2) * (a / abs(a))) - t(t(p1) * abs(a))
    result <- t(t(tang) * (theta / mag(manifold, tang)))
    result[, which(mag(manifold, tang) == 0)] <- 0 # mag(manifold, tang) == 0 case must be dealt with separately due to division by mag(manifold, tang)
  }
  return(result)
}

# Function for internal use in the geo_reg function. Similar to log_map, but vectorized and without the errors, for speed and so as to not unexpectedly stop the algorithm because of one small error that can be safely ignored without changing the result. Difference from loga1 and loga3: p1 should be a column matrix, p2 should be a matrix with the same number of rows as p1. Needed for the sphere, hyperbolic space and Kendall's shape space.
loga2 <- function(manifold, p1, p2) {
  if (manifold == 'sphere') {
    p1 <- t(t(p1) / mag(manifold, p1)) # reprojects p1 onto the manifold, for precision
    p2 <- t(t(p2) / mag(manifold, p2)) # reprojects p2 onto the manifold, for precision
    a <- pmax(pmin(ip(manifold, p2, as.vector(p1)), 1), -1) # ensures a is in [-1, 1]
    theta <- acos(a)
    tang <- p2 - p1 %*% a
    t <- mag(manifold, tang)
    result <- t(t(tang) * (theta / t))
    result[, which(t == 0)] <- 0
  } else if (manifold == 'hyperbolic') {
    p1 <- t(t(p1) / sqrt(-ip(manifold, p1, p1))) # reprojects p1 onto the manifold, for precision
    p2 <- t(t(p2) / sqrt(-ip(manifold, p2, p2))) # reprojects p2 onto the manifold, for precision
    a <- pmin(ip(manifold, p2, as.vector(p1)), -1)
    theta <- acosh(-a)
    tang <- p2 + p1 %*% a
    t <- mag(manifold, tang)
    result <- t(t(tang) * (theta / t))
    result[, which(t == 0)] <- 0
  } else if (manifold == 'kendall') {
    meanp1 <- colMeans(p1)
    meanp2 <- colMeans(p2)
    p1 <- t((t(p1) - meanp1) / mag(manifold, t(t(p1) - meanp1))) # reprojects p1 onto the manifold, for precision
    p2 <-  t((t(p2) - meanp2) / mag(manifold, t(t(p2) - meanp2))) # reprojects p1 onto the manifold, for precision
    a <- ip(manifold, as.vector(p1), p2)
    theta <- acos(pmax(pmin(abs(a), 1), -1))
    tang <- t(t(p2) * (a / abs(a))) - p1 %*% abs(a)
    result <- t(t(tang) * (theta / mag(manifold, tang)))
    result[, which(mag(manifold, tang) == 0)] <- 0
  }
  return(result)
}

# Function for internal use in the geo_reg function. Similar to log_map, but vectorized and without the errors, for speed and so as to not unexpectedly stop the algorithm because of one small error that can be safely ignored without changing the result. Difference from loga2 and loga3: p1 should be a matrix, p2 should be a column matrix with the same number of rows as p1. Needed for the sphere and hyperbolic space.
loga3 <- function(manifold, p1, p2) {
  if (manifold == 'sphere') {
    p1 <- t(t(p1) / mag(manifold, p1)) # reprojects p1 onto the manifold, for precision
    p2 <- t(t(p2) / mag(manifold, p2)) # reprojects p2 onto the manifold, for precision
    a <- pmax(pmin(ip(manifold, p1, as.vector(p2)), 1), -1)
    theta <- acos(a)
    tang <- as.vector(p2) - t(t(p1) * a)
    t <- mag(manifold, tang)
    result <- t(t(tang) * (theta / t))
    result[, which(t == 0)] <- 0
  } else if (manifold == 'hyperbolic') {
    p1 <- t(t(p1) / sqrt(-ip(manifold, p1, p1))) # reprojects p1 onto the manifold, for precision
    p2 <- t(t(p2) / sqrt(-ip(manifold, p2, p2))) # reprojects p2 onto the manifold, for precision
    a <- pmin(ip(manifold, p1, as.vector(p2)), -1)
    theta <- acosh(-a)
    tang <- as.vector(p2) + t(t(p1) * a)
    t <- mag(manifold, tang)
    result <- t(t(tang) * (theta / t))
    result[, which(t == 0)] <- 0
  }
  return(result)
}

# Function for internal use in the geo_reg function. Similar to geo_dis, but vectorized and without the errors, for speed and so as to not unexpectedly stop the algorithm because of one small error that can be safely ignored without changing the result. p1 and p2 should be matrices of the same dimension.
dist <- function(manifold, p1, p2) {
  return(mag(manifold, loga(manifold, p1, p2)))
}

# Function for internal use in the geo_reg function. Similar to par_trans, but vectorized and without the errors, for speed and so as to not unexpectedly stop the algorithm because of one small error that can be safely ignored without changing the result. Difference from pt2 and pt3: p1, p2, and v should be matrices of the same dimensions.
pt <- function(manifold, p1, p2, v) {
  if (manifold == 'euclidean') {
    result <- v
  } else if (manifold == 'sphere') {
    p1 <- t(t(p1) / mag(manifold, p1)) # reprojects p1 onto the manifold, for precision
    w <- loga(manifold, p1, p2)
    t <- mag(manifold, w)
    e1 <- p1
    e2 <- t(t(w) / t)
    a <- ip(manifold, v, e2)
    invar <- v - t(t(e2) * a)
    result <- t(t(e2) * (a * cos(t))) - t(t(e1) * (a * sin(t))) + invar
    index <- which(t == 0) # t == 0 case must be dealt with separately due to division by t
    result[, index] <- v[, index]
    #    p1 <- t(t(p1) / mag(manifold, p1)) # reprojects p1 onto the manifold, for precision
    #    w <- loga(manifold, p1, p2)
    #    t <- mag(manifold, w)
    #    result <- v - t(t(w + loga(manifold, p2, p1)) * (ip(manifold, w, v) / (t ^ 2)))
    #    index <- which(t == 0) # t == 0 case must be dealt with separately due to division by t
    #    result[, index] <- v[, index]
  } else if (manifold == 'hyperbolic') {
    p1 <- t(t(p1) / sqrt(-ip(manifold, p1, p1))) # reprojects p1 onto the manifold, for precision
    w <- loga(manifold, p1, p2)
    t <- mag(manifold, w)
    e1 <- p1
    e2 <- t(t(w) / t)
    a <- ip(manifold, v, e2)
    invar <- v - t(t(e2) * a)
    result <- t(t(e2) * (a * cosh(t))) + t(t(e1) * (a * sinh(t))) + invar
    index <- which(t == 0) # t == 0 case must be dealt with separately due to division by t
    result[, index] <- v[, index]
    #p1 <- t(t(p1) / sqrt(-ip(manifold, p1, p1))) # reprojects p1 onto the manifold, for precision
    #w <- loga(manifold, p1, p2)
    #t <- mag(manifold, w)
    #result <- v - t(t(w + loga(manifold, p2, p1)) * (ip(manifold, w, v) / (t ^ 2)))
    #index <- which(t == 0) # t == 0 case must be dealt with separately due to division by t
    #result[, index] <- v[, index]
  } else  if (manifold == 'kendall') {
    meanp1 <- colMeans(p1)
    meanp2 <- colMeans(p2)
    p1 <- t((t(p1) - meanp1) / mag(manifold, t(t(p1) - meanp1))) # reprojects p1 onto the manifold, for precision
    p2 <- t((t(p2) - meanp2) / mag(manifold, t(t(p2) - meanp2))) # reprojects p2 onto the manifold, for precision
    yi <- expo(manifold, p1, v)
    a <- ip(manifold, p1, p2)
    p2 <- t(t(p2) * (a / abs(a))) # optimal alignment of p2 with p1
    b <- (1 - (abs(a)) ^ 2) ^ 0.5
    p2tilde <- t(t(p2 - t(t(p1) * abs(a))) / b)
    result <- v - t(t(p1) * (ip(manifold, v, p1))) - t(t(p2tilde) * (ip(manifold, v, p2tilde))) + t(t(p1) * ((abs(a)) * (ip(manifold, v, p1)) - b * (ip(manifold, v, p2tilde)))) + t(t(p2tilde) * (b * (ip(manifold, v, p1)) + (abs(a)) * (ip(manifold, v, p2tilde))))
    result <- t(t(result) * (Conj(a / abs(a))))
    index <- which(abs(a) >= 1)
    result[, index] <- loga(manifold, p2, yi)[, index]
  }
  return(result)
}

# Function for internal use in the geo_reg function. Similar to par_trans, but vectorized and without the errors, for speed and so as to not unexpectedly stop the algorithm because of one small error that can be safely ignored without changing the result. Difference from pt1 and pt3: p1 and p2 should be column matrices, v should be a matrix with the same number of rows as p1 and p2.
pt2 <- function(manifold, p1, p2, v) {
  if (manifold == 'euclidean') {
    result <- v
  } else if (manifold == 'sphere') {
    p1 <- t(t(p1) / mag(manifold, p1))
    w <- loga(manifold, p1, p2)
    t <- mag(manifold, w)
    e1 <- p1
    e2 <- w / t
    if (t == 0) {
      result <- v
    } else {
      a <- ip(manifold, v, as.vector(e2))
      invar <- v - e2 %*% a
      result <- e2 %*% (a * cos(t)) - e1 %*% (a * sin(t)) + invar
    }
    #    p1 <- t(t(p1) / mag(manifold, p1))
    #    w <- loga(manifold, p1, p2)
    #    t <- mag(manifold, w)
    #    if (t == 0) {
    #      result <- v
    #    } else {
    #      result <- v - (w + loga(manifold, p2, p1)) %*% (ip(manifold, v, as.vector(w)) / (t ^ 2))
    #    }
  } else if (manifold == 'hyperbolic') {
    p1 <- t(t(p1) / sqrt(-ip(manifold, p1, p1)))
    w <- loga(manifold, p1, p2)
    t <- mag(manifold, w)
    e1 <- p1
    e2 <- w / t
    if (t == 0) {
      result <- v
    } else {
      a <- ip(manifold, v, as.vector(e2))
      invar <- v - e2 %*% a
      result <- e2 %*% (a * cosh(t)) + e1 %*% (a * sinh(t)) + invar
    }
    #    p1 <- t(t(p1) / sqrt(-ip(manifold, p1, p1)))
    #    w <- loga(manifold, p1, p2)
    #    t <- mag(manifold, w)
    #    if (t == 0) {
    #      result <- v
    #    } else {
    #      result <- v - (w + loga(manifold, p2, p1)) %*% (ip(manifold, v, as.vector(w)) / (t ^ 2))
    #    }
  } else if (manifold == 'kendall') {
    meanp1 <- colMeans(p1)
    meanp2 <- colMeans(p2)
    p1 <- t((t(p1) - meanp1) / mag(manifold, t(t(p1) - meanp1)))
    p2 <- t((t(p2) - meanp2) / mag(manifold, t(t(p2) - meanp2)))
    yi <- expo2(manifold, p1, v)
    a <- ip(manifold, p1, p2)
    if (abs(a) >= 1) {
      result <- loga2(manifold, p2, yi)
    } else {
      p2 <- p2 * (a / abs(a)) # optimal alignment of p2 with p1
      b <- (1 - (abs(a)) ^ 2) ^ 0.5
      p2tilde <- (p2 - p1 * abs(a)) / b
      p1 <- as.vector(p1)
      p2tilde <- as.vector(p2tilde)
      result <- v - p1 %*% t(ip(manifold, v, p1)) - p2tilde %*% t(ip(manifold, v, p2tilde)) + p1 %*% t((abs(a)) * (ip(manifold, v, p1)) - b * (ip(manifold, v, p2tilde))) + p2tilde %*% t(b * (ip(manifold, v, p1)) + (abs(a)) * (ip(manifold, v, p2tilde)))
      result <- result * Conj(a / abs(a))
    }
  }
  return(result)
}

# Function for internal use in the geo_reg function. Similar to par_trans, but vectorized and without the errors, for speed and so as to not unexpectedly stop the algorithm because of one small error that can be safely ignored without changing the result. Difference from pt1 and pt2: p2 should be a column matrix, p1 and v should be matrices of the same dimensions with the same number of rows as p2.
pt3 <- function(manifold, p1, p2, v) {
  if (manifold == 'euclidean') {
    result <- v
  } else if (manifold == 'sphere') {
    p1 <- t(t(p1) / mag(manifold, p1))
    w <- loga3(manifold, p1, p2)
    e1 <- p1
    e2 <- t(t(w) / mag(manifold, w))
    a <- ip(manifold, v, e2)
    invar <- v - t(t(e2) * a)
    t <- mag(manifold, w)
    result <- t(t(e2) * (a * cos(t))) - t(t(e1) * (a * sin(t))) + invar
    index <- which(mag(manifold, w) == 0)
    result[, index] <- v[, index]
    #    p1 <- t(t(p1) / mag(manifold, p1))
    #    w <- loga3(manifold, p1, p2)
    #    t <- mag(manifold, w)
    #    result <- v - t(t(w + loga2(manifold, p2, p1)) * (ip(manifold, w, v) / (t ^ 2)))
    #    index <- which(t == 0)
    #    result[, index] <- v[, index]
  } else if (manifold == 'hyperbolic') {
    p1 <- t(t(p1) / sqrt(-ip(manifold, p1, p1)))
    w <- loga3(manifold, p1, p2)
    e1 <- p1
    e2 <- t(t(w) / mag(manifold, w))
    a <- ip(manifold, v, e2)
    invar <- v - t(t(e2) * a)
    t <- mag(manifold, w)
    result <- t(t(e2) * (a * cosh(t))) + t(t(e1) * (a * sinh(t))) + invar
    index <- which(mag(manifold, w) == 0)
    result[, index] <- v[, index]
    #    p1 <- t(t(p1) / sqrt(-ip(manifold, p1, p1)))
    #    w <- loga3(manifold, p1, p2)
    #    t <- mag(manifold, w)
    #    result <- v - t(t(w + loga2(manifold, p2, p1)) * (ip(manifold, w, v) / (t ^ 2)))
    #    index <- which(t == 0)
    #    result[, index] <- v[, index]
  } else if (manifold == 'kendall') {
    meanp1 <- colMeans(p1)
    meanp2 <- colMeans(p2)
    p1 <- t((t(p1) - meanp1) / mag(manifold, t(t(p1) - meanp1)))
    p2 <- t((t(p2) - meanp2) / mag(manifold, t(t(p2) - meanp2)))
    yi <- expo(manifold, p1, v)
    a <- ip(manifold, p1, as.vector(p2))
    p2 <- p2 %*% (a / abs(a)) # optimal alignment of p2 with p1
    b <- (1 - (abs(a)) ^ 2) ^ 0.5
    p2tilde <- t(t(p2 - t(t(p1) * abs(a))) / b)
    result <- v - t(t(p1) * (ip(manifold, v, p1))) - t(t(p2tilde) * (ip(manifold, v, p2tilde))) + t(t(p1) * ((abs(a)) * (ip(manifold, v, p1)) - b * (ip(manifold, v, p2tilde)))) + t(t(p2tilde) * (b * (ip(manifold, v, p1)) + (abs(a)) * (ip(manifold, v, p2tilde))))
    result <- t(t(result) * (Conj(a / abs(a))))
    index <- which(abs(a) >= 1)
    result[, index] <- loga(manifold, p2, yi)[, index]
  }
  return(result)
}

# Loss function for M-type estimators. t should be a vector of real numbers.
rho <- function(t, estimator, cutoff = NULL) {
  if (estimator == 'l2') {
    result <- 0.5 * t ^ 2
  } else if (estimator == 'l1') {
    result <- abs(t)
  } else if (estimator == 'huber') {
    result <- 0.5 * t ^ 2
    index <- which(abs(t) >= cutoff)
    result[index] <- cutoff * abs(t[index]) - 0.5 * cutoff ^ 2
  } else if (estimator == 'tukey') {
    result <- ((cutoff ^ 2) / 6) * (1 - (1 - (t / cutoff) ^ 2) ^ 3)
    result[which(abs(t) >= cutoff)] <- (cutoff ^ 2) / 6
  }
  return(result)
}

# Derivative of the loss function for M-type estimators. t should be a vector of real numbers.
rho_prime <- function(t, estimator, cutoff = NULL) {
  if (estimator == 'l2') {
    result <- t
  } else if (estimator == 'l1') {
    result <- sign(t)
  } else if (estimator == 'huber') {
    result <- t
    index <- which(abs(t) >= cutoff)
    result[index] <- cutoff * sign(t[index])
  } else if (estimator == 'tukey') {
    result <- t * ((1 - (t / cutoff) ^ 2) ^ 2) * sign(t)
    result[which(abs(t) >= cutoff)] <- 0
  }
  return(result)
}

# Calculates the residual vector for each data point. p, V, x, and y should all be matrices of appropriate dimensions.
eps <- function(manifold, p, V, x, y) {
  shifts <- V %*% t(x)
  predictions <- expo2(manifold, p, shifts)
  result <- loga(manifold, predictions, y)
  if (manifold == 'sphere') {
    result[, which((mag(manifold, y - t(t(predictions) * max(min(abs(ip(manifold, predictions, y)), 1), -1))) == 0) & (mag(manifold, predictions + y) < 0.000001))] <- 0
  } # we are ignoring cases where p2 is approximately -p1 in order to avoid an error
  return(result)
}

# Move tangent vectors v2 at expo(p, v1) to stangent vectors at p1 using Jacobi fields and adjoint operators; used in gradient calculations. p should be a column matrix, v1 and v2 should be matrices of the same dimension. The columns of v1 should be tangent to p, and the columns of v2 should be tangent to exp_map(p, 'corresponsing column in v1').
jacobi <- function(manifold, p, v1, v2) {
  result <- vector('list')
  L <- mag(manifold, v1)
  if (manifold == 'euclidean') {
    result$p <- v2
    result$V <- v2
  } else if (manifold == 'sphere') {
    v2_0 <- pt3(manifold, expo2(manifold, p, v1), p, v2)
    v2_tan <- t((t(v1) / L) * (ip(manifold, v2_0, t(t(v1) / L))))
    v2_orth <- v2_0 - v2_tan
    result$p <- t(t(v2_orth) * cos(L)) + v2_tan
    result$V <- t(t(v2_orth) * ((sin(L)) / L)) + v2_tan
    index <- which(L == 0)
    result$p[, index] <- v2[, index]
    result$V[, index] <- v2[, index]
  } else if (manifold == 'hyperbolic') {
    v2_0 <- pt3(manifold, expo2(manifold, p, v1), p, v2)
    v2_tan <- t((t(v1) / L) * (ip(manifold, v2_0, t(t(v1) / L))))
    v2_orth <- v2_0 - v2_tan
    result$p <- t(t(v2_orth) * cosh(L)) + v2_tan
    result$V <- t(t(v2_orth) * ((sinh(L)) / L)) + v2_tan
    index <- which(L == 0)
    result$p[, index] <- v2[, index]
    result$V[, index] <- v2[, index]
  } else if (manifold == 'kendall') {
    j <- (0 + 1i) * v1
    v2_0 <- pt3(manifold, expo2(manifold, p, v1), p, v2)
    w_0 <- t((t(j) / L) * (Re(ip(manifold, v2_0, t(t(j) / L)))))
    u_0 <- v2_0 - w_0
    w_tan <- t((t(v1) / L) * (Re(ip(manifold, w_0, t(t(v1) / L)))))
    w_orth <- w_0 - w_tan
    u_tan <- t((t(v1) / L) * (Re(ip(manifold, u_0, t(t(v1) / L)))))
    u_orth <- u_0 - u_tan
    result$p <- t(t(u_orth) * cos(L)) + t(t(w_orth) * cos(2 * L)) + u_tan + w_tan
    result$V <- t(t(u_orth) * ((sin(L)) / L)) + t(t(w_orth) * ((sin(2 * L)) / (2 * L))) + u_tan + w_tan
    index <- which(L == 0)
    result$p[, index] <- v2[, index]
    result$V[, index] <- v2[, index]
  }
  return(result)
}

# Calculates the gradient of the loss function at a given p, V. p, V, x, y and w should all be matrices or appropriate dimensions. resids is redundant: it is simply eps(manifold, p, V, x, y). However, including it as an argument to the function quickens the calculation.
grad <- function(manifold, p, V, x, y, w, resids, estimator, cutoff = NULL) {
  k <- dim(x)[2]
  result <- vector("list")
  mags <- mag(manifold, resids)
  shifts <- V %*% t(x)
  multiplier <- rho_prime(mags, estimator, cutoff)
  unit_resids<- t(t(resids) / mags)
  unit_resids[, which(mags == 0)] <- 0
  jf <- jacobi(manifold, p, shifts, unit_resids)
  result$p <- t(t(jf$p) * multiplier)
  result$p <- t(t(result$p) * w)
  result$V <- aperm(replicate(k, jf$V), c(1, 3, 2)) * aperm(replicate(dim(p)[1], w * x * multiplier), c(3, 2, 1))
  index <- which(mags <= 0.000001) # to avoid division by a small number
  result$p[, index] <- 0
  result$V[, , index] <- 0
  result$p <- as.matrix(-rowSums(result$p))
  result$V <- -rowSums(result$V, dims = 2)
  return(result)
}

#' Manifold check and projection
#'
#' Checks whether each data point in \eqn{y} is on the given manifold, and if
#' not, provides a modified version of \eqn{y} where each column has been
#' projected onto the manifold.
#'
#' @param manifold Type of manifold (\code{'euclidean'}, \code{'sphere'},
#'   \code{'hyperbolic'}, or \code{'kendall'}).
#' @param y A vector, matrix, or data frame whose columns should represent
#'   points on the manifold.
#' @return A named list containing \item{on}{a logical vector describing whether
#'   or not each column of \code{y} is on the manifold.} \item{data}{a matrix of
#'   data frame of the same dimensions as \code{y}; each column of \code{y} has
#'   been projected onto the manifold.}
#' @author Ha-Young Shin
#' @examples
#' y1 <- matrix(rnorm(10), ncol = 2)
#' y1 <- y1[, 1] + (1i) * y1[, 2]
#' y2 <- matrix(rnorm(10), ncol = 2)
#' y2 <- y2[, 1] + (1i) * y2[, 2]
#' y3 <- matrix(rnorm(10), ncol = 2)
#' y3 <- y3[, 1] + (1i) * y3[, 2]
#' y3 <- (y3 - mean(y3)) / norm(y3 - mean(y3), type = '2') # project onto preshape space
#' y <- matrix(c(y1, y2, y3), ncol = 3)
#' onmanifold('kendall', y)
#'
#' @export
onmanifold <- function(manifold, y) {
  y <- as.matrix(y)
  if (any(is.nan(y))) {
    stop('y should not contain NaN values')
  }
  sample_size <- dim(y)[2]
  result <- vector("list")
  if (manifold == 'euclidean') {
    result$on <- !logical(sample_size)
    result$data <- y
  } else if (manifold == 'sphere') {
    ons <- !logical(sample_size)
    mags <- mag(manifold, y)
    ons[which(abs(mags - 1) > 1e-6)] <- FALSE
    y <- t(t(y) / mags)
    result$on <- ons
    result$data <- y
  } else if (manifold == 'hyperbolic') {
    ons <- !logical(sample_size)
    mags <- sqrt(-ip(manifold, y, y))
    ons[which((abs(mags - 1) > 1e-6) | (y[1, ] < 0))] <- FALSE
    y <- t(t(y) / mags)
    result$on <- ons
    result$data <- y
  } else if (manifold == 'kendall') {
    ons <- !logical(sample_size)
    mags <- mag(manifold, y)
    means <- colMeans(y)
    ons[which((abs(mags - 1) > 1e-6) | (abs(means * Conj(means)) > 1e-6))] <- FALSE
    y <- t((t(y) - means) / mag(manifold, t(t(y) - means)))
    result$on <- ons
    result$data <- y
  } else {
    stop('the manifold must be one of euclidean, sphere, hyperbolic, or kendall')
  }
  return(result)
}

#' Gradient descent for (robust) geodesic regression
#'
#' Finds \eqn{\mathrm{argmin}_{(p,V)\in M\times (T_pM) ^ n}\sum_{i=1} ^ {N}
#' \rho(d(\mathrm{Exp}(p,Vx_i),y_i))} through a gradient descent algorithm.
#'
#' Each column of \code{x} should be centered to have an average of 0 for the
#' quickest and most accurate results. If all of the elements of a column of
#' \code{x} are equal, the resulting vector will consist of \code{NA}s. In the
#' case of the \code{'sphere'}, an error will be raised if all points are on a
#' pair of antipodes.
#'
#' @param manifold Type of manifold (\code{'euclidean'}, \code{'sphere'},
#'   \code{'hyperbolic'}, or \code{'kendall'}).
#' @param x A vector, matrix, or data frame of independent variables; for
#'   matrices and data frames, the rows and columns represent the subjects and
#'   independent variables, respectively.
#' @param y A matrix or data frame whose columns represent points on the
#'   manifold.
#' @param w A vector or matrix of weights
#' @param estimator M-type estimator (\code{'l2'}, \code{'l1'}, \code{'huber'},
#'   or \code{'tukey'}).
#' @param c Multiplier of \eqn{\sigma}, the square root of the variance, used in
#'   the cutoff parameter for the \code{'huber'} and \code{'tukey'} estimators;
#'   should be \code{NULL} for the \code{'l2'} or \code{'l1'} estimators.
#' @param p_tol Termination condition for the distance between consecutive
#'   updates of \code{p}.
#' @param V_tol Termination condition for the distance between columns of
#'   consecutive updates of \code{V}, parallel transported to be in the same
#'   tangent space. Can be a vector of positive real numbers for each
#'   independent variable or a single positive number.
#' @param max_iter Maximum number of gradient descent steps before ending the
#'   algorithm.
#' @return A named list containing \item{p}{a vector representing the estimate
#'   of the initial point on the manifold} \item{V}{a matrix representing the
#'   estimate of the initial velocities for each independent variable; the
#'   columns represent the independent variables.} \item{iteration}{number of
#'   gradient descent steps taken.}
#' @references Fletcher, P. T. (2013). Geodesic regression and the theory of
#'   least squares on Riemannian manifolds. International Journal of Computer
#'   Vision, 105, 171-185.
#'
#'   Kim, H. J., Adluru, N., Collins, M. D., Chung, M. K., Bendin, B. B.,
#'   Johnson, S. C., Davidson, R. J. and Singh, V. (2014). Multivariate general
#'   linear models (MGLM) on Riemannian manifolds with applications to
#'   statistical analysis of diffusion weighted images. 2014 IEEE Conference on
#'   Computer Vision and Pattern Recognition, 2705-2712.
#'
#'   Shin, H.-Y. and Oh H.-S. (2020). Robust Geodesic Regression.
#'   <arXiv:2007.04518>
#' @author Ha-Young Shin
#' @seealso \code{\link{intrinsic_location}}.
#' @examples
#' # an example of multiple regression with two independent variables, with 64
#' # data points
#'
#' x <- matrix(runif(2 * 64), ncol = 2)
#' x <- t(t(x) - colMeans(x))
#' y <- matrix(0L, nrow = 4, ncol = 64)
#' for (i in 1:64) {
#'   y[, i] <- exp_map('sphere', c(1, 0, 0, 0), c(0, runif(1), runif(1),
#'       runif(1)))
#' }
#' w <- rep(1/dim(y)[2],dim(y)[2])
#' geo_reg('sphere', x, y, w, 'tukey', c = are_nr('tukey', 2, 6))
#'
#' @export
#' @importFrom stats median
geo_reg <- function(manifold, x, y, w = rep(1/dim(y)[2],dim(y)[2]), estimator, c = NULL, p_tol = 1e-5, V_tol = 1e-5, max_iter = 100000) {
  if (((estimator == 'huber') | (estimator == 'tukey')) & is.null(c)) {
    stop('a c value must be provided if the chosen m-estimator is huber or tukey')
  }
  if (abs(sum(w)-1)>1e-9) {
    stop('weights must add to 1')
  }
  if (!is.null(c)) {
    if ((estimator == 'l2') | (estimator == 'l1')) {
      warning('l2 and l1 do not use a c value')
    }
    if (c <= 0) {
      stop('c must be positive')
    }
  }
  if (any(is.nan(x))) {
    stop('x should not contain NaN values')
  }
  ondata <- onmanifold(manifold, y)
  if (any(!(ondata$on))) {
    warning('all data points in y must lie on the manifold')
  }
  y <- ondata$data
  embedded <- dim(y)[1]
  sample_size <- dim(y)[2]
  x <- as.matrix(x)
  w <- as.vector(w)
  k <- dim(x)[2]
  allequal <- c()
  for (var in 1:k) { # Deals with case where all of the data for one of the independent variables are equal
    if (length(unique(x[, var])) == 1) {
      x[, var] <- numeric(sample_size)
      allequal <- c(allequal, var)
    }
  }
  if (manifold == 'sphere') {
    if (sample_size > 1) {
      antipode <- sum(mag(manifold, y[, 1] + y) == 0)
      same <- sum(mag(manifold, y[, 1] - y) == 0)
      if ((antipode + same == sample_size) & (antipode > 0)) {
        stop('there is no unique solution as all the data points are on antipodes of the sphere, with at least one data point on either antipode')
      }
    }
  }
  if (any(abs(colMeans(x)) > 0.000001)) {
    warning('the mean of the data for at least one of your independent variables is not zero; the x data should be centered for best/quickest results')
  }
  if  (length(allequal) == k) {
    current_p <- as.matrix(y[, 1])
  } else {
    current_p <- geo_reg(manifold, t(t(numeric(sample_size))), y, w, estimator, c, p_tol, V_tol, max_iter)$p
  }
  current_V <- matrix(0L, nrow = embedded, ncol = k)
  old_p <- current_p
  old_V <- current_V
  count <- 0
  alt_count <- 0
  cutoff <- NULL
  if ((estimator == 'huber') | (estimator == 'tukey')) {
    if (manifold == 'euclidean') {
      dimension <- embedded
    } else if ((manifold == 'sphere') | (manifold == 'hyperbolic')) {
      dimension <- embedded - 1
    } else if (manifold == 'kendall') {
      dimension <- 2 * embedded - 4
    }
    xi <- (2 * Pinv(dimension / 2, 0.5)) ^ 0.5
    current_shifts <- current_V %*% t(x)
    deviations <- dist(manifold, expo2(manifold, current_p, current_shifts), y)
    mad <- median(deviations)
    sigma <- mad / xi
    cutoff <- c * sigma
  }
  current_resids <- eps(manifold, current_p, current_V, x, y)
  current_loss <- sum(rho(mag(manifold, current_resids), estimator, cutoff) * w)
  step <- grad(manifold, current_p, current_V, x, y, w, current_resids, estimator, cutoff)
  V_diffs <- mag(manifold, (pt2(manifold, old_p, current_p, old_V) - current_V))
  lambda <- 0.1
  while (((count == 0) | ((count < max_iter) & ((dist(manifold, old_p, current_p) > p_tol) | (any(V_diffs > V_tol))))) & (alt_count < 100000)) {
    new_p <- tryCatch(expo(manifold, current_p, -lambda * step$p), warning = function(w) 'warning')
    while ((new_p[1] == 'warning') | (mag(manifold, -lambda * step$p) > 10)) {
      lambda <- lambda / 2
      new_p <- tryCatch(expo(manifold, current_p, -lambda * step$p), warning = function(w) 'warning')
    }
    new_V <- pt2(manifold, current_p, new_p, current_V - lambda * step$V)
    new_resids <- tryCatch(eps(manifold, new_p, new_V, x, y), warning = function(w) 'warning')
    while ((new_resids[1] == 'warning') | any(Re(ip(manifold, new_V, as.vector(new_p))) > 0.000001)) {
      lambda <- lambda / 2
      new_p <- expo(manifold, current_p, -lambda * step$p)
      new_V <- pt2(manifold, current_p, new_p, current_V - lambda * step$V)
      new_resids <- tryCatch(eps(manifold, new_p, new_V, x, y), warning = function(w) 'warning')
    }
    new_loss <- sum(rho(mag(manifold, new_resids), estimator, cutoff) * w)
    if (current_loss >= new_loss) {
      alt_count <- 0
      old_p <- current_p
      old_V <- current_V
      current_p <- new_p
      current_V <- new_V
      if ((estimator == 'huber') | (estimator == 'tukey')) {
        current_shifts <- current_V %*% t(x)
        deviations <- dist(manifold, expo2(manifold, current_p, current_shifts), y)
        mad <- median(deviations)
        sigma <- mad / xi
        cutoff <- c * sigma
      }
      current_resids <- new_resids
      current_loss <- sum(rho(mag(manifold, current_resids), estimator, cutoff) * w)
      step <- grad(manifold, current_p, current_V, x, y, w, current_resids, estimator, cutoff)
      V_diffs <- mag(manifold, (pt2(manifold, old_p, current_p, old_V) - current_V))
      if ((manifold == 'euclidean') | (manifold == 'sphere') | (manifold == 'kendall')) {
        lambda <- 8 * lambda
        #} else if (manifold == 'hyperbolic') {
        #  lambda <- 1.1 * lambda
      }
      count <- count + 1
    } else {
      lambda <- lambda / 2
      alt_count <- alt_count + 1
    }
  }
  if (count == max_iter) {
    warning('issues with convergence; make sure your independent variables are centered, and try adjusting max_iter, p_tol, or V_tol')
  }
  result <- vector("list")
  result$p <- current_p
  current_V[, allequal] <- NA
  result$V <- current_V
  result$iteration <- count
  return(result)
}

#' Gradient descent for location based on M-type estimators
#'
#' Finds \eqn{\mathrm{argmin}_{p\in M}\sum_{i=1} ^ {N} \rho(d(p,y_i))} through a
#' gradient descent algorithm.
#'
#' In the case of the \code{'sphere'}, an error will be raised if all points are
#' on a pair of antipodes.
#'
#' @param manifold Type of manifold (\code{'euclidean'}, \code{'sphere'},
#'   \code{'hyperbolic'}, or \code{'kendall'}).
#' @param y A matrix or data frame whose columns represent points on the
#'   manifold.
#' @param w A vector or matrix of weights.
#' @param estimator M-type estimator (\code{'l2'}, \code{'l1'}, \code{'huber'},
#'   or \code{'tukey'}).
#' @param c Multiplier of \eqn{\sigma}, the square root of the variance, used in
#'   the cutoff parameter for the \code{'huber'} and \code{'tukey'} estimators;
#'   should be \code{NULL} for the \code{'l2'} or \code{'l1'} estimators.
#' @param p_tol Termination condition for the distance between consecutive
#'   updates of \code{p}.
#' @param V_tol Termination condition for the distance between columns of
#'   consecutive updates of \code{V}, parallel transported to be in the same
#'   tangent space. Can be a vector of positive real numbers for each
#'   independent variable or a single positive number.
#' @param max_iter Maximum number of gradient descent steps before ending the
#'   algorithm.
#' @return A vector representing the location estimate
#' @references Fletcher, P. T. (2013). Geodesic regression and the theory of
#'   least squares on Riemannian manifolds. International Journal of Computer
#'   Vision, 105, 171-185.
#'
#'   Kim, H. J., Adluru, N., Collins, M. D., Chung, M. K., Bendin, B. B.,
#'   Johnson, S. C., Davidson, R. J. and Singh, V. (2014). Multivariate general
#'   linear models (MGLM) on Riemannian manifolds with applications to
#'   statistical analysis of diffusion weighted images. 2014 IEEE Conference on
#'   Computer Vision and Pattern Recognition, 2705-2712.
#'
#'   Shin, H.-Y. and Oh H.-S. (2020). Robust Geodesic Regression.
#'   <arXiv:2007.04518>
#' @author Ha-Young Shin
#' @seealso \code{\link{geo_reg}}, \code{\link[RiemBase]{rbase.mean}},
#'   \code{\link[RiemBase]{rbase.median}}.
#' @examples
#' y <- matrix(runif(100, 1000, 2000), nrow = 10)
#' intrinsic_location('euclidean', y, 'l2')
#'
#' @export
intrinsic_location <- function(manifold, y, w = rep(1/dim(y)[2],dim(y)[2]), estimator, c = NULL, p_tol = 1e-5, V_tol = 1e-5, max_iter = 100000) {
  sample_size <- dim(y)[2]
  return(geo_reg(manifold, numeric(sample_size), y, w, estimator, c, p_tol, V_tol, max_iter)$p)
}
