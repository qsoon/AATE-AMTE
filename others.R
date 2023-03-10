# Inner product. v1 and v2 should be vectors of the same length.
innerp <- function(manifold, v1, v2) {
  if ((manifold == 'euclidean') | (manifold == 'sphere')) {
    result <- sum(v1 * v2)
  } else if (manifold == 'hyperbolic') {
    v1[1] <- -v1[1]
    result <- sum(v1 * v2)
  } else if (manifold == 'kendall') {
    result <- sum(v1 * Conj(v2))
  }
  return(result)
}

# Magnitude of a vector
magnitude <- function(manifold, v) {
  return(Re(sqrt(innerp(manifold, v, v))))
}

#' Exponential map
#'
#' Performs the exponential map \eqn{\textrm{Exp}_p(v)} on the given manifold.
#'
#' @param manifold Type of manifold (\code{'euclidean'}, \code{'sphere'},
#'   \code{'hyperbolic'}, or \code{'kendall'}).
#' @param p A vector (or column matrix) representing a point on the manifold.
#' @param v A vector (or column matrix) tangent to \code{p}.
#' @return A vector representing a point on the manifold.
#' @references Fletcher, P. T. (2013). Geodesic regression and the theory of
#'   least squares on Riemannian manifolds. International Journal of Computer
#'   Vision, 105, 171-185.
#'
#'   Cornea, E., Zhu, H., Kim, P. and Ibrahim, J. G. (2017). Regression models
#'   on Riemannian symmetric spaces. Journal of the Royal Statistical Society:
#'   Series B, 79, 463-482.
#'
#'   Calinon, S. (2020). Gaussians on Riemannian manifolds: Applications for
#'   robot learning and adaptive control. IEEE Robotics & Automation Magazine,
#'   27, 33-45.
#'
#'   Shin, H.-Y. and Oh H.-S. (2020). Robust Geodesic Regression. <arXiv:2007.04518>
#' @author Ha-Young Shin
#' @seealso \code{\link{log_map}}.
#' @examples
#' exp_map('hyperbolic', c(1, 0, 0, 0, 0), c(0, 0, pi / 4, 0, 0))
#'
#' @export
exp_map <- function(manifold, p, v) {
  p <- as.vector(p)
  v <- as.vector(v)
  embedded <- length(p)
  if (embedded != length(v)) {
    stop('p and v must have the same length')
  }
  if (manifold == 'euclidean') {
    result <- p + v
  } else if (manifold == 'sphere') {
    magp <- magnitude(manifold, p)
    if (abs((magp - 1)) > 1e-6) {
      stop('p must be a unit vector')
    }
    if (abs(innerp(manifold, p, v)) > 1e-4) {
      stop('v must be tangent at p')
    }
    p <- p / magp # reprojects p onto the manifold, for precision
    theta <- magnitude(manifold, v)
    if (theta == 0) { # theta == 0 case must be dealt with separately due to division by theta
      result <- p
    } else {
      e1 <- p
      e2 <- v / theta
      result <- cos(theta) * e1 + sin(theta) * e2
      result <- result / magnitude(manifold, result) # reprojects result onto the manifold, for precision
    }
  } else if (manifold == 'hyperbolic') {
    magp <- sqrt(-innerp(manifold, p, p))
    if ((abs((magp - 1)) > 1e-6) | (p[1] < 0)) {
      stop('p must lie on the upper sheet of the unit hyperboloid')
    }
    if (abs(innerp(manifold, p, v)) > 1e-5) {
      stop('v must be tangent at p')
    }
    p <- p / magp # reprojects p onto the manifold, for precision
    theta <- magnitude(manifold, v)
    if (theta == 0) { # theta == 0 case must be dealt with separately due to division by theta
      result <- p
    } else {
      e1 <- p
      e2 <- v / theta
      result <- cosh(theta) * e1 + sinh(theta) * e2
      result <- result / sqrt(-innerp(manifold, result, result)) # reprojects result onto the manifold, for precision
    }
  } else if (manifold == 'kendall') {
    meanp <- sum(p) / embedded
    if (abs((magnitude(manifold, p) - 1)) > 1e-6) {
      stop('p must be a unit vector')
    }
    if (abs(innerp(manifold, p, v)) > 1e-6) {
      stop('v must be tangent at p')
    }
    if ((magnitude(manifold, meanp) > 1e-6) | (magnitude(manifold, sum(v) / embedded) > 1e-6)) {
      stop('p and v must be centered')
    }
    p <- (p - meanp) / magnitude(manifold, p - meanp) # reprojects p onto the manifold, for precision
    theta <- magnitude(manifold, v)
    if (theta == 0) { # theta == 0 case must be dealt with separately due to division by theta
      result <- p
    } else {
      e1 <- p
      e2 <- v / theta
      result <- cos(theta) * e1 + sin(theta) * e2
      result <- (result - sum(result) / embedded) / magnitude(manifold, result - sum(result) / embedded) # reprojects result onto the manifold, for precision
    }
  } else {
    stop('the manifold must be one of euclidean, sphere, hyperbolic, or kendall')
  }
  return(result)
}

#' Logarithm map
#'
#' Performs the logarithm map \eqn{\textrm{Log}_{p_1}(p_2)} on the given
#' manifold, provided \eqn{p_2} is in the domain of \eqn{\textrm{Log}_{p_1}}.
#'
#' On the sphere, \eqn{-p_1} is not in the domain of \eqn{\textrm{Log}_{p_1}}.
#'
#' @param manifold Type of manifold (\code{'euclidean'}, \code{'sphere'},
#'   \code{'hyperbolic'}, or \code{'kendall'}).
#' @param p1 A vector (or column matrix) representing a point on the manifold.
#' @param p2 A vector (or column matrix) representing a point on the manifold.
#' @return A vector tangent to \code{p1}.
#' @references Fletcher, P. T. (2013). Geodesic regression and the theory of
#'   least squares on Riemannian manifolds. International Journal of Computer
#'   Vision, 105, 171-185.
#'
#'   Cornea, E., Zhu, H., Kim, P. and Ibrahim, J. G. (2017). Regression models
#'   on Riemannian symmetric spaces. Journal of the Royal Statistical Society:
#'   Series B, 79, 463-482.
#'
#'   Calinon, S. (2020). Gaussians on Riemannian manifolds: Applications for
#'   robot learning and adaptive control. IEEE Robotics & Automation Magazine,
#'   27, 33-45.
#'
#'   Shin, H.-Y. and Oh H.-S. (2020). Robust Geodesic Regression.
#'   <arXiv:2007.04518>
#' @author Ha-Young Shin
#' @seealso \code{\link{exp_map}}, \code{\link{geo_dist}}.
#' @examples
#' log_map('sphere', c(0, 1, 0, 0), c(0, 0, 1, 0))
#'
#' @export
log_map <- function(manifold, p1, p2) {
  p1 <- as.vector(p1)
  p2 <- as.vector(p2)
  embedded <- length(p1)
  if (embedded != length(p2)) {
    stop('p1 and p2 must have the same length')
  }
  if (manifold == 'euclidean') {
    result <- p2 - p1
  } else if (manifold == 'sphere') {
    magp1 <- magnitude(manifold, p1)
    magp2 <- magnitude(manifold, p2)
    if ((abs((magp1 - 1)) > 1e-6) | (abs((magp2 - 1)) > 1e-6)) {
      stop('p1 and p2 must be unit vectors')
    }
    p1 <- p1 / magp1 # reprojects p1 onto the manifold, for precision
    p2 <- p2 / magp2 # reprojects p2 onto the manifold, for precision
    a <- max(min(innerp(manifold, p1, p2), 1), -1) # ensures a is in [-1, 1]
    theta <- acos(a)
    tang <- p2 - a * p1
    t <- magnitude(manifold, tang)
    if (t == 0) { # t == 0 case must be dealt with separately due to division by t
      if (magnitude(manifold, p1 - p2) < 1e-6) { # determining whether t == 0 because of p1 = p2 or p1 = -p2
        result <- numeric(embedded)
      } else if (magnitude(manifold, p1 + p2) < 1e-6) { # determining whether t == 0 because of p1 = p2 or p1 = -p2
        stop('p2 is the antipode of p1 and is therefore not in the domain of the log map at p1')
      }
    } else {
      result <- theta * (tang / t)
    }
  } else if (manifold == 'hyperbolic') {
    magp1 <- sqrt(-innerp(manifold, p1, p1))
    magp2 <- sqrt(-innerp(manifold, p2, p2))
    if ((abs((magp1 - 1)) > 1e-6) | (abs((magp2 - 1)) > 1e-6) | (p1[1] < 0) | (p2[1] < 0)) {
      stop('p1 and p2 must lie on the upper sheet of the unit hyperboloid')
    }
    p1 <- p1 / magp1 # reprojects p1 onto the manifold, for precision
    p2 <- p2 / magp2 # reprojects p2 onto the manifold, for precision
    a <- min(innerp(manifold, p1, p2), -1) # ensures -a is at least 1
    theta <- acosh(-a)
    tang <- p2 + a * p1
    t <- magnitude(manifold, tang)
    if (t == 0) { # t == 0 case must be dealt with separately due to division by t
      result <- numeric(embedded)
    } else {
      result <- theta * (tang / t)
    }
  } else if (manifold == 'kendall') {
    meanp1 <- sum(p1) / embedded
    meanp2 <- sum(p2) / embedded
    if ((abs((magnitude(manifold, p1) - 1)) > 1e-6) | (abs((magnitude(manifold, p2) - 1)) > 1e-6)) {
      stop('p1 and p2 must be unit vectors')
    }
    if ((magnitude(manifold, meanp1) > 1e-6) | (magnitude(manifold, meanp2) > 1e-6)) {
      stop('p1 and p2 must be centered')
    }
    p1 <- (p1 - meanp1) / magnitude(manifold, p1 - meanp1) # reprojects p1 onto the manifold, for precision
    p2 <-  (p2 - meanp2) / magnitude(manifold, p2 - meanp2) # reprojects p2 onto the manifold, for precision
    a <- innerp(manifold, p1, p2)
    theta <- acos(max(min(abs(a), 1), -1))  # ensures argument is in [-1, 1]
    tang <- (a / abs(a)) * p2 - abs(a) * p1
    t <- magnitude(manifold, tang)
    if (t == 0) { # t == 0 case must be dealt with separately due to division by t
      result <- numeric(embedded)
    } else {
      result <- theta * (tang / t)
    }
  } else {
    stop('the manifold must be one of euclidean, sphere, hyperbolic, or kendall')
  }
  return(result)
}

#' Geodesic distance between two points on a manifold
#'
#' Finds the Riemannian distance
#' \eqn{d(p_1,p_2)=||\textrm{Log}_{p_1}(p_2)||} between two points on
#' the given manifold, provided \eqn{p_2} is in the domain of
#' \eqn{\textrm{Log}_{p_1}}.
#'
#' On the sphere, \eqn{-p_1} is not in the domain of \eqn{\textrm{Log}_{p_1}}.
#'
#' @inheritParams log_map
#' @return Riemannian distance between \code{p1} and \code{p2}.
#' @author Ha-Young Shin
#' @seealso \code{\link{log_map}}.
#' @examples
#' p1 <- matrix(rnorm(10), ncol = 2)
#' p1 <- p1[, 1] + (1i) * p1[, 2]
#' p1 <- (p1 - mean(p1)) / norm(p1 - mean(p1), type = '2')
#' p2 <- matrix(rnorm(10), ncol = 2)
#' p2 <- p2[, 1] + (1i) * p2[, 2]
#' p2 <- (p2 - mean(p2)) / norm(p2 - mean(p2), type = '2')
#' geo_dist('kendall', p1, p2)
#'
#' @export
geo_dist <- function(manifold, p1, p2) {
  return(magnitude(manifold, log_map(manifold, p1, p2)))
}

#' Parallel transport
#'
#' Performs \eqn{\Gamma_{p_1 \rightarrow p_2}(v)}, parallel transport along the
#' unique minimizing geodesic connecting \eqn{p_1} and \eqn{p_2}, if it exists,
#' on the given manifold.
#'
#' On the sphere, there is no unique minimizing geodesic connecting \eqn{p_1}
#' and \eqn{-p_1}.
#'
#' @param manifold Type of manifold (\code{'euclidean'}, \code{'sphere'},
#'   \code{'hyperbolic'}, or \code{'kendall'}).
#' @param p1 A vector (or column matrix) representing a point on the manifold.
#' @param p2 A vector (or column matrix) representing a point on the manifold.
#' @param v A vector (or column matrix) tangent to \code{p1}.
#' @return A vector tangent to \code{p2}.
#' @references Fletcher, P. T. (2013). Geodesic regression and the theory of
#'   least squares on Riemannian manifolds. International Journal of Computer
#'   Vision, 105, 171-185.
#'
#'   Cornea, E., Zhu, H., Kim, P. and Ibrahim, J. G. (2017). Regression models
#'   on Riemannian symmetric spaces. Journal of the Royal Statistical Society:
#'   Series B, 79, 463-482.
#'
#'   Calinon, S. (2020). Gaussians on Riemannian manifolds: Applications for
#'   robot learning and adaptive control. IEEE Robotics & Automation Magazine,
#'   27, 33-45.
#'
#'   Shin, H.-Y. and Oh H.-S. (2020). Robust Geodesic Regression. <arXiv:2007.04518>
#' @author Ha-Young Shin
#' @examples
#' p1 <- matrix(rnorm(10), ncol = 2)
#' p1 <- p1[, 1] + (1i) * p1[, 2]
#' p1 <- (p1 - mean(p1)) / norm(p1 - mean(p1), type = '2') # project onto pre-shape space
#' p2 <- matrix(rnorm(10), ncol = 2)
#' p2 <- p2[, 1] + (1i) * p2[, 2]
#' p2 <- (p2 - mean(p2)) / norm(p2 - mean(p2), type = '2') # project onto pre-shape space
#' p3 <- matrix(rnorm(10), ncol = 2)
#' p3 <- p3[, 1] + (1i) * p3[, 2]
#' p3 <- (p3 - mean(p3)) / norm(p3 - mean(p3), type = '2') # project onto pre-shape space
#' v <- log_map('kendall', p1, p3)
#' par_trans('kendall', p1, p2, v)
#'
#' @export
par_trans <- function(manifold, p1, p2, v) {
  p1 <- as.vector(p1)
  p2 <- as.vector(p2)
  v <- as.vector(v)
  embedded <- length(p1)
  if ((embedded != length(p2)) | (embedded != length(v))) {
    stop('p1, p2, and v must have the same length')
  }
  if (manifold == 'euclidean') {
    result <- v
  } else if (manifold == 'sphere') {
    if (abs(innerp(manifold, p1, v)) > 1e-4) {
      stop('v must be tangent at p1')
    }
    p1 <- p1 / magnitude(manifold, p1) # reprojects p1 onto the manifold, for precision
    w <- log_map(manifold, p1, p2)
    t <- magnitude(manifold, w)
    if (t == 0) { # t == 0 case must be dealt with separately due to division by t
      result <- v
    } else {
      e1 <- p1
      e2 <- w / t
      a <- innerp(manifold, v, e2)
      invar <- v - a * e2
      result <- a * (cos(t) * e2 - sin(t) * e1) + invar
    }
    #p1 <- p1 / magnitude(manifold, p1) # reprojects p1 onto the manifold, for precision
    #w <- log_map(manifold, p1, p2)
    #t <- magnitude(manifold, w)
    #if (t == 0) { # t == 0 case must be dealt with separately due to division by t
    #  result <- v
    #} else {
    #  result <- v - (innerp(manifold, w, v) / (t ^ 2)) * (w + log_map(manifold, p2, p1))
    #}
  } else if (manifold == 'hyperbolic') {
    if (abs(innerp(manifold, p1, v)) > 1e-5) {
      stop('v must be tangent at p1')
    }
    p1 <- p1 / sqrt(-innerp(manifold, p1, p1)) # reprojects p1 onto the manifold, for precision
    w <- log_map(manifold, p1, p2)
    t <- magnitude(manifold, w)
    if (t == 0) { # t == 0 case must be dealt with separately due to division by t
      result <- v
    } else {
      e1 <- p1
      e2 <- w / t
      a <- innerp(manifold, v, e2)
      invar <- v - a * e2
      result <- a * (cosh(t) * e2 + sinh(t) * e1) + invar
    }
    #p1 <- p1 / sqrt(-innerp(manifold, p1, p1)) # reprojects p1 onto the manifold, for precision
    #w <- log_map(manifold, p1, p2)
    #t <- magnitude(manifold, w)
    #if (t == 0) { # t == 0 case must be dealt with separately due to division by t
    #  result <- v
    #} else {
    #  result <- v - (innerp(manifold, w, v) / (t ^ 2)) * (w + log_map(manifold, p2, p1))
    #}
  } else if (manifold == 'kendall') {
    meanp1 <- sum(p1) / embedded
    meanp2 <- sum(p2) / embedded
    if ((abs((magnitude(manifold, p2) - 1)) > 1e-6) | (magnitude(manifold, meanp2) > 1e-6)) {
      stop('p2 must be a centered unit vector')
    }
    p1 <- (p1 - meanp1) / magnitude(manifold, p1 - meanp1) # reprojects p1 onto the manifold, for precision
    p2 <- (p2 - meanp2) / magnitude(manifold, p2 - meanp2) # reprojects p2 onto the manifold, for precision
    yi <- exp_map(manifold, p1, v)
    a <- innerp(manifold, p1, p2)
    p2 <- (a / abs(a)) * p2 # optimal alignment of p2 with p1
    if (abs(a) >= 1) { # i.e. p1 == p2
      result <- log_map(manifold, p2, yi)
    } else {
      b <- sqrt(1 - (abs(a)) ^ 2)
      p2tilda <- (p2 - abs(a) * p1) / b
      result <- v - (innerp(manifold, v, p1)) * p1 - (innerp(manifold, v, p2tilda)) * p2tilda + ((abs(a)) * (innerp(manifold, v, p1)) - b * (innerp(manifold, v, p2tilda))) * p1 + (b * (innerp(manifold, v, p1)) + (abs(a)) * (innerp(manifold, v, p2tilda))) * p2tilda
      result <- (Conj(a / abs(a))) * result
    }
  } else {
    stop('the manifold must be one of euclidean, sphere, hyperbolic, or kendall')
  }
  return(result)
}

#' Loss
#'
#' Loss for a given \code{p} and \code{V}.
#'
#' @param manifold Type of manifold (\code{'euclidean'}, \code{'sphere'},
#'   \code{'hyperbolic'}, or \code{'kendall'}).
#' @param p A vector (or column matrix) on the manifold.
#' @param V A matrix (or vector) where each column is a vector in the tangent
#'   space at \code{p}.
#' @param x A matrix or data frame of independent variables; for matrices and
#'   data frames, the rows and columns represent the subjects and independent
#'   variables, respectively.
#' @param y A matrix or data frame (or vector) whose columns represent points on
#'   the manifold.
#' @param estimator M-type estimator (\code{'l2'}, \code{'l1'}, \code{'huber'},
#'   or \code{'tukey'}).
#' @param cutoff Cutoff parameter for the \code{'huber'} and \code{'tukey'}
#'   estimators; should be \code{NULL} for the \code{'l2'} or \code{'l1'}
#'   estimators.
#' @return Loss.
#' @author Ha-Young Shin
#' @examples
#' y <- matrix(0L, nrow = 3, ncol = 64)
#' for (i in 1:64) {
#'   y[, i] <- exp_map('hyperbolic', c(1, 0, 0), c(0, runif(1), runif(1)))
#' }
#' intrinsic_mean <- intrinsic_location('hyperbolic', y, 'l2')
#' loss('hyperbolic', intrinsic_mean, numeric(3), numeric(64), y, 'l2')
#'
#' @export
loss <- function(manifold, p, V, x, y, estimator, cutoff = NULL) {
  p <- as.matrix(p)
  V <- as.matrix(V)
  x <- as.matrix(x)
  y <- as.matrix(y)
  sample_size <- dim(y)[2]
  k <- dim(x)[2]
  if ((dim(y)[1] != dim(p)[1]) | (dim(V)[1] != dim(p)[1])) {
    stop('p, each vector in V, and each data point in y must have the same length')
  }
  if (dim(x)[1] != sample_size) {
    stop('the sample sizes according to x and y do not match')
  }
  if (manifold == 'sphere') {
    if (abs((magnitude(manifold, p) - 1)) > 1e-6) {
      stop('p must be a unit vector')
    }
    for (h in 1:k) {
      if (abs(innerp(manifold, p, V[, h])) > 1e-4) {
        stop('each column of V must be tangent at p')
      }
    }
  } else if (manifold == 'hyperbolic') {
    if (abs((sqrt(-innerp(manifold, p, p)) - 1)) > 1e-6) {
      stop('p must lie on the upper sheet of the unit hyperboloid')
    }
    for (h in 1:k) {
      if (abs(innerp(manifold, p, V[, h])) > 1e-5) {
        stop('each column of V must be tangent at p')
      }
    }
  } else if (manifold == 'kendall') {
    if (abs((magnitude(manifold, p) - 1)) > 1e-6) {
      stop('p must be a unit vector')
    }
    if (magnitude(manifold, mean(p)) > 1e-6) {
      stop('p must be centered')
    }
    for (l in 1:k) {
      if (abs(innerp(manifold, p, V[, l])) > 1e-6) {
        stop('each column of V must be tangent at p')
      }
      if (magnitude(manifold, sum(V[, l]) / dim(p)[1]) > 1e-6) {
        stop('each column of V must be centered')
      }
    }
  }
  shifts <- V %*% t(x)
  predictions <- expo2(manifold, p, shifts)
  resids <- loga(manifold, predictions, y)
  result <- sum(rho(magnitude(manifold, resids), estimator, cutoff))
  return(result)
}
