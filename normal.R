# used in the Newton-Raphson process in the inflection_point function
xtangentx <- function(manifold, x) {
  if (manifold == 'sphere') {
    return(x * tan(x))
  } else if (manifold == 'hyperbolic') {
    return(x * tanh(x))
  }
}

# used in the Newton-Raphson process in the inflection_point function
deriv_xtangentx <- function(manifold, x) {
  if (manifold == 'sphere') {
    return (tan(x) + x / (cos(x) ^ 2))
  } else if (manifold == 'hyperbolic') {
    return (tanh(x) + x / (cosh(x) ^ 2))
  }
}

# inflection point of H for dimension m given sigma_sq which will be in
# inv_func, found using the Newton-Raphson method
inflection_point <- function(manifold, m, sigma_sq) {
  if (manifold == 'sphere') {
    new_x <- pi / 2 - 0.1
    while (xtangentx(manifold, new_x) < m * sigma_sq) {
      new_x <- pi / 2 - (pi / 2 - new_x) / 10
    }
  } else if (manifold == 'hyperbolic') {
    new_x <- 1.199678640257733833916369849 # inflection point of x tanh(x), is equal to the point where x tanh(x) = 1
  }
  old_x <- new_x + 5
  count <- 0
  while (abs(old_x - new_x) > 0.000000001) {
    old_x <- new_x
    new_x <- old_x - (xtangentx(manifold, old_x) - m * sigma_sq) / deriv_xtangentx(manifold, old_x)
    count <- count + 1
    if (count > 1000) {
      return ('fail')
    }
  }
  return (new_x)
}

# the error function from the pracma package gives an error when z = 0; this
# erfz has been modified to remove the error
erfz <- function(z)
{
  if (is.null(z))
    return(NULL)
  else if (!is.numeric(z) && !is.complex(z))
    stop("Argument 'z' must be a numeric or complex scalar or vector.")
  a0 <- abs(z)
  c0 <- exp(-z * z)
  z1 <- ifelse(Re(z) < 0, -z, z)
  i <- a0 <= 5.8
  work.i <- i
  cer <- rep(NA, length = length(z))
  if (sum(work.i) > 0) {
    cs <- z1
    cr <- cs
    for (k in 1:120) {
      cr[work.i] <- cr[work.i] * z1[work.i] * z1[work.i]/(k +
                                                            0.5)
      cs[work.i] <- cs[work.i] + cr[work.i]
    }
    cer[i] <- 2 * c0[i] * cs[i]/sqrt(pi)
  }
  work.i <- !i
  if (sum(work.i) > 0) {
    cl <- 1/z1
    cr <- cl
    for (k in 1:13) {
      cr[work.i] <- -cr[work.i] * (k - 0.5)/(z1[work.i] *
                                               z1[work.i])
      cl[work.i] <- cl[work.i] + cr[work.i]
    }
    cer[!i] <- 1 - c0[!i] * cl[!i]/sqrt(pi)
  }
  cer[Re(z) < 0] <- -cer[Re(z) < 0]
  return(cer)
}

# the distribution function (after removing a constant scale factor and shift)
# of r = d(y, mu) for dimension m given sigma_sq
distR <- function(manifold, m, sigma_sq, r) {
  sum <- 0
  if (manifold == 'sphere') {
    for (j in 0:m) {
      sum <- sum + (factorial(m) / (factorial(j) * factorial(m - j))) * ((-1) ^ j) * exp((-sigma_sq * (m - 2 * j) ^ 2) / 2) * erfz(r / ((2 * sigma_sq) ^ 0.5) + ((sigma_sq / 2) ^ 0.5) * (m - 2 * j) * 1i)
    }
    sum <- ((1i) ^ m) * sum
  } else if (manifold == 'hyperbolic') {
    for (j in 0:m) {
      sum <- sum + (factorial(m) / (factorial(j) * factorial(m - j))) * ((-1) ^ j) * exp((sigma_sq * (m - 2 * j) ^ 2) / 2) * erfz(r / ((2 * sigma_sq) ^ 0.5) - ((sigma_sq / 2) ^ 0.5) * (m - 2 * j))
    }
  }
  return(Re(sum))
}

# limit of distR as r goes to infinity in the hyperbolic case
lim_distR <- function(manifold, m, sigma_sq) {
  sum <- 0
  for (j in 0:m) {
    sum <- sum + (factorial(m) / (factorial(j) * factorial(m - j))) * ((-1) ^ j) * exp((sigma_sq * (m - 2 * j) ^ 2) / 2)
  }
  return(sum)
}

# derivative of distR
deriv_distR <- function(manifold, m, sigma_sq, r) {
  if (manifold == 'sphere') {
    return (exp((-r ^ 2) / (2 * sigma_sq)) * (sin(r) ^ m) * (2 ^ m) * ((2 / (pi * sigma_sq)) ^ 0.5))
  } else if (manifold == 'hyperbolic') {
    return (exp((-r ^ 2) / (2 * sigma_sq)) * (sinh(r) ^ m) * (2 ^ m) * ((2 / (pi * sigma_sq)) ^ 0.5))
  }
}

# inverse of distR at t, using the Newton-Raphson method with distR and
# deriv_distR
inv_func_distR <- function(manifold, m, sigma_sq, t) {
  new_x <- inflection_point(manifold, m, sigma_sq) # use the inflection point of H as the starting point for the Newton-Raphson method
  old_x <- new_x + 5
  count <- 0
  while (any(abs(old_x - new_x) > 0.0000000001)) {
    old_x <- new_x
    new_x <- old_x - (distR(manifold, m, sigma_sq, old_x) - t) / deriv_distR(manifold, m, sigma_sq, old_x)
    count <- count + 1
    if (count > 1000) {
      return ('fail')
    }
  }
  return (new_x)
}

#' Random generation of tangent vectors from the Riemannian normal distribution
#'
#' Random generation of tangent vectors from the Riemannian normal distribution
#' on the \code{n}-dimensional sphere or hyperbolic space at mean \code{(1, 0,
#' ..., 0)}, a vector of length \code{n+1}.
#'
#' Tangent vectors are of the form \eqn{\mathrm{Log}(\mu, y)} in the tangent
#' space at the Fr\'echet mean \eqn{\mu} = \code{(1, 0, ..., 0)}, which is
#' isomorphic to \code{n}-dimensional Euclidean space, where \eqn{y} has a
#' Riemannian normal distribution. The first element of these vectors
#' will always be 0 at this \eqn{\mu}. These vectors can be
#' transported to any other \eqn{\mu} on the manifold.
#'
#' @param manifold Type of manifold (\code{'sphere'} or \code{'hyperbolic'}).
#' @param N Number of points to generate.
#' @param n Dimension of the manifold.
#' @param sigma_sq A scale parameter.
#' @return An \code{(n+1)}-by-\code{N} matrix where each column represents a random
#'   tangent vector at \code{(1, 0, ..., 0)}.
#' @references Fletcher, P. T. (2013). Geodesic regression and the theory of
#'   least squares on Riemannian manifolds. International Journal of Computer
#'   Vision, 105, 171-185.
#'
#'   Fletcher, T. (2020). Statistics on manifolds. In \emph{Riemannian Geometric
#'   Statistics in Medical Image Analysis}. 39--74. Academic Press.
#'
#'   Shin, H.-Y. and Oh H.-S. (2020). Robust Geodesic Regression.
#'   <arXiv:2007.04518>
#' @author Ha-Young Shin
#' @examples
#'
#' sims <- rnormtangents('hyperbolic', N = 4, n = 2, sigma_sq = 1)
#'
#' @export
#' @importFrom stats runif
rnormtangents <- function(manifold, N, n, sigma_sq) {
  if (manifold == 'sphere') {
    sup_distR <- distR(manifold, n - 1, sigma_sq, pi)
  } else if (manifold == 'hyperbolic') {
    sup_distR <- lim_distR(manifold, n - 1, sigma_sq)
  } else {
    stop('the manifold must be either sphere or hyperbolic')
  }
  u <- runif(N, 0, 1)
  t <- u * (sup_distR - distR(manifold, n - 1, sigma_sq, 0)) + distR(manifold, n - 1, sigma_sq, 0)
  R <- inv_func_distR(manifold, n - 1, sigma_sq, t)
  while (R[1] == 'fail') {
    u <- runif(N, 0, 1)
    t <- u * (sup_distR - distR(manifold, n - 1, sigma_sq, 0)) + distR(manifold, n - 1, sigma_sq, 0)
    R <- inv_func_distR(manifold, n - 1, sigma_sq, t)
  }
  direction <- cbind(numeric(N), matrix(MASS::mvrnorm(n = N, mu = integer(n), Sigma = diag(x = 10, nrow = n)), nrow = N))
  direction <- direction / (rowSums(direction ^ 2) ^ 0.5)
  return (t(R * direction))
}
