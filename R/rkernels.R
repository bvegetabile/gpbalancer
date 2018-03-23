################################################################################
#
# Gaussian Process Approximation via Expectation Propagation
# - Speed gains provided by Rcpp and Rcpp Armadillo
#
# Author: Brian Vegetabile
#
# References: Gaussian Processes for Machine Learning
#             - Carl Rasmussen, Christopher Williams
################################################################################

#-------------------------------------------------------------------------------
# Kernel Functions -------------------------------------------------------------
#-------------------------------------------------------------------------------

# Squared Exponential-----------------------------------------------------------

sqexp <- function(X, theta = rep(1, ncol(X)+1), sig_noise=0.000001){
  # Main input is a variable called 'theta'.  The organization of the parameters
  # of of the theta vector is as follows
  # --- theta[1] - Outer Variance Term
  # --- theta[2:(num_params+1)] - Length scales for each dimension of X.

  # Finding the number of dimensions and the number of observations
  if (class(X) == "matrix" | class(X) == "data.frame") {
    n_obs <- dim(X)[1]
    num_params <- dim(X)[2]
  } else {
    n_obs <- length(X)
    num_params <- 1
    X <- matrix(X, nrow=n_obs, ncol=num_params)
  }

  sig = theta[1]
  ls = theta[2:(num_params+1)]

  if(length(ls)==num_params){
    M_inv = diag(c(ls^2), num_params, num_params)
  } else{
    M_inv = diag(rep(ls^2, num_params))
  }

  cov_mat <- construct_sqexp(X, M_inv, sig, sig_noise)
  return(cov_mat)
}

sqexp_par <- function(X, theta = rep(1, ncol(X)+1), sig_noise=0.000001){
  # Main input is a variable called 'theta'.  The organization of the parameters
  # of of the theta vector is as follows
  # --- theta[1] - Outer Variance Term
  # --- theta[2:(num_params+1)] - Length scales for each dimension of X.

  # Finding the number of dimensions and the number of observations
  if (class(X) == "matrix" | class(X) == "data.frame") {
    n_obs <- dim(X)[1]
    num_params <- dim(X)[2]
  } else {
    n_obs <- length(X)
    num_params <- 1
    X <- matrix(X, nrow=n_obs, ncol=num_params)
  }

  cov_mat <- par_sqexp(X, theta, sig_noise)
  return(cov_mat)
}

sqexp_poly <- function(X, theta, noise = 1e-4){
  scale0 <- 1
  sig0 <- theta[1]
  scale1 <- 1
  ls1 <- theta[2]

  K1 <- polykernel(X, sig0, pwr = 1,
                   scale = scale0, noise=noise)
  K2 <- sqexp_common(X, lengthscale = ls1,
                     scale = scale1, noise = noise)

  return(K1 + K2)
}

sqexp_poly2 <- function(X, theta, noise = 1e-4){
  scale0 <- 1
  sig0 <- theta[1]
  scale1 <- 1
  ls1 <- theta[2]

  K1 <- polykernel(X, sig0, pwr = 2,
                   scale = scale0, noise=noise)
  K2 <- sqexp_common(X, lengthscale = ls1,
                     scale = scale1, noise = noise)

  return(K1 + K2)
}


