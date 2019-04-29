#' Compute posterior approximation given observed treatment assignments and a fixed covariance matrix
#'
#' @param y Set of observed treatment assignments (y in (0,1))
#' @param cov_matrix Covariance matrix; for examples, see \code{sqexp} or similar
#' @param tol Tolerance of algorithms.  Difference between the latent scores at each iteration - default 1e-2
#' @param max_iters Maximum number of iterations of the algorithm - default 20
#' @param verbose Decision to print progress to screen - default TRUE
#' @param approx_method Approximation method for posterior: 'ep' or 'laplace'
#' @param ep_vers 'Sequential' or 'Parallel' EP Algorithm - default \code{parallel}, alternative \code{sequential}
#' @return Object that contains the weights obtained from the balancing procedure and parameters from the optimization procedure
#'
#' The object that is returned is a list that contains the following entries
#' \itemize{
#' \item{ \code{Number_Iters} - Number of iterations for algorithm}
#' \item{ \code{PosteriorMean} - Posterior mean of latent scores}
#' \item{ \code{PosteriorVar} - Posterior covariance of latent scores}
#' \item{ \code{tilde_nu} - }
#' \item{ \code{tilde_tau} - }
#' \item{ \code{log_Z_ep} - EP Approximation to Log Likelihood}
#' \item{ \code{ComputationTime} - Runtime of EP algorithm for fixed covariance matrix}
#' \item{ \code{ps} - Probit transformed posterior mean}
#' }
#' @examples
#' n_obs <- 500
#' X1 <- rnorm(n_obs)
#' X2 <- rnorm(n_obs)
#' p <- pnorm( 0.5 * X1 + 0.5 * X2 )
#' TA <- rbinom(n_obs, 1, p)
#' dat <- data.frame(X1 = X1, X2 = X2, TA = TA)
#' covmat <- sqexp(cbind(X1, X2))
#' system.time(res <- gpbal_fixed(TA, covmat))
#' plot(res$ps, p, pch = 19, col = rgb(0,0,0,0.5))
#'
gpbal_fixed <- function(y,
                        cov_matrix,
                        tol=1e-2,
                        max_iters=20,
                        verbose = T,
                        approx_method = 'ep',
                        ep_vers = 'parallel'){

  ##############################################################################
  # Expectation Propagation Algorithm from Rasmussen & Williams
  # - Algorithms: 3.5, 3.6
  ##############################################################################
  start_time <- Sys.time()
  if(verbose){
    message('Starting: Expectation Propagation for Gaussian Process Classification')
    message('Start Time: ', start_time)
  }

  # Data checking.  Finding the classes and the number of observations
  classes = sort(unique(y))
  n_obs = length(y)

  # Data checking to ensure correct dimensions
  if(n_obs != nrow(cov_matrix)){
    message(paste('Error: Vec y of length:', n_obs,
                  'not equal to CovMat dimension:', nrow(cov_matrix), ncol(cov_matrix)))
  }

  if(length(classes) != 2){
    message("Error: Algorithm requires two classes exactly")
    return()
  }

  if(classes[1] == 0 && classes[2] == 1){
    y[y==0] = -1
    classes = sort(unique(y))
  } else if (classes[1] != -1 || classes[2] != 1){
    message("Error: Requires class labels -1 and 1 for Expectation Propagation")
    message(paste('Classes found: ', toString(classes)))
    return()
  }

  # Converting target inputs to a column vector
  y = matrix(y, nrow=n_obs, ncol=1)

  #-----------------------------------------------------------------------------
  # Running the Expectation Propagation Algorthim using C++ vvvvvvvvvvvvvvvvvvvv
  #-----------------------------------------------------------------------------
  if(tolower(approx_method) == 'ep'){
    if(length(grep(ep_vers, 'parallel'))){
      results <- par_ep(y, cov_matrix, tol, max_iters, verbose)
    } else if(length(grep(ep_vers, 'sequential'))) {
      results <- seq_ep(y, cov_matrix, tol, max_iters, verbose)
    } else {
      message("Error: ep_vers not compatible")
    }
  } else if (tolower(approx_method) == 'laplace') {
    results <- list()
    ests <- la_probit(y, cov_matrix, tol, max_iters)
    results[['PosteriorMean']] <- ests
  }

  #-----------------------------------------------------------------------------
  # Optimized C++ Code ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  #-----------------------------------------------------------------------------

  end_time <- Sys.time()
  dur_type <- 'mins'
  dur = difftime(end_time, start_time, units=dur_type)
  if(dur < 1){
    dur_type <- 'secs'
    dur = difftime(end_time, start_time, units=dur_type)
  }
  if(verbose){
    message('\nEnd Time: ', end_time)
    message('Duration: ', round(dur,3), paste(' ', dur_type, sep = ''))
    message('Approximate log Marginal Likelihood: ', round(results[['log_Z_ep']],3))
    message(paste('Number of Iterations:',results[['Number_Iters']]))
  }

  results[['ComputationTime']] = dur
  results[['ps']] <- pnorm(results$PosteriorMean)

  return(results)
}
