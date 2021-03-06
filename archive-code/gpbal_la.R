#' Compute posterior approximation given observed treatment assignments and a fixed covariance matrix
#'
#' @param y Set of observed treatment assignments (y \in (0,1))
#' @param cov_matrix Covariance matrix; for examples, see \code{sqexp} or similar
#' @param tol Tolerance of EP Algorithm.  Difference between the latent scores at each iteration - default 1e-2
#' @param max_iters Maximum number of iterations of the EP Algorithm - default 20
#' @param verbose Decision to print progress to screen - default TRUE
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
#' system.time(res <- gpbal_la_fixed(TA, covmat))
#' plot(res$ps, p, pch = 19, col = rgb(0,0,0,0.5))
gpbal_la_fixed <- function(y,
                           cov_matrix,
                           tol=1e-2,
                           max_iters=20,
                           verbose = T){
  ##############################################################################
  # Expectation Propagation Algorithm from Rasmussen & Williams
  # - Algorithms: 3.5, 3.6
  ##############################################################################
  start_time <- Sys.time()
  if(verbose){
    message('Starting: Laplace Approximation for Gaussian Process Classification')
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
    message("Error: Requires class labels -1 and 1 for Algorithm")
    message(paste('Classes found: ', toString(classes)))
    return()
  }

  # Converting target inputs to a column vector
  y = matrix(y, nrow=n_obs, ncol=1)

  #-----------------------------------------------------------------------------
  # Running the Laplace Approximation Algorthim using C++ vvvvvvvvvvvvvvvvvvvv
  #-----------------------------------------------------------------------------
  ests <- la_probit(y, cov_matrix, tol, max_iters)
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
    # message('Approximate log Marginal Likelihood: ', round(results[['log_Z_ep']],3))
    # message(paste('Number of Iterations:',results[['Number_Iters']]))
  }
  results <- list()
  results[['ComputationTime']] = dur
  results[['ps']] <- pnorm(ests)

  return(results)
}


gpbal_la <- function(X, y,
                     cov_function,
                     init_theta,
                     verbose = F,
                     balance_metric = 'mom_sq',
                     wts_vers='ATE'){

  objective_function <- function(theta){
    ncov <- length(theta)
    cov_matrix <- cov_function(as.matrix(X), theta)
    ps_res <- gpbal_la_fixed(y, cov_matrix, verbose = F, tol = 1e-2)
    ps_est <- ps_res$ps

    if(tolower(wts_vers) =='ate'){
      ps_wts <- ifelse(y==1, 1/ps_est, 1/(1-ps_est))
    } else if(tolower(wts_vers) == 'att') {
      ps_wts <- ifelse(y==1, 1, ps_est/(1-ps_est))
    } else {
      message('invalid weighting scheme')
      return(NULL)
    }

    if(balance_metric == 'mom_sq'){
      cb_bal <- .mom_sq_bal(data.frame(X), 1:ncol(X), y==1, ps_wts)
    } else if(balance_metric == 'mom'){
      cb_bal <- .mom_bal(data.frame(X), 1:ncol(X), y==1, ps_wts)
    } else if(balance_metric == 'ks'){
      cb_bal <- 0
      for(i in 1:ncol(X)){
        cb_bal <- cb_bal + .ks_avg_test(X[,i], y, ps_est, 500)
      }
    } else(
      return(NULL)
    )
    return(cb_bal)
  }


  start_time <- Sys.time()
  if(verbose){
    message(paste('Starting Optimization  @  ', start_time))
  }
  opt_theta <- minqa::bobyqa(par = init_theta,
                             fn = objective_function,
                             lower = rep(0, length(init_theta)))
  end_time <- Sys.time()
  if(verbose){

    message(paste('Finished Optimization  @  ', end_time))
    message(paste('Time Difference          :', round(difftime(end_time, start_time, units='secs'), 4)))

    message(paste('Optimal Covariate Balance:', opt_theta$fval))
  }

  opt_matrix <- cov_function(as.matrix(X), opt_theta$par)
  opt_ps <- list()
  opt_ps$ps <- gpbal_la_fixed(y, opt_matrix, verbose = F, tol = 1e-2)$ps
  opt_ps$ComputationTime <- difftime(end_time, start_time, units='secs')
  opt_ps$thetas <- opt_theta$par
  return(opt_ps)
}

