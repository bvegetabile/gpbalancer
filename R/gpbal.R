#' Compute optimally balanced Gaussian process propensity scores
#'
#' @param X Matrix of covariates to be included in the analysis and balanced on
#' @param y Set of observed treatment assignments (y in (0,1))
#' @param cov_function Covariance matrix; for examples, see \code{sqexp} or similar
#' @param tol Tolerance of algorithms.  Difference between the latent scores at each iteration - default 1e-2
#' @param max_iters Maximum number of iterations of algorithm - default 20
#' @param verbose Decision to print progress to screen - default TRUE
#' @param ep_vers Sequential or Parallel EP Algorthim - default \code{parallel}, alternative \code{sequential}
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
#' \item{ \code{thetas} - optimal theta from optimization routine}
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
gpbal <- function(X, y,
                  cov_function,
                  init_theta = NULL,
                  verbose = F,
                  balance_metric = 'mom_sq',
                  approx_method = 'ep',
                  ep_vers = 'parallel',
                  estimand='ATE',
                  eta = 0.05,
                  tol = 1e-2){


  if(eta < 0 | eta > 0.5) stop("Error: Set eta to values in (0, 0.5)")
  if( !(tolower(approx_method) == 'ep' | tolower(approx_method) == 'laplace') ){
    stop('Potential Correction: Set approx_method to "ep" or "laplace"')
  }

  X <- as.matrix(X)

  if(is.null(init_theta)){
    if(is.null(attributes(cov_function)[['ntheta']])){
      stop('See supplied covariance function to learn how many parameters are needed for init_theta')
    } else {
      if(attributes(cov_function)[['ntheta']] == 'ndim+1'){
        init_theta <- rep(1, ncol(X) + 1)
      } else {
        init_theta <- rep(1, attr(cov_function, 'ntheta'))
      }
    }
  }


  objective_function <- function(theta){
    ncov <- length(theta)
    cov_matrix <- cov_function(as.matrix(X), theta)
    ps_res <- gpbal_fixed(y, cov_matrix,
                          verbose = F,
                          tol = tol,
                          approx_method = approx_method,
                          ep_vers=ep_vers)
    ps_est <- as.vector(pnorm(ps_res$PosteriorMean))

    if(tolower(estimand) =='ate'){
      ps_wts <- ifelse(y==1, 1/ps_est, 1/(1-ps_est))
    }
    else if(tolower(estimand) == 'att') {
      ps_wts <- ifelse(y==1, 1, ps_est/(1-ps_est))
    }
    else if(tolower(estimand) == 'tatt'){
      ps_mask <- ifelse((ps_est > eta) & (ps_est < 1-eta), 1, 0)
      ps_wts <- ifelse(y==1, ps_mask, ps_mask * ps_est/(1-ps_est) )
    }
    # else if(tolower(estimand) == 'ato'){
    #   ps_wts <- ifelse(y==1, 1 - ps_est, ps_est)
    # }
    else if(tolower(estimand) == 'matching'){
      min_ps <- apply(cbind(ps_est, 1-ps_est), 1, min)
      ps_wts <- ifelse(y==1, min_ps/ps_est, min_ps/(1-ps_est))
    }
    else {
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
  opt_ps <- gpbal_fixed(y, opt_matrix, verbose = F, tol = tol, approx_method = approx_method, ep_vers=ep_vers)
  opt_ps$ComputationTime <- difftime(end_time, start_time, units='secs')
  opt_ps$init_theta <- init_theta
  opt_ps$thetas <- opt_theta$par
  opt_ps$minval <- opt_theta$fval

  # Final weights
  if(tolower(estimand) =='ate'){
    opt_ps$wts <- ifelse(y==1, 1/opt_ps$ps, 1/(1-opt_ps$ps))
  }
  else if(tolower(estimand) == 'att') {
    opt_ps$wts <- ifelse(y==1, 1, opt_ps$ps/(1-opt_ps$ps))
  }
  else if(tolower(estimand) == 'tatt'){
    ps_mask <- ifelse((opt_ps$ps > eta) & (opt_ps$ps < 1-eta), 1, 0)
    opt_ps$wts <- ifelse(y==1, ps_mask, ps_mask * opt_ps$ps/(1-opt_ps$ps) )
  }
  # else if(tolower(estimand) == 'ato'){
  #   ps_wts <- ifelse(y==1, 1 - ps_est, ps_est)
  # }
  else if(tolower(estimand) == 'matching'){
    min_ps <- apply(cbind(opt_ps$ps, 1-opt_ps$ps), 1, min)
    opt_ps$wts <- ifelse(y==1, min_ps/opt_ps$ps, min_ps/(1-opt_ps$ps))
  }
  else {
    message('invalid weighting scheme')
    return(NULL)
  }

  return(opt_ps)
}
