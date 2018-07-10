

gpbal.mc.wts <- function(ps, y, estimand='ATE'){
  n_obs <- length(y)
  found_classes <- as.vector(unique(y))
  found_classes <- found_classes[order(found_classes)]
  n_classes <- length(found_classes)

  wts_vec <- matrix(NA, nrow=length(y), ncol=1)
  for(ta in 1:n_classes){
    z_ind <- ifelse(y==found_classes[ta], 1, 0)
    wts_vec <- ifelse(z_ind == 1,
                      (z_ind/ps[,ta]) / sum(z_ind / ps[,ta]),
                      wts_vec)
  }
  if(estimand == 'ATZ'){
    wts_vec[y==which_ta,] <- 1 / sum(y==which_ta)
  }
  wts_vec
}


gpbal.mc <- function(X,
                     y,
                     cov_function,
                     init_theta,
                     return_theta = F,
                     verbose = F,
                     lambda = NULL,
                     n_moments = 2,
                     estimand = 'ATE'){

  n_obs <- length(y)
  found_classes <- as.vector(unique(y))
  found_classes <- found_classes[order(found_classes)]
  n_classes <- length(found_classes)
  targets <- matrix(NA, nrow=n_obs, ncol=n_classes)

  for(c in 1:n_classes){
    targets[,c] <- ifelse(y == found_classes[c], 1, 0)
  }
  colnames(targets) <- found_classes

  if(!estimand %in% c('ATE', 'ATZ')){
    message('Invalid estimand: Choose ATE, or ATZ')
    return(NULL)
  } else {

    if(is.null(lambda)){
      if(estimand == 'ATE'){
        lambda <- rep(0, n_classes)
      } else if (estimand == 'ATZ'){
        lambda <- ifelse(uniq_ta == which_ta, 0, 0)
      }
    } else {
      lambda <- lambda
    }

    if(estimand == "ATE"){
      X_scale <- scale(X, center = T, scale = F)
      Xmat <- optbw::make_Xmat(as.matrix(X_scale), n_moments)
      Xmat <- as.matrix(scale(Xmat, center = T, scale=T))
      bal_targets <- apply(Xmat, 2, mean)
    } else if (estimand == "ATZ"){
      T_mean <- matrix(apply(as.matrix(X[y==which_ta, ]), 2, mean),
                       nrow = nrow(X),
                       ncol= ncol(X),
                       byrow = T)
      X_center <- X - T_mean
      Xmat <- optbw::make_Xmat(as.matrix(X_center), n_moments)
      Xmat <- as.matrix(scale(Xmat, center = T, scale=T))
      bal_targets <- apply(as.matrix(Xmat[y==which_ta, ]), 2, mean)
    }
  }


  objective_function <- function(theta){
    cov_matrix <- cov_function(as.matrix(X), rep(theta, n_classes))
    ps_res <- mcla_fixed(y, cov_matrix, verbose = F, tol = 1e-3)
    wts_mat <- gpbal.mc.wts(ps_res$ps, y, estimand)
    loss <- optbw::targeted_loss(X = Xmat, TI = targets, w = wts_mat, lambda = lambda, targets = bal_targets)
    return(loss)
  }

  start_time <- Sys.time()
  if(verbose){
    message(paste('Starting Optimization  @  ', start_time))
  }
  opt_theta <- minqa::bobyqa(par = init_theta,
                             fn = objective_function,
                             lower = rep(0, length(init_theta)),
                             control=list('maxfun'=200))

  end_time <- Sys.time()
  if(verbose){
    message(paste('Finished Optimization  @  ', end_time))
    message(paste('Time Difference          :', round(difftime(end_time, start_time, units='secs'), 4)))
    message(paste('Optimal Covariate Balance:', opt_theta$fval))
  }

  opt_matrix <- cov_function(as.matrix(X), rep(opt_theta$par, n_classes))
  opt_ps <- mcla_fixed(y, opt_matrix, verbose = F, tol = 1e-6)
  opt_ps$ComputationTime <- difftime(end_time, start_time, units='secs')
  opt_ps$thetas <- opt_theta$par

  return(opt_ps)

}

# covmat <- gpbalancer::mc_sqexp_common(as.matrix(dat), inv_ls_vec = c(0.5,0.5,0.5))
# hmph <- gpbalancer::mcla_fixed(class_label, covmat, verbose = T, tol=1e-8)
# plot(hmph$ps, p)
# abline(0,1)
# system.time(testing1 <- gpbal.mc(as.matrix(dat), class_label,
#                     gpbalancer::mc_sqexp_common, init_theta = 1))
# system.time(testing2 <- optbw::balwts.mc(dat, class_label, lambda = c(0,0,0)))
#
# plot(gpbal.mc.wts(testing1$ps, class_label, 'ATE'), testing2$wts)
# abline(0,1)
