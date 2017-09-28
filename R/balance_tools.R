################################################################################
#
# Tools for performing causal analysis
# --- Author: Brian Vegetabile (bvegetab [AT] uci [DOTEDU])
#
################################################################################

#-------------------------------------------------------------------------------
# Estimation Tools for estimating treatment effects
# --- ht_est calculates the average treatment effect by Horvitz Thompson Estimation
#
#
# Balance Tools for assessing covariate balance
# --- bal_stats calculates balance statistics (wtd or not) for single variable
# --- bal_table calculates balance statistics for dataset
# --- wtd_ecdf provides a function to calculate the a weighted empirical CDF
# --- dens_table calculates the probability density in the treatment groups
# --- bal_plt_cat_pdf plots the density function in the treatment groups to compare
# --- bal_plt_cont_pdf plots the conditional distribution of treatment groups
# --- bal_plt_cont_cdf plots the conditional emp. cdf for treatment groups
#-------------------------------------------------------------------------------

bal_stats <- function(var_data,
                      treat_ind,
                      datatype = 'continuous',
                      wts=rep(1, length(treat_ind))){
  #-----------------------------------------------------------------------------
  # bal_stats is a function which is used to assess covariate balance for a
  # single variable.
  #
  # Input variables
  # ---  var_data  : vector of data for a single variable
  # ---  treat_ind : vector of treatment indicators.  Binary or Logical
  # ---  datatype  : type of data provided.  allows 'continuous' or 'binary'
  # ---  wts       : weights after estimating propensity score or similar.
  #                  defaults to 1 for all observations unless otherwise
  #                  specified.
  #-----------------------------------------------------------------------------

  # Converting 0/1's to logicals to ensure subsetting works as desired
  treat_ind <- as.logical(treat_ind)

  #-----------------------------------------------------------------------------
  # Weighted mean in the treated and control groups

  # Weighted sample mean in treated group
  xbar_t <- sum(wts[treat_ind] * var_data[treat_ind]) / sum(wts[treat_ind])
  # Weighted sample mean in control group
  xbar_c <- sum(wts[!treat_ind] * var_data[!treat_ind]) / sum(wts[!treat_ind])

  #-----------------------------------------------------------------------------
  # Weighted sample variance calculations in both groups for binary and
  # continuous data.  Number of samples in each group also calculated.

  if(datatype=='binary'){
    # Binary Sample Variance - Treatment Group
    s2_t <- xbar_t*(1-xbar_t)
    # Binary Sample Variance - Control Group
    s2_c <- xbar_c*(1-xbar_c)

    # Sample sizes for binary data
    NT <- sum(as.logical(var_data) & treat_ind)
    NC <- sum(as.logical(var_data) & !treat_ind)

  } else if(datatype=='continuous') {
    # Weighted Sample Variance in Treatment Group
    s2_t <- sum(wts[treat_ind] * (var_data[treat_ind]- xbar_t)**2) / sum(wts[treat_ind])
    # Weighted Sample Variance in Control Group
    s2_c <- sum(wts[!treat_ind] * (var_data[!treat_ind]- xbar_c)**2) / sum(wts[!treat_ind])

    # Sample sizes for continuous data
    NT <- sum(treat_ind)
    NC <- sum(!treat_ind)
  }

  #-----------------------------------------------------------------------------
  # Calculating relevant statistics.
  # -- Standardized Difference in Means for Balance.  Want |std_diff| < 0.1
  # -- Log of ratio of Standard Deviation.  Checks balance of second moment.
  #    Should be close to zero.

  std_diff <- (xbar_t - xbar_c) / sqrt((s2_t + s2_c) / 2)
  log_var_ratio <- log(sqrt(s2_t)) - log(sqrt(s2_c))

  return(c(NT, round(xbar_t,4), round(s2_t,4),
           NC, round(xbar_c,4), round(s2_c,4),
           round(std_diff,4), round(log_var_ratio, 4)))
}


bal_table <- function(dataset,
                      col_ind,
                      treat_ind,
                      wts = rep(1, length(treat_ind)),
                      max_uniq=5,
                      plot_balance = FALSE){
  #-----------------------------------------------------------------------------
  # bal_table is a function which provides a table of covariate balance stats.
  #
  # Input variables
  # --- dataset   : matrix or dataframe of data
  # --- col_ind   : indices of columns to check balance. vector of integers
  # --- treat_ind : indicator of treatment assignment. binary or logical
  # --- wts       : weights for use after estimating the propensity score
  # --- max_uniq  : defines threshold for what to consider a continuous or
  #                 categorical variable.  If the number of unique members of a
  #                 category is less than this number it will be converted to a
  #                 factor if not already a factor.
  #-----------------------------------------------------------------------------
  var_names <- colnames(dataset)[col_ind]

  treat_ind <- as.logical(treat_ind)
  outtable <- c()
  counter <- 1

  if(plot_balance){
    nplots <- length(col_ind)
    if(nplots <= 3){
      par(mfrow=c(1, nplots))
    } else{
      nr <- ceiling(sqrt(nplots))
      nc <- ceiling(sqrt(nplots))
      par(mfrow=c(nr, nc))
    }
  }

  # Iteration based on the order of column numbers provided
  for(i in 1:length(col_ind)){
    c <- col_ind[i]
    col_data <- dataset[, c]
    # Conversion of variables to categorical if it has less than the max_uniq
    # categories.
    if(length(unique(col_data)) <= max_uniq){
      col_data <- as.factor(col_data)
    }
    if(is.factor(col_data)){
      obs_lvls <- levels(col_data)
      outtable <- rbind(outtable, rep(NA, 8))
      row.names(outtable)[counter] <- var_names[i]
      counter <- counter + 1
      for(lvl in obs_lvls){
        stddiff <- bal_stats(col_data==lvl, treat_ind, 'binary', wts)
        outtable <- rbind(outtable, stddiff)
        row.names(outtable)[counter] <- lvl
        counter <- counter + 1
      }
      if(plot_balance){
        bal_plt_cat_pdf(col_data, treat_ind, toptitle = var_names[i], vertoffset = 0.5)
      }
    } else {
      stddiff <- bal_stats(col_data, treat_ind, 'continuous', wts)
      if(plot_balance){
        bal_plt_cont_cdf(col_data, treat_ind, wts = wts,
                         toptitle = var_names[i], var_name = var_names[i])
      }
      outtable <- rbind(outtable, stddiff)
      row.names(outtable)[counter] <- var_names[i]
      counter <- counter + 1
    }
  }
  colnames(outtable) <- c('NT', 'MeanT', 'VarT',
                          'NC', 'MeanC', 'VarC',
                          'StdDiff', 'LogRatio')
  return(outtable)
}


.mom_bal <- function(dataset,
                     col_ind,
                     treat_ind,
                     wts = rep(1, length(treat_ind)),
                     max_uniq=5){

  treat_ind <- as.logical(treat_ind)
  outtable <- c()
  outbal <- 0.0
  counter <- 1
  # Iteration based on the order of column numbers provided
  for(i in 1:length(col_ind)){
    c <- col_ind[i]
    col_data <- dataset[, c]
    # Conversion of variables to categorical if it has less than the max_uniq
    # categories.
    if(length(unique(col_data)) <= max_uniq){
      col_data <- as.factor(col_data)
    }
    if(is.factor(col_data)){
      obs_lvls <- levels(col_data)
      outtable <- c()
      for(lvl in obs_lvls){
        stddiff <- bal_stats(col_data==lvl, treat_ind, 'binary', wts)
        outtable <- rbind(outtable, stddiff)
      }
      outbal <- outbal + sum(abs(outtable[2:nrow(outtable),7])) + sum(abs(outtable[2:nrow(outtable),8]))
    } else {
      stddiff <- bal_stats(col_data, treat_ind, 'continuous', wts)
      outbal <- outbal + abs(stddiff[7]) + abs(stddiff[8])
    }
  }
  return(outbal)
}

.mom_sq_bal <- function(dataset,
                        col_ind,
                        treat_ind,
                        wts = rep(1, length(treat_ind)),
                        max_uniq=5){

  treat_ind <- as.logical(treat_ind)
  outtable <- c()
  outbal <- 0.0
  counter <- 1
  # Iteration based on the order of column numbers provided
  for(i in 1:length(col_ind)){
    c <- col_ind[i]
    col_data <- dataset[, c]
    # Conversion of variables to categorical if it has less than the max_uniq
    # categories.
    if(length(unique(col_data)) <= max_uniq){
      col_data <- as.factor(col_data)
    }
    if(is.factor(col_data)){
      obs_lvls <- levels(col_data)
      outtable <- c()
      for(lvl in obs_lvls){
        stddiff <- bal_stats(col_data==lvl, treat_ind, 'binary', wts)
        outtable <- rbind(outtable, stddiff)
      }
      outbal <- outbal + sum(abs(outtable[2:nrow(outtable),7])^2) + sum(abs(outtable[2:nrow(outtable),8])^2)
    } else {
      stddiff <- bal_stats(col_data, treat_ind, 'continuous', wts)
      outbal <- outbal + abs(stddiff[7]^2) + abs(stddiff[8]^2)
    }
  }
  return(outbal)
}

.wtd_ecdf <- function (var_data, wts) {
  #-----------------------------------------------------------------------------
  # wtd_ecdf is a modification of the ecdf() function in base R.  It modifies
  # the function to be able to incorporate weights.  This is to visualize
  # balance using the empirical cumulative distribution function for continuous
  # covariates after weighting by the inverse of the propensity score (IPTW)
  #
  # Input variables
  # --- var_data : covariate values - vector of data
  # --- wts      : weights for assessing cov balance by IPTW - vector of data.
  #-----------------------------------------------------------------------------
  ord <- order(var_data)
  var_ordered <- var_data[ord]
  wts_ordered <- wts[ord]
  n <- length(var_data)
  if (n < 1)
    stop("'var_data' must have 1 or more non-missing values")
  vals <- unique(var_ordered)
  matched_vals <- match(var_ordered, vals)
  weight_list <- aggregate(wts_ordered, by=list(matched_vals), sum)
  rval <- approxfun(vals, cumsum(weight_list[,2])/sum(wts_ordered),
                    method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}

.dens_table <- function(var_data, treat_data){
  #-----------------------------------------------------------------------------
  # dens_table creates a contingency table of a covariate and treatment.
  # - Provides the conditional density based upon treatment assignment after
  #   normalizing by the number of treated and control observations
  #
  # Input Variables
  # --- var_data   : vector of a categorical variable
  # --- treat_data : vector of treatment data
  #-----------------------------------------------------------------------------
  nc <- length(unique(var_data))
  nt <- length(unique(treat_data))
  outro <- table(treat_data, var_data,
                 dnn=c('Treatment Levels','Covariate Levels'))
  rs <- rowSums(outro)
  return(outro / matrix(rs,nrow=nt,ncol=nc))
}
#
# bal_plt_cat_pdf <- function(var_data,
#                             treat_data,
#                             spread=0.5,
#                             vertoffset=0.05,
#                             toptitle="Conditional Distribution of Covariate",
#                             legendpos='topright',
#                             col_0 = rgb(0,0,0.75,0.75),
#                             col_1 = rgb(0,0.75,0,0.75)){
#   discrete_dens <- dens_table(var_data, treat_data)
#   col_list <- c(col_0, col_1)
#   nt <- nrow(discrete_dens)
#   t_names <- row.names(discrete_dens)
#
#   ncov <- ncol(discrete_dens)
#   cov_names <- colnames(discrete_dens)
#
#   xspots <- seq(0.5, ncov-0.5, 1)
#   xdiffs <- seq(0, spread, length.out = nt) - spread/2
#   ymax <- min(1.1, max(discrete_dens)+vertoffset)
#
#   plot(0, xlim=c(0, ncov), ylim=c(0, ymax),
#        pch=19, col=rgb(0,0,0,0.0),
#        xaxs="i", yaxs='i',
#        xlab="Category", ylab='Density', axes=F,
#        main=toptitle)
#   for(i in 1:nt){
#     for(j in 1:ncov){
#       lines(rep(xspots[j] + xdiffs[i],2),
#             c(0, discrete_dens[i,j]),
#             lty=1, col=col_list[i])
#     }
#     points(xspots + xdiffs[i],
#            discrete_dens[i,],
#            pch=19, col=col_list[i])
#     text(xspots + xdiffs[i],
#          discrete_dens[i,],
#          labels = round(discrete_dens[i,],2), pos=3)
#   }
#   abline(v=0:ncov)
#   axis(1, xspots, cov_names)
#   axis(2, seq(0,1,0.1), las=1)
#   box()
#   legend(legendpos, legend = c('Control', 'Treated'),
#          lwd=2, lty=1, col=col_list)
# }
#
# bal_plt_cont_pdf <- function(var_data,
#                              treat_data,
#                              var_name = 'Covariate Name',
#                              toptitle = "Conditional Distributions of Covariate",
#                              legendpos = 'topright',
#                              treat_names = c('Treated', 'Control'),
#                              adj_bw = 1,
#                              col_c = rgb(0,0,0.75,0.75),
#                              col_t = rgb(0,0.75,0,0.75)){
#   treat_data <- as.logical(treat_data)
#   var_range <- range(var_data)
#   var_pts <- seq(var_range[1], var_range[2], 0.005)
#   dens_t <- approxfun(density(var_data[treat_data],
#                               adjust = adj_bw,
#                               from = var_range[1], to = var_range[2]))
#   dens_c <- approxfun(density(var_data[!treat_data],
#                               adjust = adj_bw,
#                               from = var_range[1], to = var_range[2]))
#   y_max <- max(dens_t(var_pts), dens_c(var_pts))
#   y_max <- y_max + y_max/100
#   plot(var_pts, dens_t(var_pts),
#        xlim = var_range, ylim=c(0, y_max),
#        xlab=var_name, ylab = 'Density',
#        type='l', col=col_t, lwd=3,
#        main=toptitle)
#   abline(v=var_range, lty=3)
#   lines(var_pts, dens_c(var_pts),
#        col=col_c, lwd=3, xlab=var_name,
#        main=toptitle)
#   legend(legendpos, treat_names, lty=1, lwd=3, col=c(col_t, col_c))
# }
#
# bal_plt_cont_cdf <- function(var_data,
#                              treat_data,
#                              wts = rep(1, length(treat_data)),
#                              var_name = 'Covariate Name',
#                              toptitle = "Conditional Emp. CDF of Covariate",
#                              legendpos = 'bottomright',
#                              treat_names = c('Treated', 'Control'),
#                              adj_bw = 1,
#                              col_c = rgb(0,0,0.75,0.75),
#                              col_t = rgb(0,0.75,0,0.75)){
#   treat_data <- as.logical(treat_data)
#   var_range <- range(var_data)
#   var_pts <- seq(var_range[1], var_range[2], 0.005)
#   wecdf_t <- wtd_ecdf(var_data[treat_data], wts[treat_data])
#   wecdf_c <- wtd_ecdf(var_data[!treat_data], wts[!treat_data])
#   plot(var_pts, wecdf_t(var_pts),
#        xlim = var_range, ylim=c(0, 1),
#        xlab=var_name, ylab = 'Cumulative Distribution Function',
#        type='l', col=col_t, lwd=3,
#        main=toptitle)
#   abline(h=c(0,1), lty=3)
#   abline(v=var_range, lty=3)
#   lines(var_pts, wecdf_c(var_pts),
#         col=col_c, lwd=3, xlab=var_name,
#         main=toptitle)
#   legend(legendpos, treat_names, lty=1, lwd=3, col=c(col_t, col_c))
# }


.ks_avg_test <- function(X, TA, ps=rep(1/length(X), length(X)), n_pts=1000){
  xmin <- min(X)
  xmax <- max(X)
  int_pts <- seq(xmin, xmax, length.out = n_pts)
  t_wts <- (TA / ps) / sum(TA / ps)
  c_wts <- ((1-TA) / (1-ps)) / sum(((1-TA) / (1-ps)))
  t_fn <- .wtd_ecdf(X, wts = t_wts)
  c_fn <- .wtd_ecdf(X, wts = c_wts)
  return(mean(abs(t_fn(int_pts) - c_fn(int_pts))))
}


.ks_test <- function(X, TA, ps=rep(1/length(X), length(X)), n_pts=1000){
  xmin <- min(X)
  xmax <- max(X)
  int_pts <- seq(xmin, xmax, length.out = n_pts)
  t_wts <- (TA / ps) / sum(TA / ps)
  c_wts <- ((1-TA) / (1-ps)) / sum(((1-TA) / (1-ps)))
  t_fn <- .wtd_ecdf(X, wts = t_wts)
  c_fn <- .wtd_ecdf(X, wts = c_wts)
  return(max(abs(t_fn(int_pts) - c_fn(int_pts))))
}

