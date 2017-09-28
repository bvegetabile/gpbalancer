// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]


arma::vec std_norm_density(arma::vec X){
  arma::vec Z = exp(-pow(X, 2)/2) / std::sqrt(arma::datum::pi * 2.0);
  return(Z);
}

arma::vec std_norm_cdf(arma::vec X){
  arma::vec F = 0.5 * ( 1.0 + arma::erf(X / std::sqrt(2.0)));
  return(F);
}
//
  double norm_dens(double x){
    double res = std::exp(-std::pow(x, 2.0) / 2.0)/ std::sqrt(arma::datum::pi * 2.0);
    return(res);
  }

double norm_cdf(double x){
  double res = 0.5 * ( 1.0 + std::erf(x / std::sqrt(2.0)));
  return(res);
}


struct UpdateSiteParameters : public Worker {
  /* Needed Terms:
    * - random_order: shuffled data
  * - tilde_nu:
    * - tilde_tau:
    * - mu: Current Posterior mean
  * - sigma_mat: Current Posterior Covariance
  */
    const arma::ivec& random_order;
  const arma::vec& y;
  const arma::vec& old_nu;
  const arma::vec& old_tau;
  const arma::vec& mu;
  const arma::mat& sigma_mat;


  /* Things we want to update, that do not touch other loop throughs...
  * cavity_tau
  * cavity_nu
  * small_z
  * big_z
  *
    * final site parameters: tilde_nu, tilde_tau
  */
    arma::vec& cavity_tau;
  arma::vec& cavity_nu;
  arma::vec& hat_mu;
  arma::vec& hat_sig;
  arma::vec& small_z;
  arma::vec& big_z;
  arma::vec& tilde_nu;
  arma::vec& tilde_tau;
  /*
    * Output Variables
  */

    UpdateSiteParameters(const arma::ivec& random_order,
                         const arma::vec& y,
                         const arma::vec& old_nu,
                         const arma::vec& old_tau,
                         const arma::vec& mu,
                         const arma::mat& sigma_mat,
                         arma::vec& cavity_tau,
                         arma::vec& cavity_nu,
                         arma::vec& hat_mu,
                         arma::vec& hat_sig,
                         arma::vec& small_z,
                         arma::vec& big_z,
                         arma::vec& tilde_nu,
                         arma::vec& tilde_tau)
  : random_order(random_order),
  y(y),
  old_nu(old_nu),
  old_tau(old_tau),
  mu(mu),
  sigma_mat(sigma_mat),
  cavity_tau(cavity_tau),
  cavity_nu(cavity_nu),
  hat_mu(hat_mu),
  hat_sig(hat_sig),
  small_z(small_z),
  big_z(big_z),
  tilde_nu(tilde_nu),
  tilde_tau(tilde_tau) {}
  // function call
  void operator()(std::size_t i_beg, std::size_t i_end) {
    int ind;
    double var_ratio;
    double normal_ratio;
    double delta_tau;

    for(int i = i_beg; i < i_end; i++){
      // Getting the random index from the shuffled list of indices
      ind = random_order(i);

      //# Computing Approximate Cavity Parameters

        cavity_tau(ind) = (1.0/sigma_mat(ind, ind)) - tilde_tau(ind);
      cavity_nu(ind) = (1.0/sigma_mat(ind, ind)) * mu(ind) - tilde_nu(ind);

      //# Compute Posterior Marginal Moments - hat terms from algorithm 3.58
        //# -- Converted into terms of natural cavity parameters
        var_ratio = 1.0 / (cavity_tau(ind) * std::sqrt(1.0 + 1.0/cavity_tau(ind)));
      small_z(ind) = y(ind) * cavity_nu(ind) * var_ratio;
      big_z(ind) = norm_cdf(small_z(ind));

      normal_ratio = norm_dens(small_z(ind)) / norm_cdf(small_z(ind));

      hat_mu(ind) = cavity_nu(ind)/cavity_tau(ind) + y(ind) * normal_ratio * var_ratio;
      hat_sig(ind) = (1.0/cavity_tau(ind)) - normal_ratio * std::pow(var_ratio, 2.0) * (small_z(ind) + normal_ratio);

      // Updating the site parameters
      delta_tau = (1.0/hat_sig(ind)) - cavity_tau(ind) - tilde_tau(ind);
      tilde_tau(ind) = tilde_tau(ind) + delta_tau;
      tilde_nu(ind) = hat_mu(ind) / hat_sig(ind) - cavity_nu(ind);
    }
  }
};


// [[Rcpp::export]]
List par_ep(arma::vec y,
           arma::mat cov_matrix,
           double tol,
           int max_iters,
           bool verbose) {

  /* allocate things...
  *
    */
    RNGScope rngScope;

  // Finding the total number of observations to initialize all vectors
  int n_obs = y.size();
  const arma::mat big_eye = arma::eye(n_obs, n_obs);

  //Initialization of natural site parameters
  arma::vec tilde_nu(n_obs, arma::fill::zeros);
  arma::vec tilde_tau(n_obs, arma::fill::zeros);

  // Initialization of Posterior Values
  arma::vec mu(n_obs, arma::fill::zeros);
  arma::mat sigma_mat = cov_matrix;

  // Initialization of hat parameters
  arma::vec small_z(n_obs, arma::fill::zeros);
  arma::vec big_z = std_norm_cdf(small_z);
  arma::vec hat_mu(n_obs, arma::fill::zeros);
  arma::vec hat_sig(n_obs, arma::fill::zeros);

  // Initializatio of cavity parameters
  arma::vec cavity_tau(n_obs, arma::fill::zeros);
  arma::vec cavity_nu(n_obs, arma::fill::zeros);

  // Setting Initial Values for convergence criteria
  double max_change_tau = 100;
  double max_change_nu = 100;

  // double percent_change_tau = 1;
  // double percent_change_nu = 1;



  // Initializing count for iterations
  int iters = 0;

  if(verbose){
    REprintf(".");
  }

  /*
    * Initialization outside of loop for efficiency;
  */

    arma::mat S_tilde_onehalf;
  arma::mat L;
  arma::mat V;
  arma::vec old_tau;
  arma::vec old_nu;


  while((max_change_nu > tol or max_change_tau > tol) and iters < max_iters){
    //# while(iters < max_iters){
      // Creating a Random Order to Calculate Approximations
    arma::ivec random_order(n_obs);
    std::iota(random_order.begin(), random_order.end(), 0);
    std::random_shuffle(random_order.begin(), random_order.end());

    old_tau = tilde_tau;
    old_nu = tilde_nu;

    UpdateSiteParameters updateSiteParameters(random_order,
                                              y,
                                              old_nu,
                                              old_tau,
                                              mu,
                                              sigma_mat,
                                              cavity_tau,
                                              cavity_nu,
                                              hat_mu,
                                              hat_sig,
                                              small_z,
                                              big_z,
                                              tilde_nu,
                                              tilde_tau);

    parallelFor(0, n_obs, updateSiteParameters);

    S_tilde_onehalf = diagmat(arma::sqrt(tilde_tau));
    L = arma::chol(big_eye + S_tilde_onehalf * cov_matrix * S_tilde_onehalf, "lower");
    V = arma::solve(trimatl(L), S_tilde_onehalf*cov_matrix);

    sigma_mat = cov_matrix - V.t() * V;
    mu = sigma_mat * tilde_nu;

    max_change_nu = max(abs(old_nu - tilde_nu));
    max_change_tau = max(abs(old_tau - tilde_tau));
    iters += 1;
    if(verbose){
      REprintf(".");
    }
  }

  arma::vec cavity_mu = cavity_nu / cavity_tau;
  arma::mat T_mat = diagmat(cavity_tau);
  arma::mat S_tilde = diagmat(tilde_tau);

  double terms_4_1 = 0.5*sum(log(1.0 + tilde_tau/cavity_tau)) - trace(log(L));
  double terms_2_5quad = 0.5 * dot(tilde_nu, (sigma_mat - inv_sympd(T_mat + S_tilde))*tilde_nu);
  double terms_5remain = 0.5 * dot(cavity_mu, T_mat * inv_sympd(S_tilde + T_mat) * (S_tilde*cavity_mu + 2.0*tilde_nu));
  double terms_3 = sum(log(big_z));
  //
    double log_Z_ep = terms_4_1 + terms_2_5quad + terms_5remain + terms_3;
  return(List::create(
    _["Number_Iters"] = iters,
    _["PosteriorMean"] = mu,
    _["PosteriorVar"] = sigma_mat,
    _["tilde_nu"] = tilde_nu,
    _["tilde_tau"] = tilde_tau,
    _["log_Z_ep"] = log_Z_ep));
};



// -----------------------------------------------------------------------------
// Sequential EP

// [[Rcpp::export]]
List seq_ep(arma::vec y,
            arma::mat cov_matrix,
            double tol,
            int max_iters,
            bool verbose) {
  RNGScope rngScope;

  // Finding the total number of observations to initialize all vectors
  int n_obs = y.size();
  arma::mat big_eye = arma::eye(n_obs, n_obs);

  //Initialization of natural site parameters
  arma::vec tilde_nu(n_obs, arma::fill::zeros);
  arma::vec tilde_tau(n_obs, arma::fill::zeros);

  // Initialization of Posterior Values
  arma::vec mu(n_obs, arma::fill::zeros);
  arma::mat sigma_mat = cov_matrix;

  // Initialization of hat parameters
  arma::vec small_z(n_obs, arma::fill::zeros);
  arma::vec big_z(n_obs, arma::fill::zeros);
  arma::vec hat_mu(n_obs, arma::fill::zeros);
  arma::vec hat_sig(n_obs, arma::fill::zeros);

  // Initializatio of cavity parameters
  arma::vec cavity_tau(n_obs, arma::fill::zeros);
  arma::vec cavity_nu(n_obs, arma::fill::zeros);

  // Setting Initial Values for convergence criteria
  double max_change_tau = 100;
  double max_change_nu = 100;

  // double percent_change_tau = 1;
  // double percent_change_nu = 1;

  // Initializing count for iterations
  int iters = 0;

  if(verbose){
    REprintf(".");
  }

  //----------------------------------------------------------------------------
  // Initialization of variables for loop

  arma::vec old_tau; arma::vec old_nu;
  int ind; double var_ratio; double normal_ratio; double delta_tau;
  double needed_constant;
  // arma::mat sigma_diff;
  arma::mat S_tilde_onehalf; arma::mat L; arma::mat V;
  arma::vec cavity_mu; arma::mat T_mat; arma::mat S_tilde;
  // arma::vec column_of_sigma;
  //----------------------------------------------------------------------------

  while((max_change_nu > tol or max_change_tau > tol) and iters < max_iters){
    //# while(iters < max_iters){
    // Creating a Random Order to Calculate Approximations
    arma::ivec random_order(n_obs);
    std::iota(random_order.begin(), random_order.end(), 0);
    std::random_shuffle(random_order.begin(), random_order.end());

    old_tau = tilde_tau;
    old_nu = tilde_nu;

    for(int i = 0; i < n_obs; i++){
      // Getting the random index from the shuffled list of indices
      ind = random_order(i);

      //# Computing Approximate Cavity Parameters

      cavity_tau(ind) = (1.0/sigma_mat(ind, ind)) - tilde_tau(ind);
      cavity_nu(ind) = (1.0/sigma_mat(ind, ind)) * mu(ind) - tilde_nu(ind);

      //# Compute Posterior Marginal Moments - hat terms from algorithm 3.58
      //# -- Converted into terms of natural cavity parameters
      var_ratio = 1.0 / (cavity_tau(ind) * std::sqrt(1.0 + 1.0/cavity_tau(ind)));
      small_z(ind) = y(ind) * cavity_nu(ind) * var_ratio;
      big_z(ind) = norm_cdf(small_z(ind));

      normal_ratio = norm_dens(small_z(ind)) / norm_cdf(small_z(ind));

      hat_mu(ind) = cavity_nu(ind)/cavity_tau(ind) + y(ind) * normal_ratio * var_ratio;
      hat_sig(ind) = (1.0/cavity_tau(ind)) - normal_ratio * std::pow(var_ratio, 2.0) * (small_z(ind) + normal_ratio);

      // Updating the site parameters
      delta_tau = (1.0/hat_sig(ind)) - cavity_tau(ind) - tilde_tau(ind);
      tilde_tau(ind) = tilde_tau(ind) + delta_tau;
      tilde_nu(ind) = hat_mu(ind) / hat_sig(ind) - cavity_nu(ind);

      needed_constant = 1.0/((1.0/delta_tau) + sigma_mat(ind, ind));
      sigma_mat -= needed_constant * (sigma_mat.col(ind) * sigma_mat.col(ind).t());
      mu = sigma_mat * tilde_nu;
    }

    S_tilde_onehalf = diagmat(arma::sqrt(tilde_tau));
    L = chol(big_eye + S_tilde_onehalf * cov_matrix * S_tilde_onehalf, "lower");
    V = solve(trimatl(L), S_tilde_onehalf*cov_matrix);

    sigma_mat = cov_matrix - V.t() * V;
    mu = sigma_mat * tilde_nu;
    max_change_nu = max(abs(old_nu - tilde_nu));
    max_change_tau = max(abs(old_tau-tilde_tau));
    iters += 1;
    if(verbose){
      REprintf(".");
    }
  }

  // # Calculating the approximate marginal likelihood q(y|X)--Equation 3.65-------
  // #  - Terms 4 and 1 of the equation

  // S_tilde_onehalf = diagmat(arma::sqrt(tilde_tau));
  // L = chol(big_eye + S_tilde_onehalf * cov_matrix * S_tilde_onehalf, "lower");
  // V = solve(trimatl(L), S_tilde_onehalf*cov_matrix);
  //
  // sigma_mat = cov_matrix - V.t() * V;
  // mu = sigma_mat * tilde_nu;

  cavity_mu = cavity_nu / cavity_tau;
  T_mat = diagmat(cavity_tau);
  S_tilde = diagmat(tilde_tau);

  double terms_4_1 = 0.5*sum(log(1.0 + tilde_tau/cavity_tau)) - trace(log(L));
  double terms_2_5quad = 0.5 * dot(tilde_nu, (sigma_mat - inv_sympd(T_mat + S_tilde))*tilde_nu);
  double terms_5remain = 0.5 * dot(cavity_mu, T_mat * inv_sympd(S_tilde + T_mat) * (S_tilde*cavity_mu + 2.0*tilde_nu));
  double terms_3 = sum(log(big_z));
  //
  double log_Z_ep = terms_4_1 + terms_2_5quad + terms_5remain + terms_3;
  return(List::create(
      _["Number_Iters"] = iters,
      _["PosteriorMean"] = mu,
      _["PosteriorVar"] = sigma_mat,
      _["tilde_nu"] = tilde_nu,
      _["tilde_tau"] = tilde_tau,
      _["log_Z_ep"] = log_Z_ep));
}


