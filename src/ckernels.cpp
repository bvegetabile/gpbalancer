// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]


// [[Rcpp::export]]
NumericMatrix construct_sqexp(NumericMatrix X, NumericMatrix M, double sig, double noise) {
  int N = X.nrow();
  int D = X.ncol();

  NumericMatrix cov_mat(N,N);

  for(int i = 0; i < N; i++){
    for(int j = i + 1; j < N; j++){
      NumericVector diff_vect = X(i, _) - X(j, _);
      double inner_sum = 0;
      for(int d1 = 0; d1 < D; d1++){
        for(int d2 = 0; d2 < D; d2++){
          inner_sum += M(d1,d2) * diff_vect(d1) * diff_vect(d2);
        }
      }
      cov_mat(i,j) = std::pow(sig,2) * std::exp(-inner_sum / 2);
      cov_mat(j,i) = cov_mat(i,j);
    }
  }
  for(int k =0; k < N; k++){
    cov_mat(k,k) = std::pow(sig, 2) + noise;
  }
  return(cov_mat);
}

// [[Rcpp::export]]
arma::mat sqexp_ard(arma::mat X,
                arma::rowvec hyperparams,
                double scale=1.0,
                double noise = 1e-6){
  arma::mat cov_mat(X.n_rows, X.n_rows, arma::fill::zeros);
  for(int i = 0; i < X.n_rows; i++){
    for(int j = i; j < X.n_rows; j++){
      cov_mat(i, j) = scale * exp(- 0.5 * arma::sum(arma::pow(X.row(i) - X.row(j), 2) % arma::pow(hyperparams,2)));
      cov_mat(j, i) = cov_mat(i, j);
    }
    cov_mat(i,i) = cov_mat(i,i) + noise;
  }
  return(cov_mat);
}

// [[Rcpp::export]]
arma::mat sqexp_common(arma::mat X,
                       double lengthscale,
                       double scale=1.0,
                       double noise = 1e-6){
  arma::mat cov_mat(X.n_rows, X.n_rows, arma::fill::zeros);
  for(int i = 0; i < X.n_rows; i++){
    for(int j = i; j < X.n_rows; j++){
      cov_mat(i, j) = scale * std::exp(- 0.5 * arma::sum(arma::pow(X.row(i) - X.row(j), 2) * std::pow(lengthscale,2)));
      cov_mat(j, i) = cov_mat(i, j);
    }
    cov_mat(i,i) = cov_mat(i,i) + noise;
  }
  return(cov_mat);
}

// [[Rcpp::export]]
arma::mat polykernel(arma::mat X,
                     double sig_zero,
                     int pwr = 1,
                     double scale=1.0,
                     double noise = 1e-6){
  arma::mat cov_mat(X.n_rows, X.n_rows, arma::fill::zeros);
  double upper;
  double low_l;
  double low_r;

  for(int i = 0; i < X.n_rows; i++){
    for(int j = i; j < X.n_rows; j++){
      upper = arma::dot(X.row(i), X.row(j)) + std::pow(sig_zero,2);
      low_l = std::sqrt(arma::dot(X.row(i), X.row(i)) + std::pow(sig_zero,2));
      low_r = std::sqrt(arma::dot(X.row(j), X.row(j)) + std::pow(sig_zero,2));
      cov_mat(i, j) = std::pow(scale,2) * std::pow(upper / low_l / low_r, pwr);
      cov_mat(j, i) = cov_mat(i, j);
    }
    cov_mat(i,i) = cov_mat(i,i) + noise;
  }
  return(cov_mat);
}
