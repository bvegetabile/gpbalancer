// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat mc_sqexp_common(arma::mat X,
                          arma::vec inv_ls_vec,
                          double scale=1.0,
                          double noise = 1e-6){
    int n_classes = inv_ls_vec.n_elem;
    int n_obs = X.n_rows;

    arma::mat cov_mat(n_classes*n_obs, n_classes*n_obs, arma::fill::zeros);
    for(int d = 0; d < n_classes; d++){
        for(int i = 0; i < X.n_rows; i++){
            for(int j = i; j < X.n_rows; j++){
                cov_mat(d*n_obs+i, d*n_obs+j) = pow(scale,2) * std::exp(- 0.5 * arma::sum(arma::pow(X.row(i) - X.row(j), 2) * std::pow(inv_ls_vec[d],2)));
                cov_mat(d*n_obs+j, d*n_obs+i) = cov_mat(d*n_obs+i, d*n_obs+j);
            }
            cov_mat(d*n_obs+i,d*n_obs+i) = cov_mat(d*n_obs+i,d*n_obs+i) + noise;
        }
    }

    return(cov_mat);
}

// [[Rcpp::export]]
arma::mat mc_normpoly_common(arma::mat X,
                             arma::vec sig_shift,
                             arma::vec sig_scale,
                             int power = 1,
                             double noise = 1e-6){
    int n_classes = sig_shift.n_elem;
    int n_obs = X.n_rows;

    arma::mat cov_mat(n_classes*n_obs, n_classes*n_obs, arma::fill::zeros);
    double top_val; double bot_l; double bot_r;
    for(int d = 0; d < n_classes; d++){
        for(int i = 0; i < X.n_rows; i++){
            for(int j = i; j < X.n_rows; j++){
                top_val = arma::dot(X.row(i), X.row(j)) + sig_shift[d];
                bot_l = std::sqrt(arma::dot(X.row(i), X.row(i)) + sig_shift[d]);
                bot_r = std::sqrt(arma::dot(X.row(j), X.row(j)) + sig_shift[d]);
                cov_mat(d*n_obs+i, d*n_obs+j) = sig_scale[d] * std::pow(top_val / bot_l / bot_r, power);
                cov_mat(d*n_obs+j, d*n_obs+i) = cov_mat(d*n_obs+i, d*n_obs+j);
            }
            cov_mat(d*n_obs+i,d*n_obs+i) = cov_mat(d*n_obs+i,d*n_obs+i) + noise;
        }
    }

    return(cov_mat);
}

