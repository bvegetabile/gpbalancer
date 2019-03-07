// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>
using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]


arma::vec snd(arma::vec X){
  arma::vec Z = exp(-pow(X, 2)/2) / std::sqrt(arma::datum::pi * 2.0);
  return(Z);
}

arma::vec snc(arma::vec X){
  arma::vec F = 0.5 * ( 1.0 + arma::erf(X / std::sqrt(2.0)));
  return(F);
}

// [[Rcpp::export]]
arma::vec la_probit(arma::vec& targets,
                    arma::mat& covmat,
                    double tol=1e-2,
                    int max_iters = 20,
                    bool verbose = false){
  double n_obs = targets.n_elem;

  arma::vec f(n_obs, arma::fill::zeros);
  arma::vec f_old(n_obs, arma::fill::ones);
  arma::mat W(n_obs, n_obs, arma::fill::zeros);
  arma::mat sqrt_W(n_obs, n_obs, arma::fill::zeros);
  arma::mat big_eye = arma::eye(n_obs, n_obs);
  arma::mat L(n_obs, n_obs);
  // Put loop here
  arma::vec hess_vec(n_obs);
  arma::vec grad_vec(n_obs);
  arma::vec b;
  arma::vec a;

  int i = 1;
  while(arma::mean(arma::abs(f-f_old)) >= tol){
    f_old = f;

    if(verbose) std::cout << "Iter :" << i << "\n";

    hess_vec = arma::pow(snd(f) / snc(targets % f), 2) + targets % f % snd(f) / snc(targets % f);
    grad_vec = targets % snd(f) / snc(targets % f);
    W.diag() = hess_vec;
    sqrt_W.diag() = arma::sqrt(hess_vec);
    L = arma::chol(big_eye + sqrt_W * covmat * sqrt_W, "lower");
    b = W * f + grad_vec;
    a = b - sqrt_W * arma::solve(arma::trimatu(L.t()), arma::solve(arma::trimatl(L), sqrt_W * covmat * b));
    f = covmat * a;
    i++;
  }
  // std::cout << "Mean Difference : " << arma::mean(arma::abs(f-f_old)) << "\n";
  return(f);
}
