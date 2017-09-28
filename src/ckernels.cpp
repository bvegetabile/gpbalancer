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

// // [[Rcpp::export]]
// NumericMatrix construct_poly(NumericMatrix X, double p, double sig0, double sig1, double noise) {
//   int N = X.nrow();
//   int D = X.ncol();
//
//   NumericMatrix cov_mat(N,N);
//   double top_sum;
//   double l_sum;
//   double r_sum;
//   double normval;
//   for(int i = 0; i < N; i++){
//     for(int j = i; j < N; j++){
//       top_sum = sig0;
//       l_sum = sig0;
//       r_sum = sig0;
//       for(int d1 = 0; d1 < D; d1++){
//         top_sum += X(i, d1) * X(j, d1);
//         l_sum += X(i, d1) * X(i, d1);
//         r_sum += X(j, d1) * X(j, d1);
//       }
//       normval = top_sum / (sqrt(l_sum)*sqrt(r_sum));
//       cov_mat(i,j) = sig1 * pow(normval,p);
//       cov_mat(j,i) = cov_mat(i,j);
//     }
//   }
//   for(int k = 0; k < N; k++){
//     cov_mat(k,k) += noise;
//   }
//   return(cov_mat);
// }
