// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace RcppParallel;

//
struct ConstructSqExpKernel : public Worker {
  /* Inputs:
   * - design_x  : design matrix for distance computations
   * - vec_ls    : vector of inverse length scale parameters
   * - sig_var   : outer varianve component
   * - sig_noise : diagonal variance component for numerical stability
   */
  const arma::mat& design_x;
  const arma::vec& vec_theta;
  const double& sig_noise;

  /* Outputs:
   * - cov_mat : covariance matrix from kernel functions
   */
  arma::mat& cov_mat;


  // initialize
  ConstructSqExpKernel(const arma::mat& design_x,
                       const arma::vec& vec_theta,
                       const double& sig_noise,
                       arma::mat& cov_mat)
    : design_x(design_x),
      vec_theta(vec_theta),
      sig_noise(sig_noise),
      cov_mat(cov_mat) {}

  // function call
  void operator()(std::size_t row_beg, std::size_t row_end) {
    double x_dim = design_x.n_cols;
    double sig_var = vec_theta[0];

    for(unsigned int i = row_beg; i < row_end; i++){
      for(unsigned int j = 0; j <= i; j++){
        double inner_sum = 0.0;
        for(int d = 0; d < x_dim; d++){
          inner_sum += pow(vec_theta[d+1], 2) * pow(design_x(i, d) - design_x(j,d), 2);
        }

        cov_mat(i, j) = pow(sig_var,2) * exp( - inner_sum / 2.0);
        cov_mat(j, i) = cov_mat(i, j);
        if(i == j){
          cov_mat(i, i) += sig_noise;
        }
      }
    }
  }
};

// [[Rcpp::export]]
arma::mat par_sqexp(arma::mat design_x,
                    arma::vec vec_theta,
                    double sig_noise = 1e-6) {

  // allocate the matrix we will return
  arma::mat cov_mat(design_x.n_rows, design_x.n_rows);

  // create the worker
  ConstructSqExpKernel constructSqExpKernel(design_x,
                                            vec_theta,
                                            sig_noise,
                                            cov_mat);

  // call it with parallelFor
  parallelFor(0, design_x.n_rows, constructSqExpKernel);

  return cov_mat;
}
