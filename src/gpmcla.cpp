// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;

// [[Rcpp::export]]
List gp_mcla(arma::mat covmat,
             arma::vec targets,
             int n_classes,
             double tol = 1e-10,
             int max_iters = 20,
             bool verbose = false){
    int N = targets.n_elem;
    int obs_per_class = N / n_classes;
    // Initializing latent scores
    arma::vec f(N, arma::fill::zeros);
    arma::vec f_old(N, arma::fill::ones);

    // Initializing probabilities
    arma::uvec classes(n_classes);
    std::iota(classes.begin(), classes.end(), 0);
    arma::uvec class_loc = classes * obs_per_class;

    // Initializing some other things
    arma::vec pi_vec(N);
    arma::vec prob_i(n_classes);
    arma::mat pi_mat(N, obs_per_class, arma::fill::zeros);
    arma::mat E(N,N, arma::fill::zeros);
    arma::mat Lc(obs_per_class, obs_per_class, arma::fill::zeros);
    arma::mat Dc(obs_per_class, obs_per_class, arma::fill::zeros);
    arma::mat M(obs_per_class, obs_per_class, arma::fill::zeros);
    arma::mat Kc(obs_per_class, obs_per_class);
    arma::vec z(n_classes);
    arma::mat I_n(obs_per_class, obs_per_class, arma::fill::eye);
    arma::vec b_vec(N);
    arma::vec c_vec(N);
    arma::vec a_vec(N);
    arma::mat R(N, obs_per_class, arma::fill::zeros);

    arma::mat out_pi(obs_per_class, n_classes, arma::fill::zeros);

    for(int i = 0; i < obs_per_class; i++){
        R(i*N + class_loc + i) += 1.0;
    }


    int start_loc;
    int end_loc;

    // This is where there while loop goes...

    int counter = 0;
    while(arma::mean(arma::pow(f_old - f, 2)) > tol){
        if(counter == max_iters){
            std::cout << "Max Iterations Reached" << "\n";
            break;
        }
        f_old = f;
        M.fill(0);

        for(int i = 0; i < obs_per_class; i++){
            prob_i = f.elem(i + class_loc);
            pi_vec(i + class_loc) = arma::exp(prob_i) / arma::sum(arma::exp(prob_i));
            pi_mat(i*N + class_loc + i) = arma::exp(prob_i) / arma::sum(arma::exp(prob_i));
        }


        for(int c = 0; c < n_classes; c++){
            start_loc = c*obs_per_class;
            end_loc = (c+1)*obs_per_class - 1;
            Kc = covmat.submat(start_loc, start_loc, end_loc, end_loc);
            Dc = arma::sqrt(arma::diagmat(pi_vec.subvec(start_loc, end_loc)));
            Lc = arma::chol(I_n + Dc * Kc * Dc, "lower");
            E.submat(start_loc,
                     start_loc,
                     end_loc,
                     end_loc) = Dc * arma::solve(arma::trimatu(Lc.t()),
                                                 arma::solve(arma::trimatl(Lc), Dc));
            z(c) = arma::sum(Lc.diag());
            M += E.submat(start_loc, start_loc, end_loc, end_loc);
        }
        M = arma::chol(M, "lower");
        b_vec = (arma::diagmat(pi_vec) - pi_mat * pi_mat.t()) * f + targets - pi_vec;
        c_vec = E * covmat * b_vec;

        a_vec = b_vec - c_vec + E * R * arma::solve(arma::trimatu(M.t()),
                                                    arma::solve(arma::trimatl(M),
                                                                R.t() * c_vec));

        f = covmat * a_vec;
        counter += 1;
        if(verbose){
            std::cout << "Iteration : " << counter << "\n";
        }
    }

    for(int i = 0; i < obs_per_class; i++){
        prob_i = f.elem(i + class_loc);
        out_pi(i + class_loc) = arma::exp(prob_i) / arma::sum(arma::exp(prob_i));
    }

    double log_marg_lik = -0.5 * arma::dot(a_vec, f) + arma::dot(targets, f) - arma::sum(z);
    for(int i = 0; i < obs_per_class; i++){
        log_marg_lik += std::log(arma::sum(arma::exp(f.elem(i + class_loc))));
    }

    return(List::create(_["ps"] = out_pi,
                        _["logmarglik"] = log_marg_lik));
}
//
// // [[Rcpp::export]]
// List gp_mcla_fast(arma::mat covmat,
//              arma::vec targets,
//              int n_classes,
//              double tol = 1e-10,
//              int max_iters = 20,
//              bool verbose = false){
//     Timer timer;
//     timer.step("start");
//
//
//     int N = targets.n_elem;
//     int obs_per_class = N / n_classes;
//     // Initializing latent scores
//     arma::vec f(N, arma::fill::zeros);
//     arma::vec f_old(N, arma::fill::ones);
//
//     // Initializing probabilities
//     arma::uvec classes(n_classes);
//     std::iota(classes.begin(), classes.end(), 0);
//     arma::uvec class_loc = classes * obs_per_class;
//
//     // Initializing some other things
//     arma::vec pi_vec(N);
//     arma::vec prob_i(n_classes);
//     arma::mat pi_mat(N, obs_per_class, arma::fill::zeros);
//     arma::mat E(N,obs_per_class, arma::fill::zeros);
//     arma::mat Lc(obs_per_class, obs_per_class, arma::fill::zeros);
//     arma::mat Dc(obs_per_class, obs_per_class, arma::fill::zeros);
//     arma::mat M(obs_per_class, obs_per_class, arma::fill::zeros);
//     arma::mat Kc(obs_per_class, obs_per_class);
//     arma::vec z(n_classes);
//     arma::mat I_n(obs_per_class, obs_per_class, arma::fill::eye);
//     arma::vec b_vec(N);
//     arma::vec c_vec(N);
//     arma::vec a_vec(N);
//     arma::mat R(N, obs_per_class, arma::fill::zeros);
//
//     arma::mat out_pi(obs_per_class, n_classes, arma::fill::zeros);
//
//
//     for(int i = 0; i < obs_per_class; i++){
//         R(i*N + class_loc + i) += 1.0;
//     }
//
//
//     int start_loc;
//     int end_loc;
//
//     // This is where there while loop goes...
//
//     timer.step("init_setup");
//
//     int counter = 0;
//     timer.step("loop_start");
//     while(arma::mean(arma::pow(f_old - f, 2)) > tol){
//
//         if(counter == max_iters){
//             std::cout << "Max Iterations Reached" << "\n";
//             break;
//         }
//         f_old = f;
//         M.fill(0);
//
//         for(int i = 0; i < obs_per_class; i++){
//             prob_i = f.elem(i + class_loc);
//             pi_vec(i + class_loc) = arma::exp(prob_i) / arma::sum(arma::exp(prob_i));
//             pi_mat(i*N + class_loc + i) = arma::exp(prob_i) / arma::sum(arma::exp(prob_i));
//         }
//         timer.step("prop_calc:" + std::to_string(counter+1));
//
//         b_vec = (arma::diagmat(pi_vec) - pi_mat * pi_mat.t()) * f + targets - pi_vec;
//         for(int c = 0; c < n_classes; c++){
//             start_loc = c*obs_per_class;
//             end_loc = (c+1)*obs_per_class - 1;
//             Kc = covmat.submat(start_loc, start_loc, end_loc, end_loc);
//             Dc = arma::sqrt(arma::diagmat(pi_vec.subvec(start_loc, end_loc)));
//             Lc = arma::chol(I_n + Dc * Kc * Dc, "lower");
//             E.rows(start_loc, end_loc) = Dc * arma::solve(arma::trimatu(Lc.t()), arma::solve(arma::trimatl(Lc), Dc));
//             z(c) = arma::sum(Lc.diag());
//             M += E.rows(start_loc, end_loc);
//         }
//         M = arma::chol(M, "lower");
//         timer.step("E/M_calc:" + std::to_string(counter+1));
//         // b_vec = (arma::diagmat(pi_vec) - pi_mat * pi_mat.t()) * f + targets - pi_vec;
//
//         for(int c = 0; c < n_classes; c++){
//             start_loc = c*obs_per_class;
//             end_loc = (c+1)*obs_per_class - 1;
//             c_vec.subvec(start_loc, end_loc) = E.rows(start_loc, end_loc) * covmat.submat(start_loc, start_loc, end_loc, end_loc) * b_vec.subvec(start_loc, end_loc);
//         }
//         timer.step("c_calc:" + std::to_string(counter+1));
//         a_vec = b_vec - c_vec + E * arma::solve(arma::trimatu(M.t()),
//                                                     arma::solve(arma::trimatl(M),
//                                                                 R.t() * c_vec));
//         timer.step("a_calc:" + std::to_string(counter+1));
//         f = covmat * a_vec;
//         counter += 1;
//         if(verbose){
//             std::cout << "Iteration : " << counter << "\n";
//         }
//     }
//
//     timer.step("while_loop");
//
//     for(int i = 0; i < obs_per_class; i++){
//         prob_i = f.elem(i + class_loc);
//         out_pi(i + class_loc) = arma::exp(prob_i) / arma::sum(arma::exp(prob_i));
//     }
//
//     double log_marg_lik = -0.5 * arma::dot(a_vec, f) + arma::dot(targets, f) - arma::sum(z);
//     for(int i = 0; i < obs_per_class; i++){
//         log_marg_lik += std::log(arma::sum(arma::exp(f.elem(i + class_loc))));
//     }
//     timer.step("log_like_calc");
//
//     NumericVector res(timer);   //
//     for (int i=0; i<res.size(); i++) {
//         res[i] = res[i];
//     }
//
//     return(List::create(_["ps"] = out_pi,
//                         _["logmarglik"] = log_marg_lik,
//                         _["some_times"] = res));
// }
//
//
//
