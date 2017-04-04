#include "probability_functions.h"
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// EM-algorithm ----------------------------------------------------------------

//' Calculate a term of the incomplete information log likelihood.
//' 
//' Calculate a term of the incomplete information log likelihood for the 
//' zero-inflated Poisson distribution.
//' @param y A non-negative integer; the observed count.
//' @param mu A positive scalar; the expected value of the count.
//' @param p A scalar between 0 and 1; the structural zero probability.
//' @param q A scalar greater than or equal to 1; the relative risk.
//' @return A non-positive scalar; the log-likelihood contribution of the 
//'    observation.
//' @keywords internal
// [[Rcpp::export]]
double incomplete_loglihood_term(int y, double mu, double p, double q) {
  if (y == 0) {
    return log(p + (1 - p) * exp(-q * mu));
  } else {
    return log(1 - p) + y * log(q * mu) - lgamma(y + 1.0) - q * mu;
  }
}

// 
// //' Calculate the incomplete information log likelihood.
// //' 
// //' Calculate the incomplete information log likelihood for the zero-inflated 
// //' Poisson distribution.
// //' @param y A non-negative integer vector; the observed counts.
// //' @param mu A vector of positive scalars; the expected values of the counts.
// //' @param p A vector of scalars between 0 and 1; the structural zero 
// //'    probabilities.
// //' @param q A scalar greater than or equal to 1; the relative risk.
// //' @return A non-positive scalar; the incomplete information ZIP 
// //'    log-likelihood.
// //' @keywords internal
// // [[Rcpp::export]]
// double incomplete_loglihood(const arma::uvec& y,
//                             const arma::vec& mu,
//                             const arma::vec& p,
//                             double q) {
//   double loglihood = 0.0;
//   for (int i = 0; i < y.n_elem; i++) {
//     loglihood += incomplete_loglihood_term(y[i], 
//                                            mu[i], 
//                                            p[i], 
//                                            q);
//   }
//   return loglihood;
// }
// 
// // [[Rcpp::export]]
// double score_zip_eb(const arma::umat counts,
//                     const arma::vec& baselines,
//                     const arma::vec& probs, 
//                     bool complete_info = false) {
//   int n_obs = counts.n_elem;
//   
//   
//   double score;
//   return score;
// }
// 
// // [[Rcpp::export]]
// arma::mat calc_all_zip_eb(const arma::umat& counts,
//                           const arma::mat& baselines,
//                           const arma::mat& probs,
//                           const std::vector<IntegerVector>& zones) {
//   int max_duration = counts.n_rows;
//   int n_zones = zones.size();
//   arma::mat scores (n_zones, max_duration);
//   
//   for (int d = 0; d < max_duration; ++d) {
//     for (int i = 0; i < n_zones; ++i) {
//       scores[i, d] = score_zip(
//         arma::vectorise(counts.submat(arma::span(0, d), zones[i])),
//         arma::vectorise(baselines.submat(arma::span(0, d), zones[i])),
//         arma::vectorise(probs.submat(arma::span(0, d), zones[i])));
//     }
//   }
//   return scores;
// }