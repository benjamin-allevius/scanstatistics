#include <cmath>
#include "probability_functions.h"
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// ZIP log-likelihood (incomplete information) ---------------------------------

//' Calculate a term of the incomplete information loglihood.
//' 
//' Calculate a term of the incomplete information loglihood for the 
//' zero-inflated Poisson distribution.
//' @param y A non-negative integer; the observed count.
//' @param mu A positive scalar; the expected value of the count.
//' @param p A scalar between 0 and 1; the structural zero probability.
//' @param q A scalar greater than or equal to 1; the relative risk.
//' @return A non-positive scalar; the loglihood contribution of the 
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

//' Calculate the incomplete information loglihood.
//'
//' Calculate the incomplete information loglihood for the zero-inflated Poisson 
//' distribution.
//' @param y A non-negative integer vector; the observed counts.
//' @param mu A vector of positive scalars; the expected values of the counts.
//' @param p A vector of scalars between 0 and 1; the structural zero
//'    probabilities.
//' @param q A scalar greater than or equal to 1; the relative risk.
//' @return A non-positive scalar; the incomplete information ZIP loglihood.
//' @keywords internal
// [[Rcpp::export]]
double incomplete_loglihood(const arma::uvec& y,
                            const arma::vec& mu,
                            const arma::vec& p,
                            double q) {
  double loglihood = 0.0;
  for (int i = 0; i < y.n_elem; i++) {
    loglihood += incomplete_loglihood_term(y[i], mu[i], p[i], q);
  }
  return loglihood;
}

// ZIP EM algorithm ------------------------------------------------------------

//' Calculate the conditional expectation of the structural zero indicator.
//' @param mu The expected values of the count (which is zero).
//' @param p The structural zero probability.
//' @param q A scalar greater than or equal to 1; the relative risk.
//' @return A scalar between 0 and 1.
//' @keywords internal
// [[Rcpp::export]]
double estimate_struc_zero(double mu, double p, double q) {
  return p / (p + (1 - p) * exp(-q * mu));
}

//' Estimate the relative risk.
//' @param y_sum A non-negative integer; the sum of the observed counts.
//' @param mu A vector of positive scalars; the expected values of the counts.
//' @param p A vector of scalars between 0 and 1; the structural zero
//'    probabilities.
//' @param d_hat A vector of estimates of the structural zero indicators.
//' @return A scalar; the relative risk.
//' @keywords internal
// [[Rcpp::export]]
double estimate_q(const int y_sum,
                  const arma::vec& mu,
                  const arma::vec& p,
                  const arma::vec& d_hat) {
  double denominator = 0;
  for (int i = 0; i < mu.n_elem; ++i) {
    denominator += mu[i] * (1.0 - d_hat[i]);
  }
  return std::max(1.0, y_sum / denominator);
}

//' Calculate the loglihood ratio statistic and the relative risk.
//' @param y A non-negative integer vector; the observed counts.
//' @param mu A vector of positive scalars; the expected values of the counts.
//' @param p A vector of scalars between 0 and 1; the structural zero
//'    probabilities.
//' @param d_hat A vector of estimates of the structural zero indicators.
//' @param rel_tol The absolute convergence criterion.
//' @return A list with three elements:
//'    \enumerate{
//'      \item The loglihood ratio statistic.
//'      \item The estimated relative risk.
//'      \item The number of iterations of the EM algorithm performed.
//'    } 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List zip_em_algo(const arma::uvec y,
                       const arma::vec& mu,
                       const arma::vec& p,
                       double rel_tol = 1e-2) {
  Rcpp::List score_q_niter (3);
  
  arma::vec d_hat = arma::zeros(y.n_elem); // Structural zero estimates
  double q_hat = 1.0; // Relative risk estimate
  
  double loglik_null = incomplete_loglihood(y, mu, p, q_hat);
  double loglik_old = loglik_null;
  double loglik_new = 0.0;
  
  // Store indices of zero counts; only loop through these when estimating 
  // structural zero indicators. Also compute the sum of all counts.
  std::vector<int> zero_idx;
  int y_sum = 0;
  for (int i = 0; i < y.n_elem; ++i) {
    y_sum += y[i];
    if (y[i] == 0) zero_idx.push_back(i);
  }
  
  // Run EM algorithm
  double diff = 2.0 * rel_tol;
  int n_iterations = 0;
  while (diff > rel_tol) {
    // Expectation-step
    for (const int& i : zero_idx) {
      d_hat[i] = estimate_struc_zero(mu[i], p[i], q_hat);
    }
    // Maximization-step
    q_hat = estimate_q(y_sum, mu, p, d_hat);
    
    // Update likelihood
    loglik_new = incomplete_loglihood(y, mu, p, q_hat);
    diff = abs(exp(loglik_new - loglik_old) - 1);
    loglik_old = loglik_new;
    n_iterations += 1;
  }
  
  score_q_niter[0] = loglik_new - loglik_null;
  score_q_niter[1] = q_hat;
  score_q_niter[2] = n_iterations;
  
  return score_q_niter;
}

// EB ZIP scan statistic -------------------------------------------------------

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