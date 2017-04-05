

#' Calculate the expectation-based ZIP scan statistic.
#' @param counts A matrix of observed counts. Rows indicate time and are ordered
#'    from most recent (row 1) to least recent. Columns indicate locations, 
#'    numbered from 1 and up.
#' @param zones A list of integer vectors. Each vector corresponds to a single
#'    zone; its elements are the numbers of the locations in that zone.
#' @param baselines A matrix of the same dimensions as \code{counts}. Holds the
#'    Poisson mean parameter of the ZIP distribution for each observed count.
#' @param probs A matrix of the same dimensions as \code{counts}. Holds the
#'    structural zero probability of the ZIP distribution for each observed 
#'    count.
#' @param population A matrix or vector of populations for each location. Only
#'    needed if \code{baselines} and \code{probs} are to be estimated and you 
#'    want to account for the different populations in each location (and time).
#'    If a matrix, should be of the same dimensions as \code{counts}. If a 
#'    vector, should be of the same length as the number of columns in 
#'    \code{counts}.
#' @param n_mcsim A non-negative integer; the number of replicate scan 
#'    statistics to generate in order to calculate a P-value.
#' @param gumbel Boolean. If \code{TRUE}, and \code{n_mcsim > 0}, then a Gumbel
#'    distribution is fit to the replicate scan statistics and a \eqn{p}-value
#'    for the observed scan statistic is calculated using the fitted 
#'    distribution. If \code{FALSE} and \code{n_mcsim > 0}, the empirical
#'    \eqn{p}-value is calculated instead, defined as 1 + the number of 
#'    replicate scan statistics that were larger than the observed scan 
#'    statistic, divided by \code{1 + n_mcsim}.
#' @param max_only Boolean. If \code{FALSE} (default) the log-likelihood ratio
#'    statistic for each zone and duration is returned. If \code{TRUE}, only the
#'    largest such statistic (i.e. the scan statistic) is returned, along with
#'    the corresponding zone and duration.
#' @param rel_tol A positive scalar. If the relative change in the incomplete
#'    information likelihood is less than this value, then the EM algorithm is
#'    deemed to have converged.
#' @details For the expectation-based zero-inflated Poisson scan statistic 
#'    (Kjellson 2015), the null hypothesis of no anomaly holds that the count 
#'    observed at each location \eqn{i} and duration \eqn{t} (the number of time 
#'    periods before present) has a zero-inflated Poisson distribution with 
#'    expected value parameter \eqn{\mu_{it}} and structural zero probability
#'    \eqn{p_{it}}:
#'    \deqn{
#'      H_0 : Y_{it} \sim \textrm{ZIP}(\mu_{it}, p_{it}).
#'    }
#'    This holds for all locations \eqn{i = 1, \ldots, m} and all durations 
#'    \eqn{t = 1, \ldots,T}, with \eqn{T} being the maximum duration considered.
#'    Under the alternative hypothesis, there is a space-time window \eqn{W}
#'    consisting of a spatial zone \eqn{Z \subset \{1, \ldots, m\}} and a time 
#'    window \eqn{D \subseteq \{1, \ldots, T\}} such that the counts in that
#'    window have their Poisson expected value parameters inflated by a factor 
#'    \eqn{q_W > 1} compared to the null hypothesis:
#'    \deqn{
#'    H_1 : Y_{it} \sim \textrm{ZIP}(q_W \mu_{it}, p_{it}), ~~(i,t) \in W.
#'    }
#'    For locations and durations outside of this window, counts are assumed to
#'    be distributed as under the null hypothesis. The sets \eqn{Z} considered 
#'    are those specified in the argument \code{zones}, while the maximum 
#'    duration \eqn{T} is taken as the maximum value in the column 
#'    \code{duration} of the input \code{table}. 
#'    
#'    For each space-time window \eqn{W} considered, (the log of) a likelihood 
#'    ratio is computed using the distributions under the alternative and null 
#'    hypotheses, and the expectation-based Poisson scan statistic is calculated 
#'    as the maximum of these quantities over all space-time windows. The 
#'    expectation-maximization (EM) algorithm is used to obtain maximum 
#'    likelihood estimates. Point estimates of the parameters \eqn{\mu_{it}} 
#'    must be specified in the column \code{mu} of the argument \code{table} 
#'    before this function is called.
#' @references 
#'    Kjellson, B. (2015), \emph{Spatiotemporal Outbreak Detection: A Scan 
#'    Statistic Based on the Zero-Inflated Poisson Distribution}, (Master 
#'    Thesis, Stockholm University),
#'    \href{http://goo.gl/6Q89ML}{Link to PDF}.
#' @examples 
#' \dontrun{
#' set.seed(1)
#' # Create location coordinates, calculate nearest neighbors, and create zones
#' geo <- matrix(rnorm(100), 50, 2)
#' knn_mat <- coords_to_knn(geo, 15)
#' zones <- knn_zones(knn_mat)
#' 
#' # Simulate data
#' baselines <- matrix(rexp(50, 1/5), 5, 50)
#' probs <- matrix(runif(50) / 4, 5, 50)
#' counts <- gamlss.dist::rZIP(1, baselines, probs)
#' 
#' # Inject outbreak/event/anomaly
#' ob_dur <- 1:3
#' ob_zone <- zones[[10]]
#' counts[ob_dur, ob_zone] <- gamlss.dist::rZIP(
#'   1, 2 * baselines[ob_dur, ob_zone], probs[ob_dur, ob_zone])
#' res <- scan_zip_eb(counts = counts, 
#'                    zones = zones,
#'                    baselines = baselines, 
#'                    probs = probs,
#'                    n_mcsim = 100,
#'                    gumbel = TRUE,
#'                    max_only = FALSE,
#'                    rel_tol = 1e-3)
#' }
scan_zip_eb <- function(counts,
                        zones,
                        baselines = NULL,
                        probs = NULL,
                        population = NULL,
                        n_mcsim = 0,
                        gumbel = TRUE, 
                        max_only = FALSE,
                        rel_tol = 1e-3) {
  
  if (is.null(baselines) || is.null(probs)) {
    pars <- estimate_zip_params(counts, population)
    baselines <- pars$baselines
    probs <- pars$probs
  }
  
  zones_flat <- unlist(zones) - 1
  zone_lengths <- unlist(lapply(zones, length))
  
  if (max_only) {
    scan <- calc_one_zip_eb(counts, baselines, probs, 
                            zones_flat, zone_lengths,
                            rel_tol)
  } else {
    scan <- calc_all_zip_eb(counts, baselines, probs, 
                            zones_flat, zone_lengths,
                            rel_tol)
  }
  
  MLC <- scan[which.max(scan$score), ]
  
  # Monte Carlo replications of the scan statistic under null hypothesis
  repl_stat <- numeric(n_mcsim)
  for (i in seq_len(n_mcsim)) {
    repl_stat[i] <- calc_one_zip_eb(
      gamlss.dist::rZIP(nrow(counts) * ncol(counts), baselines, probs), 
      baselines, 
      probs, 
      zones_flat, zone_lengths,
      rel_tol)$score
  }
  
  # Fit Gumbel distribution to Monte Carlo replicates
  gumbel_mu <- NA
  gumbel_sigma <- NA
  if (n_mcsim > 0 && gumbel) {
    gum_fit <- ismev::gum.fit(repl_stat, show = FALSE)
    gumbel_mu <- gum_fit$mle[1]
    gumbel_sigma <- gum_fit$mle[2]
  }
  
  # Get P-value
  pvalue <- NA
  if (n_mcsim > 0) {
    if (gumbel) {
      pvalue <- reliaR::pgumbel(MLC$score, gumbel_mu, gumbel_sigma, 
                                lower.tail = FALSE)
    } else {
      pvalue <- mc_pvalue(MLC$score, repl_stat)
    }
  }
  
  list(MLC = list(zone_number = MLC$zone,
                  locations = zones[[MLC$zone]],
                  duration = MLC$duration,
                  score = MLC$score,
                  relative_risk = MLC$relrisk,
                  EM_iterations = MLC$n_iter,
                  observed = counts[seq_len(MLC$duration), 
                                    zones[[MLC$zone]]],
                  baselines = baselines[seq_len(MLC$duration), 
                                        zones[[MLC$zone]]],
                  probs = probs[seq_len(MLC$duration), 
                                zones[[MLC$zone]]]),
       table = scan,
       replicate_statistics = repl_stat,
       gumbel_mu = gumbel_mu,
       gumbel_sigma = gumbel_sigma,
       pvalue = pvalue)
}
