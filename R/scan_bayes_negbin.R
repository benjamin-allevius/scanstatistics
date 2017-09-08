

#' Calculate the negative binomial bayesian scan statistic..
#' 
#' Calculate the "Bayesian Spatial Scan Statistic" by Neill et al. (2006),
#' adapted to a spatio-temporal setting. The scan statistic assumes that,
#' given the relative risk, the data follows a Poisson distribution. The 
#' relative risk is in turn assigned a Gamma distribution prior, yielding a 
#' negative binomial marginal distribution for the counts under the null 
#' hypothesis. Under the alternative hypothesis, the 
#' @param counts Either:
#'    \itemize{
#'      \item A matrix of observed counts. Rows indicate time and are ordered
#'            from least recent (row 1) to most recent (row 
#'            \code{nrow(counts)}). Columns indicate locations, numbered from 1 
#'            and up. If \code{counts} is a matrix, the optional matrix argument
#'            \code{baselines} should also be specified.
#'      \item A data frame with columns "time", "location", "count", "baseline".
#'            Alternatively, the column "baseline" can be replaced by a column
#'            "population". The baselines are the expected values of the counts.
#'    }
#' @param zones A list of integer vectors. Each vector corresponds to a single
#'    zone; its elements are the numbers of the locations in that zone.
#' @param baselines Optional. A matrix of the same dimensions as \code{counts}. 
#'    Not needed if \code{counts} is a data frame. Holds the Poisson mean 
#'    parameter for each observed count. Will be estimated if not supplied 
#'    (requires the \code{population} argument). These parameters are typically 
#'    estimated from past data using e.g. Poisson (GLM) regression.
#' @param population Optional. A matrix or vector of populations for each 
#'    location. Not needed if \code{counts} is a data frame. If \code{counts} is
#'    a matrix, \code{population} is only needed if \code{baselines} are to be 
#'    estimated and you want to account for the different populations in each 
#'    location (and time). If a matrix, should be of the same dimensions as 
#'    \code{counts}. If a vector, should be of the same length as the number of 
#'    columns in \code{counts}.
#' @param outbreak_prob A scalar; the probability of an outbreak (at any time,
#'    any place). Defaults to 0.05.
#' @param alpha_null A scalar; the shape parameter for the gamma distribution
#'    under the null hypothesis of no anomaly. Defaults to 1.
#' @param beta_null A scalar; the scale parameter for the gamma distribution
#'    under the null hypothesis of no anomaly. Defaults to 1.
#' @param alpha_alt A scalar; the shape parameter for the gamma distribution
#'    under the alternative hypothesis of an anomaly. Defaults to the same value
#'    as \code{alpha_null}.
#' @param beta_alt A scalar; the scale parameter for the gamma distribution
#'    under the alternative hypothesis of an anomaly. Defaults to the same value
#'    as \code{beta_null}.
#' @param inc_values A vector of possible values for the increase in the mean
#'    (and variance) of an anomalous count. Defaults to evenly spaced values 
#'    between 1 and 3, with a difference of 0.1 between consecutive values.
#' @param inc_probs A vector of the prior probabilities of each value in 
#'    \code{inc_values}. Defaults to 1, implying a discrete uniform 
#'    distribution.
#' @return A list which, in addition to the information about the type of scan
#'    statistic, has the following components: \code{priors} (list), 
#'    \code{posteriors} (list), \code{MLC} (list) and \code{marginal_data_prob} 
#'    (scalar). The list \code{MLC} has elements
#'    \describe{
#'      \item{zone}{The number of the spatial zone of the most likely cluster 
#'                  (MLC).}
#'      \item{duration}{The most likely event duration.}
#'      \item{log_posterior}{The posterior log probability that an event is 
#'                           ongoing in the MLC.}
#'      \item{log_bayes_factor}{The logarithm of the Bayes factor for the MLC.}
#'      \item{posterior}{The posterior probability that an event is ongoing in 
#'                       the MLC.}
#'      \item{locations}{The locations involved in the MLC.}
#'    }    
#'    The list \code{priors} has elements
#'    \describe{
#'      \item{null_prior}{The prior probability of no anomaly.}
#'      \item{alt_prior}{The prior probability of an anomaly.}
#'      \item{inc_prior}{A vectorof prior probabilities of each value in the 
#'                       argument \code{inc_values}.}
#'      \item{window_prior}{The prior probability of an outbreak in any of the
#'                          space-time windows.}
#'    }
#'    The list \code{posteriors} has elements
#'    \describe{
#'      \item{null_posterior}{The posterior probability of no anomaly.}
#'      \item{alt_posterior}{The posterior probability of an anomaly.}
#'      \item{inc_posterior}{A data frame with columns \code{inc_values} and
#'                           \code{inc_posterior}.}
#'      \item{window_posteriors}{A data frame with columns \code{zone}, 
#'                               \code{duration}, \code{log_posterior} and 
#'                               \code{log_bayes_factor}, each row corresponding
#'                               to a space-time window.}
#'      \item{space_time_posteriors}{A matrix with the posterior anomaly 
#'                                   probability of each location-time 
#'                                   combination.}
#'      \item{location_posteriors}{A vector with the posterior probability of an
#'                                 anomaly at each location.}
#'    }
#' @references 
#'    Neill, D. B., Moore, A. W., Cooper, G. F. (2006). 
#'    \emph{A Bayesian Spatial Scan Statistic}. Advances in Neural Information 
#'    Processing Systems 18.
#' @importFrom dplyr arrange
#' @importFrom magrittr %<>%
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' # Create location coordinates, calculate nearest neighbors, and create zones
#' n_locs <- 50
#' max_duration <- 5
#' n_total <- n_locs * max_duration
#' geo <- matrix(rnorm(n_locs * 2), n_locs, 2)
#' knn_mat <- coords_to_knn(geo, 15)
#' zones <- knn_zones(knn_mat)
#'
#' # Simulate data
#' baselines <- matrix(rexp(n_total, 1/5), max_duration, n_locs)
#' counts <- matrix(rpois(n_total, as.vector(baselines)), max_duration, n_locs)
#'
#' # Inject outbreak/event/anomaly
#' ob_dur <- 3
#' ob_cols <- zones[[10]]
#' ob_rows <- max_duration + 1 - seq_len(ob_dur)
#' counts[ob_rows, ob_cols] <- matrix(
#'   rpois(ob_dur * length(ob_cols), 2 * baselines[ob_rows, ob_cols]), 
#'   length(ob_rows), length(ob_cols))
#' res <- scan_bayes_negbin(counts = counts,
#'                          zones = zones,
#'                          baselines = baselines)
#' }
scan_bayes_negbin <- function(counts,
                              zones,
                              baselines = NULL,
                              population = NULL,
                              outbreak_prob = 0.05,
                              alpha_null = 1,
                              beta_null = 1,
                              alpha_alt = alpha_null,
                              beta_alt = beta_null,
                              inc_values = seq(1, 3, by = 0.1),
                              inc_probs = 1) {
  if (is.data.frame(counts)) {
    # Validate input -----------------------------------------------------------
    if (any(c("time", "location", "count") %notin% names(counts))) {
      stop("Data frame counts must have columns time, location, count")
    }
    if (all(c("baseline", "population") %notin% names(counts))) {
      stop("For data frame counts, either population or baseline must be ",
           "specified.")
    }
    counts %<>% arrange(location, -time)
    # Create matrices ----------------------------------------------------------
    if ("baseline" %notin% names(counts)) {
      baselines <- NULL
      population <- df_to_matrix(counts, "time", "location", "population")
    } else {
      baselines <- df_to_matrix(counts, "time", "location", "baseline")
    }
    counts <- df_to_matrix(counts, "time", "location", "count")
  }
  
  if (length(inc_probs) == 1) {
    inc_probs <- rep(1, length(inc_values))
  }
  inc_probs <- inc_probs / sum(inc_probs)
  
  # Validate input -------------------------------------------------------------
  if (any(as.vector(counts) != as.integer(counts))) {
    stop("counts must be integer")
  }
  if (!is.null(baselines) && any(baselines <= 0)) {
    stop("baselines must be positive")
  }
  if (length(inc_values) != length(inc_probs)) {
    stop("inc_probs must be either of length 1 or length(inc_values).")
  }
  if (any(inc_probs <= 0)) {
    stop("inc_probs must contain positive values only.")
  }
  if (any(c(alpha_null, beta_null, alpha_alt, beta_alt) <= 0)) {
    stop("The alpha and beta parameters must be positive")
  }
  if (outbreak_prob <= 0 || outbreak_prob >= 1) {
    stop("It must hold that 0 < outbreak_prob < 1.")
  }
  
  # Reshape into matrices ------------------------------------------------------
  if (is.vector(counts)) {
    counts <- matrix(counts, nrow = 1)
  }
  if (!is.null(baselines) && is.vector(baselines)) {
    baselines <- matrix(baselines, nrow = 1)
  }

  # Estimate baselines if not supplied -----------------------------------------
  if (is.null(baselines)) {
    baselines <- estimate_baselines(counts, population)
  } 
  
  # Reverse time order: most recent first --------------------------------------
  counts <- flipud(counts)
  baselines <- flipud(baselines)

  # Prepare zone arguments for C++ ---------------------------------------------
  args <- list(counts = counts, 
               baselines = baselines,
               zones = unlist(zones) - 1, 
               zone_lengths = unlist(lapply(zones, length)),
               outbreak_prob = outbreak_prob,
               alpha_null = alpha_null,
               beta_null = beta_null,
               alpha_alt = alpha_alt,
               beta_alt = beta_alt,
               inc_values = inc_values,
               inc_probs = inc_probs)
  
  # Run analysis on observed counts --------------------------------------------
  scan <- do.call(scan_bayes_negbin_cpp, args)
  
  #  Rename and reshape some elements
  scan$priors$inc_prior <- as.vector(scan$priors$inc_prior)
  
  scan$posteriors$inc_prior <- as.vector(scan$posteriors$inc_prior)
  
  scan$posteriors$location_posteriors <- 
    as.vector(scan$posteriors$location_posteriors)
  
  # Sort posterior probabilities
  scan$posteriors$window_posteriors %<>% arrange(-log_posterior)
  
  MLC <- as.list(scan$posteriors$window_posteriors[1, ])
  MLC$posterior <- exp(MLC$log_posterior)
  MLC$locations <- zones[[MLC$zone]]

  structure(
    c(list(# General
           distribution = "negative binomial (Gamma-Poisson mixture)",
           type = "Bayesian",
           setting = "univariate"),
      # MLC + analysis
      MLC = list(MLC),
      scan,
      # Configuration
      list(n_zones = length(zones),
           n_locations = ncol(counts),
           max_duration = nrow(counts))),
    class = "scanstatistic")
}
