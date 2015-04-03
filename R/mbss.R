# Hack based on Hadley Wickhams comment: 
# http://stackoverflow.com/a/12429344/897506
globalVariables(c(# Variables used unquoted inside functions
                 "duration", 
                 "event",
                 "location",
                 "stream",
                 "region",
                 "posterior_prob",
                 "posterior_logprob",
                 "event_posterior",
                 "llr",
                 "rs",
                 "count",
                 "baseline",
                 "overdispersion",
                 "effect_logprob",
                 "duration_condposterior",
                 "num",
                 "denom",
                 "timesum",
                 "locsum",
                 "prop",
                 "duration_event_posterior",
                 "regions_in_part",
                 "totalsum",
                 "loglikelihood",
                 # data.table functions
                 "data.table",
                 "is.data.table",
                 ".",
                 ":=",
                 "setkeyv",
                 ".SD",
                 # foreach
                 "%dopar%",
                 "foreach"),
                package = "scanstatistics")

#' Calculate posterior event probabilities for all regions and event durations.
#' 
#' Calculate the posterior probability of an event for each event type,
#' region (group of one or more locations), and event duration.
#' Assumes that the prior probability of a region being affected by a given
#' event is uniform over all regions, and that events are mutually exclusive
#' and affect only a single region.
#' 
#' @param loglikelihoods A \code{data.table} in long format,
#'        containing at least the columns
#'        \code{event, stream, location, time, null_loglikelihood} and
#'        \code{event_loglikelihood}. The column \code{null_loglikelihood}
#'        contains the predicted log-likelihood for each observation 
#'        (i.e.\ the observation made for a particular data stream, 
#'        at a particular location, at a particular time) 
#'        under the null hypothesis of no event. These values will be duplicated
#'        across the different event types.
#'        The column \code{event_loglikelihood} contains the predicted 
#'        log-likelihood for each observation given that an event of the 
#'        specified type occurs.
#' @param regions A \code{set} of regions, each region itself a set containing
#'        one or more locations of those found in the table 
#'        \code{loglikelihoods}.
#' @param null_prior The prior probability that no event has taken place.
#' @param event_priors The prior probabilities of each event, given as a vector.
#' @param duration_condpriors The prior conditional probabilities of the event
#'        durations, given the event types. Either the string 'uniform',
#'        meaning that the event durations are uniformly distributed for all
#'        event types, or a matrix in which the row numbers correspond to the
#'        event duration, and the columns to the different events.
MBSS <- function(loglikelihoods,
                 regions,
                 null_prior,
                 event_priors,
                 duration_condpriors = "uniform") {
  # Input checking -------------------------------------------------------------
  if (!is.data.table(loglikelihoods)) {
    stop("The log-likelihoods must be supplied as a data.table.")
  }
  if (!all(c("event", "stream", "location", "time", 
             "null_loglikelihood", "event_loglikelihood") %in% 
             names(loglikelihoods))) {
    stop("The data.table containing the log-likelihood must ",
         "contain the columns 'event', 'stream', 'location', 'time', ",
         "'null_loglikelihood', 'event_loglikelihood'.")
  }
  if (any(is.na(loglikelihoods))) {
    stop("No support for missing values yet. Consider imputation methods.")
  }
  if (length(null_prior) != 1 || !is.numeric(null_prior) || null_prior < 0 ||
        null_prior > 1) {
    stop("Prior null hypothesis probability must be a single numeric value ",
         "between 0 and 1.")
  }
  if (!all.equal(1, null_prior + sum(event_priors))) {
    stop("Prior probabilities for events or no events must sum to 1.")
  }
  if (is.character(duration_condpriors) && 
        (length(duration_condpriors) !=1 ||
           duration_condpriors != "uniform")) {
    stop("Specify conditional probabilities for event durations given event ",
         "type either as the string 'uniform' or as a matrix.")
  }
  if (!is.character(duration_condpriors) && 
      !all.equal(rep(1, length(event_priors)), 
                 colSums(duration_condpriors))) {
    stop("Conditional probabilities for event durations given event type ",
         "must sum to the event priors.")
  }
  if (!is.character(duration_condpriors) && 
        ncol(duration_condpriors) != length(event_priors)) {
    stop("The number of event prior probabilities must equal the ",
         "number of events found in the event duration priors.")
  }
  # Proceed with calculations --------------------------------------------------
  times_durations <- times_and_durations(loglikelihoods)
  n_regions <- length(regions)
  n_events <- length(event_priors)
  max_duration <- max(times_durations[, duration])
  
  
  if (is.character(duration_condpriors)) {
    duration_condpriors <- matrix(1 / max_duration, 
                                  nrow = max_duration,
                                  ncol = n_events)
  }
  
  null_loglikelihood <- loglikelihoods[event == loglikelihoods[1, event],  
                                       sum(null_loglikelihood)]

  # Calculate log-likelihood ratios for all events, regions, and durations
  keys_used_for_stllr <- c("region", "event", "duration", "stream", "location")
  spacetime_output <- 
    loglikelihoods[, .(llr = event_loglikelihood - null_loglikelihood),
      keyby = .(location, event, stream, time)] %>%
    add_duration %>%
    region_joiner(regions = regions, keys = keys_used_for_stllr) %>%
    spacetime_llr
  
  # Calculate the marginal probability of the data
  data_logratio <- 
    spacetime_output %>%
    spatial_llr(dur_given_event_logprobs = log(duration_condpriors)) %>%
    data_to_nulldata_logratio(event_logpriors = log(event_priors),
                              null_prior = null_prior,
                              n_regions = n_regions)
  
  marginal_logprob_of_data <- data_logratio + null_loglikelihood

  # Add a column to \code{spacetime_output} for posterior probability
  spacetime_logposterior(spacetime_output,
                         event_logpriors = log(event_priors),
                         dur_given_event_logprobs = log(duration_condpriors),
                         n_regions = n_regions,
                         data_logratio = data_logratio)

  # Calculate probability maps
  event_logpmap <- 
    spacetime_output %>%
    event_logprobability_map(
      region_table = region_table_creator(regions, key = "region"))

  event_pmap <- event_probability_map(event_logpmap)
  pmap <- probability_map(event_logpmap)

  # Calculate the posterior probabilities
  spacetime_output[, posterior_prob := exp(posterior_logprob)]
  null_posterior <- exp(null_loglikelihood + 
                          log(null_prior) - marginal_logprob_of_data)
  event_duration_joint <- posterior_duration_event_jdist(spacetime_output)
  event_posteriors <- posterior_event_probabilities(event_duration_joint)
  duration_condposteriors <- posterior_duration_givn_event(event_duration_joint)

  total_posterior <- null_posterior + event_posteriors[, sum(event_posterior)]
  if (!all.equal(1, total_posterior)) {
    warning("Posterior probabilities don't sum to 1.")
  }
  
  structure(list(loglikelihoods = loglikelihoods,
                 max_duration = max_duration,
                 time_and_durations = times_durations,
                 n_regions = n_regions,
                 regions = regions,
                 null_loglikelihood = null_loglikelihood,
                 null_prior = null_prior,
                 event_priors = event_priors,
                 duration_condpriors = duration_condpriors,
                 posteriors = spacetime_output,
                 null_posterior = null_posterior,
                 event_posteriors = event_posteriors,
                 duration_condposteriors = duration_condposteriors,
                 event_probability_map = event_pmap,
                 probability_map = pmap), 
            class = "MBSS")
}

#' Extract the posterior conditional event duration probabilities from an MBSS
#' object.
#' 
#' @param MBSS_object An object of class "MBSS".
#' @param duration_column A logical scalar. Should the output matrix
#'        contain a column for the event duration (elements being integers
#'        equal to the row numbers)?
get_duration_condposteriors <- function(MBSS_obj, duration_column = TRUE) {
  out <- tidyr::spread(MBSS_obj$duration_condposteriors, 
                       key = event, 
                       value = duration_condposterior)
  out <- as.matrix(out)
  if (duration_column) {
    return(out)
  } else {
    return(out[, -1])
  }
}

#' Get the k most probable region-event type-duration combinations.
#' 
#' Extracts the k most probable combinations of region, event type and event
#' duration. One can also specify a restriction to a subset of the event,
#' and/or a subset of the durations.
get_k_most_probable <- function(MBSS_obj, 
                                k = 10, 
                                events = "all", 
                                durations = "all") {
  if (is.character(events) && events == "all" && 
      is.character(durations) && durations == "all") {
    MBSS_obj$posteriors[order(-posterior_prob), 
      .(region, event, duration, posterior_prob)][seq(k)]
  } else if (is.character(durations) && durations == "all") {
    MBSS_obj$posteriors[event %in% events,
      .(region, event, duration, posterior_prob)][
      order(-posterior_prob)][seq(k)]
  } else if (is.character(events) && events == "all") {
    MBSS_obj$posteriors[duration %in% durations,
      .(region, event, duration, posterior_prob)][order(-posterior_prob)][seq(k)]
  } else {
    MBSS_obj$posteriors[event %in% events & duration %in% durations,
      .(region, event, duration, posterior_prob)][
      order(-posterior_prob)][seq(k)]
  }
}