#' Calculate posterior event probabilities for all regions and event durations.
#' 
#' Calculate the posterior probability of an event for each event type,
#' region (group of one or more locations), and event duration.
#' Assumes that the prior probability of a region being affected by a given
#' event is uniform over all regions, and that events are mutually exclusive
#' and affect only a single region.
#' @param loglikelihoods A \code{data.table} in long format,
#'    containing at least the columns \code{event, stream, location, time} and 
#'    \code{loglikelihood}. See details below for full specification.
#' @param regions A \code{list} or \code{set} of regions, 
#'    each region itself a set containing one or more locations of those found 
#'    in \code{loglikelihoods}.
#' @param null_prior The prior probability that no event has taken place.
#' @param event_priors The prior probabilities of each event type, or numbers 
#'    that express how often the events occur compared to each other. For 
#'    example, \code{c(0.02, 0.05)} and \code{c(2, 5)} express the same thing, 
#'    given the \code{null_prior}. If event types are given as strings, should 
#'    be a named vector or \code{data.frame} with the different events as 
#'    columns.
#' @param duration_condpriors The prior conditional probabilities of the event
#'    durations, given the event types. Either the string 'uniform', meaning 
#'    that the event durations are uniformly distributed for all event types, 
#'    or a matrix in which the row numbers correspond to the event duration, and 
#'    the columns to the different events. If events types are given as strings, make it a \code{data.frame} with
#'    the events types as column names instead.
#' @param n_partitions Integer parameter that specifies how many parts, 
#'    approximately equal in size in the number of locations, the set of all
#'    regions should be partitioned into to avoid memory issues. The more 
#'    regions, and size of these regions, the larger this parameter should be.
#'    Can be left as \code{"auto"} (default) to be figured out for itself.
#' @details 
#'    More details on the input arguments:
#'    \itemize{
#'      \item{loglikelihoods}
#'        {The column \code{event} contains identifiers (integer or character) 
#'        for both the null hypothesis of no event and for the event types
#'        under the alternative hypothesis. If integer, the null hypothesis
#'        should be specified as \code{0L}, and if character, 
#'        should be specified as \code{"null"}.
#'        The column \code{loglikelihood} contains the predicted log-likelihood 
#'        for each observation and each event type (including the null).}
#'    }
#' @importFrom magrittr %>%
#' @export
MBSS <- function(loglikelihoods,
                 regions,
                 null_prior,
                 event_priors,
                 duration_condpriors = "uniform",
                 n_partitions = "auto") {
  # Input checking -------------------------------------------------------------
  if (!is.data.table(loglikelihoods)) {
    stop("The log-likelihoods must be supplied as a data.table.")
  }
  if (!all(c("event", "stream", "location", "time", "loglikelihood") %in% 
             names(loglikelihoods))) {
    stop("The data.table containing the log-likelihoods must contain the ",
         "columns 'event', 'stream', 'location', 'time', 'loglikelihood'.")
  }
  if (any(is.na(loglikelihoods))) {
    stop("No support for missing values yet. Consider imputation methods.")
  }
  if (any(loglikelihoods[, loglikelihood > 0])) {
    stop("Log-likelihoods should be non-positive numbers. ",
         "Did you really take the logarithm?")
  }
  if (length(null_prior) != 1 || !is.numeric(null_prior) || null_prior < 0 ||
        null_prior > 1) {
    stop("Prior null hypothesis probability must be a single numeric value ",
         "between 0 and 1.")
  }
  if (any(event_priors <= 0)) {
    stop("Event priors must be positive numers, given either as probabilities ",
         "that together with the prior null hypothesis probability sum to 1, ",
         "or as numbers that express how often each event occurs compared ",
         "to the others.")
  }
  if (is.character(duration_condpriors) && 
        (length(duration_condpriors) !=1 ||
           duration_condpriors != "uniform")) {
    stop("Specify conditional probabilities for event durations given event ",
         "type either as the string 'uniform' or as a matrix.")
  }
  n1 <- length(unique(loglikelihoods[, event])) - 1
  n2 <- length(event_priors)
  n3 <- ifelse(!is.character(duration_condpriors), 
               ncol(duration_condpriors), n2)
  if (n1 != n2 || n1 != n3 || n2 != n3) {
    stop("The number of event types must be the same.")
  }
  if (!is.character(duration_condpriors) && 
      !isTRUE(all.equal(rep(1, length(event_priors)), 
                 unname(colSums(duration_condpriors))))) {
    stop("Conditional probabilities for event durations given event type ",
         "must sum to 1.")
  }
  # No errors; proceed with definitions-----------------------------------------

  events_are_strings <- typeof(loglikelihoods$event) == "character"
  if (events_are_strings) {
    event_names <- names(event_priors)
  }
  null_name <- ifelse(events_are_strings, "null", 0)
  
  # Normalize event priors to sum to 1 together with the null prior
  event_priors <- (1 - null_prior) * event_priors / sum(event_priors)
  
  times_durations <- times_and_durations(loglikelihoods)
  n_regions <- length(regions)
  n_events <- length(event_priors)
  max_duration <- max(times_durations[, duration])
  
  # If duration_condpriors is "uniform", create the appropriate matrix
  if (is.character(duration_condpriors)) {
    duration_condpriors <- matrix(1 / max_duration, 
                                  nrow = max_duration,
                                  ncol = n_events)
    if (events_are_strings) {
      duration_condpriors <- as.data.frame(duration_condpriors)
      names(duration_condpriors) <- event_names
    }
  }
  
  null_loglikelihood <- loglikelihoods[event == null_name, sum(loglikelihood)]
  
  if (n_partitions == "auto") {
    n_partitions <- min(length(regions), 
                        floor(log(sum(vapply(regions, length, integer(1))))))
  }
  
  region_partition <- partition_regions(regions, n_parts = n_partitions)

  # Set columns in correct order for adding log-likelihood ratios (add_llr)
  setkeyv(loglikelihoods, c("event", "location", "stream", "time"))
  
  # Calculate log-likelihood ratios for all events, regions, and durations
  keys_used_for_stllr <- c("region", "event", "duration", "stream", "location")
  spacetime_output <- 
    loglikelihoods %>%
    add_llr(null_name = null_name) %>%
    add_duration(keys = c("location", "stream", "duration", "event")) %>%
    region_apply(region_partition = region_partition, 
                 f = spacetime_llr, 
                 keys = keys_used_for_stllr)
  
  # Calculate the marginal probability of the data
  data_logratio <- 
    spacetime_output %>%
    spatial_llr(dur_given_event_logprobs = log(duration_condpriors)) %>%
    data_to_nulldata_logratio(event_logpriors = log(event_priors),
                              null_prior = null_prior,
                              n_regions = n_regions)
  marginal_logprob_of_data <- data_logratio + null_loglikelihood

  # Calculate posterior event log-probabilities for all space-time windows
  # Adds a column to 'spacetime_output' for posterior probability
  spacetime_logposterior(spacetime_output,
                         event_logpriors = log(event_priors),
                         dur_given_event_logprobs = log(duration_condpriors),
                         n_regions = n_regions,
                         data_logratio = data_logratio)
  
  # Calculate posterior event probabilities for all spatial windows
  spatial_posterior <- 1

  # Calculate probability maps
  event_logpmap <- 
    spacetime_output %>%
    event_logprobability_map(
      region_table = region_table_creator(regions, keys = "region"))

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
  setkey(spacetime_output, "posterior_prob")
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

#' Calculate the posterior event probability for each region and event type.
#' 
#' Calculate the posterior probability for each region and event type, 
#' marginalizing over the event durations.
#' @inheritParams k_most_probable
#' @return A \code{data.table} with columns \code{region, event} and 
#'    \code{posterior_prob}.
spatial_posteriors <- function(MBSS_obj) {
  MBSS_obj$posteriors[, .(posterior_prob = exp(logsumexp(posterior_logprob))),
                      by = .(region, event)]
}


#' Extract the posterior conditional event duration probabilities from an MBSS
#' object.
#' 
#' @param MBSS_obj An object of class "MBSS".
#' @param duration_column A logical scalar. Should the output matrix
#'    contain a column for the event duration?
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

#' Get the k most probable region-event-duration combinations.
#' 
#' Extracts the k most probable combinations of region, event type and event
#' duration from an MBSS object. One can also specify a restriction to a 
#' subset of the different event types, and a subset of the durations.
#' @param MBSS_obj An object of class "MBSS".
#' @param k An integer specifying how many of the most probable 
#'        region-event-duration combinations to output.
#' @param events Either a subset of all event types, or the string "all".
#'        Specifies which event types should be included in the output.
#' @param durations Either a subset of all event durations, or the string 
#'        "all". Specifies which event durations should be included in the 
#'        output.
k_most_probable <- function(MBSS_obj, 
                            k = 10, 
                            events = "all", 
                            durations = "all") {
  if (is.character(events) && events == "all") {
    events <- MBSS_obj$event_posteriors[, event]
  } 
  if (is.character(durations) && durations == "all") {
    durations <- seq(MBSS_obj$max_duration)
  }
  MBSS_obj$posteriors[event %in% events & duration %in% durations,
    .(region, event, duration, posterior_prob)][.N + 1 - seq(k)]
}