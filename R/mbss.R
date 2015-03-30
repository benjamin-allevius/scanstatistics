# TODO: documentation
# Change input to densities instead of likelihood ratios
# Work out how to convert back to actual times 
# - pass in observed times and counts?

MBSS <- function(counts_and_logdensities,
                 regions,
                 null_prior,
                 event_priors,
                 duration_condpriors = "uniform") {
  # Input checking -------------------------------------------------------------
  if (!is.data.table(counts_and_logdensities)) {
    stop("The log-likelihood ratios must be supplied as a data.table.")
  }
  if (!all(c("event", "stream", "location", "time", "count", 
             "null_logdensity", "event_logdensity") %in% 
             names(counts_and_logdensities))) {
    stop("The data.table containing the log-likelihood ratios must ",
         "contain the columns 'event', 'stream', 'location', 'time', 'count', ",
         "'null_logdensity', 'event_logdensity'.")
  }
  if (any(is.na(counts_and_logdensities))) {
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
  # Proceed with calculations --------------------------------------------------
  times_durations <- times_and_durations(counts_and_logdensities)
  n_regions <- length(regions)
  n_events <- length(event_priors)
  max_duration <- max(times_durations[, duration])
  
  
  if (is.character(duration_condpriors)) {
    duration_condpriors <- matrix(1 / max_duration, 
                                  nrow = max_duration,
                                  ncol = n_events)
  }
  
  null_loglikelihood <- 
    counts_and_logdensities[event == counts_and_logdensities[1, event],
                            sum(null_logdensity)]

  # Calculate log-likelihood ratios for all events, regions, and durations
  keys_used_for_stllr <- c("region", "event", "duration", "stream", "location")
  spacetime_output <- 
    counts_and_logdensities[, .(llr = event_logdensity - null_logdensity),
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
  
  structure(list(counts_and_logdensities = counts_and_logdensities,
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