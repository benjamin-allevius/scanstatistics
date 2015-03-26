# TODO: documentation
# Change input to densities instead of likelihood ratios
# Work out how to convert back to actual times 
# - pass in observed times and counts?

#   if (any(grepl("POSIX", class(event_loglikelihood_ratios$time)))) {
#     actual_times <- 
#   }

MBSS <- function(event_loglikelihood_ratios,
                 null_loglikelihood,
                 regions,
                 null_prior,
                 event_priors,
                 duration_condpriors) {
  # Input checking -------------------------------------------------------------
  if (!is.data.table(event_loglikelihood_ratios)) {
    stop("The log-likelihood ratios must be supplied as a data.table.")
  }
  if (!all(names(event_loglikelihood_ratios) %in% 
            c("event", "stream", "location", "time", "llr"))) {
    stop("The data.table containing the log-likelihood ratios must ",
         "contain the columns 'stream', 'location', 'time', 'llr'.")
  }
  if (any(is.na(event_loglikelihood_ratios))) {
    stop("No support for missing values yet.")
  }
  if (length(null_prior) != 1 || !is.numeric(null_prior) || null_prior < 0 ||
        null_prior > 1) {
    stop("Prior null hypothesis probability must be a single numeric value ",
         "between 0 and 1.")
  }
  if (!all.equal(1, null_prior + sum(event_priors))) {
    stop("Prior probabilities for events or no events must sum to 1.")
  }
  if (!all.equal(rep(1, length(event_priors)), 
                 colSums(duration_condpriors))) {
    stop("Conditional probabilities for event durations given event type ",
         "must sum to the event priors.")
  }
  # Proceed with calculations --------------------------------------------------
  n_regions <- length(regions)
  max_duration <- length(unique(event_loglikelihood_ratios[, time]))
  
  setkeyv(event_loglikelihood_ratios, c("location"))
  
  # Calculate log-likelihood ratios for all events, regions, and durations
  spacetime_output <- 
    event_loglikelihood_ratios %>%
    region_joiner(regions = regions) %>%
    spacetime_llr
  
  # Calculate the marginal probability of the data
  data_logratio <- 
    spacetime_output %>%
    spatial_llr(dur_given_event_logprobs = log(duration_condpriors)) %>%
    data_to_nulldata_logratio(event_logpriors = log(event_priors),
                              null_prior = null_prior,
                              n_regions = n_regions)
  marginal_prob_of_data <- exp(data_logratio + null_loglikelihood)

  # Add a column to \code{spacetime_output} for posterior probability
  spacetime_logposterior(spacetime_output,
                         event_logpriors = log(event_priors),
                         dur_given_event_logprobs = log(duration_condpriors),
                         n_regions = n_regions,
                         data_logratio = data_logratio)

  setkeyv(spacetime_output, c("region", "event"))

  # Calculate probability maps
  event_logpmap <- 
    spacetime_output %>%
    event_logprobability_map(
      region_table = region_table_creator(regions, key = "region"))

  event_pmap <- event_probability_map(event_logpmap)
  pmap <- probability_map(event_logpmap)

  # Calculate the posterior probabilities
  spacetime_output[, posterior_prob := exp(posterior_logprob)]
  null_posterior <- exp(null_loglikelihood) * null_prior / marginal_prob_of_data
  event_duration_joint <- posterior_duration_event_jdist(spacetime_output)
  event_posteriors <- posterior_event_probabilities(event_duration_joint)
  duration_condposteriors <- posterior_duration_givn_event(event_duration_joint)

  if (!all.equal(1, null_posterior + 
                    event_posteriors[, sum(event_posterior)])) {
    warning("Posterior probabilities don't sum to 1.")
  }
  
  structure(list(event_loglikelihood_ratios = event_loglikelihood_ratios,
                 max_duration = max_duration,
                 n_regions = n_regions,
                 regions = regions,
                 null_loglikelihood = null_loglikelihood,
                 null_prior = null_prior,
                 event_priors = event_priors,
                 duration_condpriors = duration_condpriors,
                 marginal_prob_of_data = marginal_prob_of_data,
                 posteriors = spacetime_output,
                 null_posterior = null_posterior,
                 event_posteriors = event_posteriors,
                 duration_condposteriors = duration_condposteriors,
                 event_probability_map = event_pmap,
                 probability_map = pmap), 
            class = "MBSS")
}