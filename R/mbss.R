

MBSS <- function(event_loglikelihood_ratios,
                 null_loglikelihood,
                 regions,
                 null_prior,
                 event_priors,
                 duration_condpriors) {
  # Input checking -------------------------------------------------------------
  if (!is.data.table(loglikelihood_ratios)) {
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
  if (length(null_prior) != 1 || !is.numeric(null_prior) || null_prior >= 0 ||
        null_prior <= 1) {
    stop("Prior null hypothesis probability must be a single numeric value ",
         "between 0 and 1.")
  }
  if (!all.equal(1, null_prior + sum(event_priors))) {
    stop("Prior probabilities for events or no events must sum to 1.")
  }
  if (!all.equal(rep(1, length(event_priors)), 
                 rowSums(duration_condpriors))) {
    stop("Conditional probabilities for event durations given event type ",
         "must sum to the event priors.")
  }
  
  # 
  n_regions <- length(regions)
  max_duration <- length(unique(event_loglikelihood_ratios[, time]))
  
#   if (any(grepl("POSIX", class(event_loglikelihood_ratios$time)))) {
#     actual_times <- 
#   }
  
  # Assume everything is correct below
  
  
  spacetime_output <- 
    event_loglikelihood_ratios %>%
    region_joiner(regions = regions) %>%
    spacetime_llr
  
  data_logratio <- 
    spacetime_output %>%
    spatial_llr(dur_given_event_logprobs = log(duration_condpriors)) %>%
    data_to_nulldata_logratio(event_logpriors = log(event_priors),
                              null_prior = null_prior,
                              n_regions = n_regions)
  
  marginal_prob_of_data <- exp(data_logratio + null_loglikelihood)

  # Modify spacetime_output: add colum for posterior probability
  spacetime_logposterior(spacetime_output,
                         event_logpriors = log(event_priors),
                         dur_given_event_logprobs = log(duration_condpriors),
                         n_regions = n_regions,
                         data_logratio = data_logratio)

  setkeyv(spacetime_output, c("region", "event"))

  event_logpmap <- 
    spacetime_output %>%
    event_logprobability_map(
      region_table = region_table_creator(regions, key = "region"))
  
  event_pmap <- event_probability_map(event_logpmap)
  pmap <- probability_map(event_logpmap)
  
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
                 event_probability_map = event_pmap,
                 probability_map = pmap), 
            class = "MBSS")
}