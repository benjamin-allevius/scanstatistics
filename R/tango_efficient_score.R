




score_numerator <- Vectorize(function(x, W) sum(W:1 * x[1:W]), 
                       vectorize.args = "W")

score_denominator <- Vectorize(function(x, W) sum((W:1)^2 * x[1:W]), 
                         vectorize.args = "W")


# factors in numerator and denominator of eqn(21) Tango 2011 
# if using negative binomial distribution
efficient_score_terms_nbin <- function(table) {
  table[, .(num = sum((count - baseline) / overdispersion),
        denom = sum(baseline / overdispersion)),
        by = .(region, duration)]
}

# numerator and denominator of eqn(21) Tango 2011 if using binomial dist
efficient_score_terms_binom <- function(table) {
  table[,
    .(num = sum((count - baseline)),
      denom = sum(baseline)),
    by = .(region, duration)]
}

# components of eqn (21) Tango 2011
outbreak_efficient_score <- function(table) {
  table[,
    .(duration = duration,
      score = score_numerator(num, duration)
      / sqrt(score_denominator(denom, duration))), 
    by = .(region)]
}

# components of eqn (23) Tango 2011
hotspot_efficient_score <- function(table) {
  table[,
    .(duration = duration,
      score = cumsum(num) / sqrt(cumsum(denom))),
    by = .(region)]
}