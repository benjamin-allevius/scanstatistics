




num_score <- Vectorize(function(x, W) sum(W:1 * x[1:W]), 
                       vectorize.args = "W")

denom_score <- Vectorize(function(x, W) sum((W:1)^2 * x[1:W]), 
                         vectorize.args = "W")


# numerator and denominator of eqn(21) Tango 2011 if using NBin dist
efficient_score_terms_nbin <- function(d) {
  d[, .(num = sum((count - baseline) / overdispersion),
        denom = sum(baseline / overdispersion)),
    by = .(region, time)]
}

# numerator and denominator of eqn(21) Tango 2011 if using binomial dist
efficient_score_terms_binom <- function(d) {
  d[, .(num = sum((count - baseline)),
        denom = sum(baseline)),
    by = .(region, time)]
}

# components of eqn (21) Tango 2011
outbreak_efficient_score <- function(d) {
  d[,
    .(time = time,
      score = num_score(num, time + 1)
      / sqrt(denom_score(denom, time + 1))), 
    by = .(region)]
}

# components of eqn (23) Tango 2011
hotspot_efficient_score <- function(d) {
  d[,
    .(time = time,
      score = cumsum(num) / sqrt(cumsum(denom))),
    by = .(region)]
}