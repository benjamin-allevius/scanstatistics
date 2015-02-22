

# log-sum-exp trick
logsumexp <- function(x) {
    A <- max(x)
    A + log(sum(exp(x - A)))
}