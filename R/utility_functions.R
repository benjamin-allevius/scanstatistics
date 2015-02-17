

#' log-sum-exp trick
#' 
#' @param
#' @return
logsumexp <- function(x) {
    A <- max(x)
    A + log(sum(exp(x - A)))
}