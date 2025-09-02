#' Calculate the number of target copies from a uniplex ddPCR output
#'
#' Given the number of positive droplets and the total number of droplets,
#' the function computes an estimate of the number of copies in the sample.
#' If the percentage of positives is very low such that all droplets contain
#' a single copy of the target only, this is trivial. However, at a higher
#' percentage of positives, we need to account for the fact that some droplets
#' will contain more than one copy. Hence, we need to use the Poisson
#' distribution.
#'
#' @param n_positive Number of positive droplets among \code{n_total} (integer,
#'   possibly a vector representing multiple measurements).
#' @param n_droplets Total number of droplets (integer vector of the same length
#'   as \code{n_positive}).
#'
#' @return The estimated number of target copies (numeric).
#'
#' @note The plain solution is
#'
#'    \deqn{copies = -1 * ln(1 - n_positive / n_droplets) * n_droplets}
#' 
#' with \eqn{ln} being the natural logarithm. See comments in the code for the
#' derivation.
#'
#' @author David Kneis \email{david.kneis@@tu-dresden.de}
#'
#' @export
#'
#' @examples
#' n_droplets <- 1000
#' n_positive <- seq(from=1, to=200, by=10)
#' copies <- sapply(n_positive, ddPCR_uniplex, n_droplets=n_droplets)
#' plot(n_positive, copies, xlab="Positive droplets", ylab="Copies in sample")
#' legend("bottomright", bty="n", lty=c(1,NA), pch=c(NA,1),
#'   legend=c("1:1 relation","estimates based Poisson distribution"))
#' abline(0, 1)

ddPCR_uniplex <- function(
  n_positive,        # number of droplets being positive for the target
  n_droplets         # total number of droplets
) {
  if (!all(sapply(c(n_positive, n_droplets), is.numeric))) {
    stop("droplet counts must be numeric")
  }
  if (!identical(length(n_positive), length(n_droplets))) {
    stop("input vectors must be of equal length")
  }
  if (any(n_positive > n_droplets)) {
    stop("number of positive droplets cannot exceed total number of droplets") 
  }
  # observed empirical probability
  p_positive_obs <- n_positive / n_droplets
  # find parameter Lambda of the Poisson distribution such that the expected
  # probability of positives (p_positive_exp) equals the observed empirical
  # probability p_positive_obs
  # notes:
  # * we use the principle: p_positive = 1 - p_negative
  # * p_negative_exp can be found using ppois(q=0, lambda) but we can directly
  #   use the formula of the CDF for k=0 which is p_negative = exp(-lambda)
  #   (see formula of the CDF on wikipedia)
  # * to find the lambda, we simply set p_positive_obs == p_positive_exp and
  #   solve for lambda
  lambda <- -1 * log(1 - p_positive_obs)
  # lambda represents the most probably true average of copies per droplet
  copies <- lambda * n_droplets
  round(copies)
}
