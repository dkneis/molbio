#' Analyze a duplex ddPCR output regarding linked and unlinked targets
#'
#' The function analyzes the output of a ddPCR run containing targets 'A' and 'B',
#' but also 'ab'. Here, 'ab' represents cases where A and B are present on the same
#' fragment of DNA. The tricky part is to properly resolve the double-positive
#' droplets into those containing the physically linked target (ab) and those
#' containing the two individual targets (A and B) without physical linkage.
#' Likewise, when estimating the number of copies of a single target, e.g. A,
#' we need to account for mixed droplets (e.g., A and B, A and ab) which are
#' equally observed as double-positive. Finally, there may be multiple copies
#' of the same target in a single droplet if the proportion of empty droplets
#' is low (which is generally undesired but may happen in real-worlds settings).
#'
#' @param x A named numeric vector with elements 'A', 'B', 'AandB', 'N' or a
#' matrix with these column names. If a matrix, each row represents a measurement.
#' In the input,
#' 'A' is the number of droplets being detected as positive for the first target
#' only. Likewise, 'B' is the number of droplets being positive
#' for the second target only. 'AandB' is the observed number of
#' droplets being double-positive for any reason, i.e. because they contain
#' 'ab' or 'A' and 'B', including mixtures thereof like, e.g., 'ab' and 'A'.
#' Finally, 'N' represents the total number of droplets.
#'
#' @return A numeric vector or matrix, depending on the input, with (column)
#'   names 'A', 'B', 'ab'. Values represent estimates of
#'   the respective numbers of copies in the sample. Here, 'ab' refers to
#'   copies of physically linked targets while 'A' and 'B' refer to copies
#'   of the physically isolated targets. Note that the output values of 'A' and
#'   'B' will generally be larger than the corresponding input values
#'   due to the occurrence of multiple targets per droplet. For example, target
#'   A is contained in droplets having A only, but also in those having A and B,
#'   A and ab, A and B and ab, and there may even be cases where a droplet holds
#'   multiple copies of a particular target (e.g. A and A, A and A and B, ...).
#'
#' @note
#' The solution is currently computed numerically although an analytic solution
#' of the underlying system of equations is possible. In fact, the
#' analytic solution can be found in Eq. 4 of this supplement
#' \url{https://doi.org/10.1371/journal.pone.0118270.s005} where they compute
#' directly the average number of linked targets per droplet, i.e. lambda,
#' the parameter of the Poisson distribution. To get the number of linked copies,
#' one just needs to multiply that lambda with the total number of droplets.
#'
#' @author David Kneis \email{david.kneis@@tu-dresden.de}
#'
#' @export
#'
#' @examples
#' 
#' # a simple call
#' ddPCR_linked_duplex(c(N=1000, A=30, B=50, AandB=100))
#' 
#' # run a series of virtual ddPCR experiments and compare the estimated
#' # copy numbers to the known true numbers
#' out <- NULL
#' for (replicate_number in 1:30) {
#'   
#' # random numbers for current experiment
#'   n_droplets <- 1000
#'   copies_A <- runif(n=1, min=5, max=500)
#'   copies_B <- runif(n=1, min=5, max=500)
#'   copies_ab <- runif(n=1, min=5, max=500)
#'   
#'   # generate droplets containing 0, 1, or more copies of the targets
#'   A <- rpois(n=n_droplets, lambda=copies_A/n_droplets)
#'   B <- rpois(n=n_droplets, lambda=copies_B/n_droplets)
#'   ab <- rpois(n=n_droplets, lambda=copies_ab/n_droplets)
#'   
#'   # observations corresponding to droplets
#'   pos <- function(x) {x > 0}
#'   neg <- function(x) {x == 0}
#'   obs <- c(
#'     A=sum(pos(A) & neg(B) & neg(ab)),
#'     B=sum(pos(B) & neg(A) & neg(ab)),
#'     AandB=sum(pos(ab) | (pos(A) & pos(B)))
#'   )
#'   
#'   # estimate the number of copies
#'   est <- ddPCR_linked_duplex(c(N=n_droplets, obs))
#'   
#'   # save estimates together with known values
#'   out <- rbind(out, data.frame(
#'     A.true=copies_A, B.true=copies_B, ab.true=copies_ab, as.list(est)
#'   ))
#' }
#' 
#' # plot estimates vs. truth
#' par(mfrow=c(3,1))
#' for (x in c("A","B","ab")) {
#'   plot(out[,paste0(x,".true")], out[,x],
#'     xlab="true number of copies", ylab="estimated number of copies")
#'   abline(0,1)
#'   legend("topleft", bty="n", legend=paste("Target:",x))
#' }
#' par(mfrow=c(1,1))

ddPCR_linked_duplex <- function(x) {
  
  # workhorse function for evaluating a single experiment
  eval_single <- function(x) {
    obs <- x[c("A", "B", "AandB")]
    N <- x["N"]
    rm(x)
    # use the observed numbers as initial guess; this will be close to truth if
    # the proportion of empty droplets is large
    guess <- c(A=obs[["A"]], B=obs[["B"]], ab=obs[["AandB"]])
    # solve the characteristic system of equations originating from set theory
    #           +------------------+
    #   +-------|--------+         |
    #   |   (A) | A&B    |   (B)   |
    #   |       |        |         |
    #   |   +---|--------|------+  |
    #   |   |   | A&B&ab | ab&B |  |
    #   |   |   +--------|------|--+
    #   |   | A&ab       |      |
    #   +---|------------+      |
    #       |        (ab)       |
    #       +-------------------+
    # Considering three intersecting circles, each representing one of the targets
    # (A, B, ab) one can derive the following equations. The left hand side
    # represents the observable signals, other symbols represent the actual numbers
    # of droplets containing the targets. For example "AB" means "A" and "B" and
    # "Aab" means "A and ab". Note that both would be contained in "AandB.obs"
    # because they are double positive.
    #
    # A.obs     = A        - AB - Aab       + abAB
    # B.obs     =    B     - AB       - Bab + abAB
    # AandB.obs =       ab + AB             - abAB
    #
    # Using simple theory of combined probabilities (p(a & b) = p(a) * p(b)),
    # the following holds, with N being the total number of droplets.
    #
    # AB = A * B / N
    # Aab = A * ab / N
    # Bab = B * ab / N
    # abAB = A * B * ab / N^2
    
    objf <- function(p, obs) {
      with(as.list(p),
        sum(c(
          A      - A*B/N - A*ab/N          + A*B*ab/(N^2) - obs["A"],
            B    - A*B/N          - B*ab/N + A*B*ab/(N^2) - obs["B"],
              ab + A*B/N                   - A*B*ab/(N^2) - obs["AandB"]
        )^2)
      )
    }
    opt <- stats::optim(par=guess, fn=objf, obs=obs)
    if (opt$convergence != 0) {
      stop("failed to solve the characteristic system of equations")
    }
    est <- stats::setNames(opt$par, names(opt$par))
    # convert positive droplets to copies using the Poisson distribution;
    # --> could also call 'ddPCR_uniplex' here; see there for theory
    est <- -1 * log(1 - est / N) * N
    round(est)
  }
  
  # apply the workhorse function after checking inputs
  
  # always use a matrix as input, even if a vector
  if (is.vector(x)) {
    x <- t(as.matrix(x))
  }
  if (is.null(colnames(x))) {
    stop("input must have names (if vector) or column names (if matrix)")
  }
  if (!identical(sort(colnames(x)), c("A","AandB","B","N"))) {
    stop("elements of input must be named 'A', 'B', 'AandB', 'N'")
  }
  if (!all(is.numeric(x))) {
    stop("input must be numeric")
  }
  out <- NULL
  for (i in 1:nrow(x)) {
    out <- rbind(out, eval_single(x[i,]))
  }
  if (nrow(out) == 1) {
    out <- out[1,]
  }
  out
}
