#' Break a sequence into random fragments
#'
#' The function simulates the shearing of, e.g., a DNA molecule. The molecule
#' is broken at random locations into a set of fragments.
#'
#' @param x The original sequence as a character string (scalar).
#' @param avglen The expected average fragment length (mean value). Should
#'   be a single number greater than zero.
#' @param circular Logical. If \code{TRUE}, the input sequence is assumed to
#'   represent a circular molecule like, e.g., a fully intact plasmid and
#'   the ends of \code{x} are joined before fragmentation starts.
#'   If \code{FALSE}, a linear molecule with given front and tail is assumed.
#'
#' @return A vector of type character, each element of which is a random
#'   sub-section of the input character string passed in argument \code{x}.
#'
#' @note
#' The actual number of generated fragments is based on Poisson-distributed
#' random numbers derived from the provided value of \code{avglen}.
#'
#' @author David Kneis \email{david.kneis@@tu-dresden.de}
#'
#' @export
#'
#' @examples
#'
#' # to create sample sequences
#'
#' makeSeq <- function(len) {
#'   bases <- c("A","C","G","T")
#'   paste(bases[sample.int(length(bases), size=len, replace=TRUE)],
#'     collapse="")
#' }
 
#' # shear a single sequence
#'
#' s <- makeSeq(50)
#' print(s)
#' frag <- seq_fragment(s, nchar(s)/3, FALSE)
#' print(frag)
#' 
#' # shear multiple sequences, applying identical settings for each
#' 
#' s <- c(s1=makeSeq(50), s2=makeSeq(30), s3=makeSeq(60))
#' print(s)
#' frag <- do.call(c, sapply(s, seq_fragment,
#'   avglen=mean(sapply(s, nchar))/3, simplify=FALSE))
#' print(frag)
#' 
#' # analyze fragment length distribution
#' 
#' s <- makeSeq(64000)   # original sequence (typical plasmid size)
#' avglen <- 1e4         # average fragment length desired
#' frag <- c()           # fragment 100 replicates
#' for (i in 1:100) {
#'   frag <- c(frag, seq_fragment(s, avglen=avglen, TRUE))
#' }
#' hist(log10(sapply(frag, nchar)), xaxt="n", breaks=20,
#'   main="Fragment length histogram", xlab="bp (log scaled)")
#' axis(side=1, at=1:6, labels=10^(1:6))
#' abline(h=avglen, col="red")
#' legend("topleft", bty="n", legend=c(
#'   paste0("Average length requested: ",avglen),
#'   paste0("Average length achieved: ",round(mean(sapply(frag, nchar))))
#' ))

seq_fragment <- function(x, avglen, circular=TRUE) {
  if (!identical(length(x), 1L) || !is.character(x))
    stop("expecting 'x' to be a scalar character string")
  if (!identical(length(avglen), 1L) || avglen <= 0)
    stop("expecting 'avglen' to be a scalar number > 0")
  if (!identical(length(circular), 1L) || !is.logical(circular))
    stop("expecting 'circular' to be a scalar logical value")
  # utility to account for the absence of natural ends in circular molecules
  break_and_rejoin <- function(x) {
    at <- sample.int(nchar(x), size=1)
    if (at > 1) {
      x <- paste0(substr(x, at, nchar(x)), substr(x, 1, at-1))
    } else {
      x 
    }
  }
  # get the expected number of breaks
  nbr <- stats::rpois(1, nchar(x) / avglen) - 1
  # return original sequence in case of no breaks, build fragments otherwise
  if (nbr <= 0) {                
    x
  } else {
    x <- if (circular) break_and_rejoin(x) else x  # ensure random start point
    # convert to vector
    x <- unlist(strsplit(x, split=""))
    # perform fragmentation
    cuts <- sample.int(n=length(x), size=nbr, replace=TRUE) # set cut positions
    cuts <- sort(unique(c(0, cuts, length(x))))             # register ends too
    len <- diff(cuts)                            # length of fragments
    i <- rep(1:length(len), len)                 # vector of fragment indices
    tapply(x, i, paste, collapse="")             # build fragments
  }
}

