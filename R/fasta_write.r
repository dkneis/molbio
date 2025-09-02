#' Write sequences in FASTA format
#'
#' Exports a vector of sequences in FASTA format.
#'
#' @param x A named vector of type character holding sequence information.
#' @param f File path as a character string.
#'
#' @return NULL
#'
#' @note Sequences are written to the file without line breaks. Sequence headers
#'   are constructed from element names by appending the usual ">" as a suffix.
#'
#' @author David Kneis \email{david.kneis@@tu-dresden.de}
#'
#' @export
#'
#' @examples
#' # export sequence
#' f <- tempfile()
#' fasta_write(c(`seq 1`="ATG", `seq2`="CAT"), f)
#' 
#' # re-import
#' s <- fasta_read(f)
#' 
#' rm(f)
#' print(s)

fasta_write <- function(x, f) {
  for (i in 1:length(x)) {
    write(file=f, x=paste0(">",names(x)[i]), append=i>1) 
    write(file=f, x=x[i], append=T) 
  }
  invisible(NULL)
}
