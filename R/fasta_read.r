#' Read sequences in FASTA format
#'
#' Reads sequences from a file in FASTA format. The function is applicable to
#' files with fixed-width formatting (sequences wrapped by newlines) as well as
#' to files without line breaks (single line per sequence).
#'
#' @param f File path as a character string.
#'
#' @return A named vector of type character, each element of which holds as
#'   sequence from the input file. Element names correspond to the sequence
#'    identifiers after removal of the initial ">" character.
#'
#' @note The function may fail on very large input files (e.g. metagenomes)
#'   not fitting into memory as a whole. Typical cases of file corruption
#'   are recognized and an error is thrown.
#'
#' @author David Kneis \email{david.kneis@@tu-dresden.de}
#'
#' @export
#'
#' @examples
#' # create sample file with and without line breaks
#' f <- tempfile()
#' write(c(">seq 1", "ATGCCTA", ">seq2", "AGATA", "CAT"), file=f)
#' 
#' # import sequences
#' s <- fasta_read(f)
#' 
#' rm(f)
#' print(s)

fasta_read <- function(f) {
  s <- readLines(f)
  block <- rep(0, length(s))
  for (i in 1:length(s)) {
    if (grepl(s[i], pattern="^>")) {
      if (i == 1) {
        block[i] <- 1
      } else if (i == length(s)) {
        stop(paste0("final line in file ",f," is a sequence header"))
      } else {
        block[i] <- block[i-1] + 1
      }
    } else {
      if (i == 1) {
        stop(paste0("initial line in file ",f," is not a sequence header"))
      } else {
        block[i] <- block[i-1]
      }
    }
  }
  keep <- !grepl(s, pattern="^>")
  if (any(diff(which(!keep)) == 1L))
    stop(paste0("sequence header immediately followed another sequence",
    " header in file ",f,""))
  ids <- s[!keep]
  s <- s[keep]
  block <- block[keep]
  s <- tapply(s, block, paste, collapse="")
  if (!identical(length(s), length(ids)))
    stop(paste0("number of sequence headers (",length(ids),") in file ",f,
      " does not match number of sequences (",length(s),")"))
  names(s) <- gsub(ids, pattern="^>", replacement="")
  s
}
