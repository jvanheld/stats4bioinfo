#' @title Compute the parameters of a Markov model (transition matrix) from k-mer occurrences
#' 
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' 
#' @description Compute a transition matrix from a list of k-mer occurrences. 
#' 
#' @details
#' First version: 2003-09
#' Last modification: 2015-02
#'
#' @param kmer.seq A vector of k-mer sequences. 
#' @param kmer.occ A vector of the same length as seq.
#' @param alphabet Exhaustive list of residues to take in consideration as suffixes for the transition matrix.
#' 
#' @return A list with the parameters of the model (k, m) + the transition matrix.
#' @examples 
#' 
#' ## Compute Markov model from dinucleotide occurrences measures in Escherichia coli genome
#' org <- "Escherichia_coli_K12"
#' k <- 3
#' freq.file.path <- file.path("markov_models", org, paste(sep="", k, "nt_genomic_",org,"-1str-ovlp.freq.gz"))
#' freq.file <- system.file("extdata", freq.file.path, package = "stats4bioinfo")
#' message("Loading k-mer frequency file", freq.file)
#' freq.table <- read.delim(freq.file, row.names = 1, comment.char=";")
#' head(freq.table)
#' tail(freq.table)
#' markov.model <- kmerFrequenciesToMarkov(
#'     kmer.seq = freq.table$id, 
#'     kmer.occ = freq.table$occ)
#' transitions <- markov.model[["transition"]]
#' heatmap.simple(round(digits=3, as.matrix(transitions)),main=org, zlim=c(0,1))
#' 
#' @export
kmerFrequenciesToMarkov <- function(kmer.seq, 
                                    kmer.occ,
                                    alphabet = c("A", "C", "G", "T")) {
  
  kmer.seq <- as.vector(kmer.seq)
  kmer.occ <- as.vector(kmer.occ)
  k <- nchar(kmer.seq[1]) ## Length of the k-mers
  m <- k - 1 ## Markov order
  
  # define a function for the sum, because for some genomes sums exceed max integer value    
  sumnum <- function(x) {
    sum(as.numeric(x))
  }
  
  kmer.table <- data.frame(row.names=kmer.seq, Occurrences=kmer.occ)
  kmer.table$Frequency <- kmer.table$Occurrences / sumnum(kmer.table$Occurrences)
  
  ## Compute a transition matrix (order m=k-1) from k-mer frequencies
  count.matrix <- data.frame(matrix(nrow=0, ncol=4))
  names(count.matrix) <- toupper(alphabet)
  i <- 1
  for (i in 1:nrow(kmer.table)) {
    seq <- toupper(kmer.seq[i])
    occ <- kmer.occ[i]
    if (m == 0) {
      prefix <- "."
    } else {
      prefix <- substring(text=seq, first = 1, last = nchar(seq)-1)
    }
    suffix <- substring(text=seq, first = nchar(seq))
    count.matrix[prefix, paste(sep="", suffix)] <- occ
  }
  
  count.matrix["Total", ] <- apply(count.matrix[, alphabet], 2, sumnum)
  count.matrix$Total <- apply(count.matrix[, alphabet], 1, sumnum)
  
  ## Derive a frequency matrix from the count matrix
  transition.matrix <- count.matrix[,alphabet] / count.matrix[, "Total"]
  transition.matrix$Prefix <- 2*count.matrix[, "Total"]/sumnum(count.matrix[,"Total"])
  row.names(transition.matrix) <- sub(row.names(transition.matrix), pattern = "Total", replacement = "Suffix")
  transition.matrix["Suffix", "Prefix"] <- NA

  model <- list(m=m,k=k, counts = count.matrix, transition=transition.matrix)
  return(model)
}
