#' @title Build a position-specific scoring matrix (PSSM) from a set of sequences. 
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Given a vector of sequences, built a position-specific scoring matrix (PSSM)
#'  with different derived statistics (counts, frequencies, probabilities, weights, 
#' information content).
#'
#' @details
#' First version: 2016-12-23
#' Last modification: 2016-12
#'
#' @param sequences vector of strings corresponding to biological sequences (DNA, RNA, proteins)
#' @param prior=NULL vector of residue prior probabilities (names must correspond to residues)
#' @param pseudo.count=2 pseudo-count
#' @param IC.log.base=2 Logarithmic base for the information content
#' @param case.sensitive=FALSE by default residues are considered case-insensitive and converted to uppercases.
#'
#' @examples
#' ## Define the sequences of yeast Met31p binding sites
#' sequences <- c(
#'   "MET28"="cgcccAAAACTGTGGtgttag",
#'   "MET3"="gttgtAAAACTGTGGCTTTGT",
#'   "MUP3"="cggaaAAAACTGTGGcgtcgc",
#'   "SAM1"="acaggAAAACTGTGGtggcgc",
#'   "SAM2"="gcttgAAAACTGTGGcgtttt",
#'   "MET6"="gtcgcAAAACTGTGGtagtca",
#'   "MET30"="ccgcgCAAACTGTGGcttccc",
#'   "ZWF1"="ataagCAAACTGTGGgttcat",
#'   "MET14"="cctcaAAAAATGTGGcaatgg",
#'   "MET17"="tcatgAAAACTGTGTaacata",
#'   "MET2"="tgcaaAAAATTGTGGatgcac",
#'   "MET8"="ggaaaAAAAATGTGAaaatcg",
#'   "MET1"="cataaTAAACTGTGAacggac")
#' 
#' ## Chose priors based on yeast non-coding sequences
#' prior <- c("A"=0.32, "C"=0.18, "G"=0.18, "T"=0.32)
#' 
#' ## Build the PSSM
#' pssm <- seqToPSSM(seq=sequences, prior = prior)
#' 
#' ## Print count table
#' print(pssm$counts)
#' 
#' ## Print weight matrix
#' signif(pssm$weights, digits=2)
#' 
#' ## Plot a heatmap with the weights
#' heatmap.simple(pssm$counts, auto.margins=FALSE, xlab="Position", 
#'      ylab="Residues", main="Yeast Met13p count matrix", las=1)
#' 
#' @export
seqToPSSM <- function(sequences, 
                      prior = NULL,
                      pseudo.count = 2, # Pseudo-count
                      IC.log.base = 2, # Logarithmic base for the information content
                      case.sensitive = FALSE) {

  # Convert sequences to upper cases if required, 
  # to ensure case insensitivity
  if (!case.sensitive) {
    sequences <- toupper(sequences)
  }
  
  # Build a residue table from the sequence vector
  residue.table <- as.data.frame(strsplit(sequences, split = ""))
  
  # Define alphabet a prior residue probabilities
  residues <- sort(unique(as.vector(as.matrix(residue.table))))
  if (is.null(prior)) {
    alphabet <- residues
    prior <- rep(x = 1/length(alphabet), length.out=length(alphabet))
    names(prior) <- alphabet
  } else {
    alphabet <- names(prior)
    missing.residues <- setdiff(residues, alphabet)
    if (length(missing.residues) > 0) {
      message("Residues in sequence not defined in priors: ", setdiff(residues, alphabet))
      stop("Prior names should include all residue symbols.")
    }
  }
  
  
  ## Build a count matrix
  counts <- data.frame(matrix(ncol=nrow(residue.table), nrow=length(alphabet)))
  colnames(counts) <- 1:ncol(counts)
  rownames(counts) <- alphabet
  for (residue in alphabet) {
    counts[residue,] <- apply(residue.table==residue, 1, sum)
  }
  # apply(counts, 2, sum)
  
  ## Compute frequencies
  freq <- counts / apply(counts, 2, sum)
  # apply(freq, 2, sum)
  
  ## Compute probabilities, as frequencies smoothed by a pseudo-count
  counts.corrected <- counts + pseudo.count * prior
  proba <- counts.corrected / apply(counts.corrected, 2, sum)
  
  ## Compute weights, i.e. log-ratios between position-specific probabilities and priors
  weights <- log(proba / prior)
  
  ## Compute information content
  info <- weights * proba
  info.colsum <- apply(info, 2, sum)
  
  # apply(counts.corrected, 2, sum)
  result <- list(
    alphabet = alphabet,
    prior = prior,
    nb.sites = max(apply(counts, 2, sum)),
    sequences = sequences,
    IC.log.base = IC.log.base,
    counts = counts,
    freq = freq,
    counts.corrected = counts.corrected,
    proba = proba,
    weights = weights,
    info = info,
    info.colsum = info.colsum
  )
  return(result)
}
