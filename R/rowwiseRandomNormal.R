#' @title Generate a data frame where each row follows a custom random normal distribution.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' 
#' @description Generate a data frame where each row contains random normal values
#' with the same mean and standard deviation as the corresponding row in the 
#' input data frame.
#'
#' @details
#' First version: 2015-03
#' Last modification: 2015-03
#'
#' @param x       A matrix or data frame
#'
#' @return
#' A data frame of same dimensions as the input matrix/data frame, 
#' with random normal values generated in a row-wise way.
#'
#' @examples
#' ## Load example data set from Den Boer, 2009
#' library(denboer2009)
#' data(denboer2009.expr)     ## Load expression table
#' data(denboer2009.pheno)    ## Load phenotypic data
#'
#' ## Define group labels and group of interest
#' g <- denboer2009.pheno$sample.labels
#' goi <- "Bh" ## Select one cancer subtype as group of interest
#'
#' ## Generate row-wise normal values
#' denboer2009.rnorm <- rowwiseRandomNormal(x=denboer2009.expr)
#'
#' ## Run Welch test on the row-wise permuted values
#' rnorm.welch <- meanEqualityTests(denboer2009.rnorm, 
#'                                  g=denboer2009.pheno$sample.labels, 
#'                                  goi="Bh",
#'                                  selected.tests="welch"
#'                                  )
#'                
#' ## Draw volcano plot of Welch test result with the random values.
#' ## This should show more or less no significant features.
#' meanEqualityTests.plotVolcano(rnorm.welch, test="welch", main="Random normal values, Welch volcano")
#'
#' @export
rowwiseRandomNormal <- function(x) {
  verbose(paste(sep="", "Generating row-wise random normal matrix (", 
                nrow(x), " rows x ", ncol(x), " columns)"), 1)
  return.rnorm.one.row <- function(x){ rnorm(n=length(x), mean=mean(x), sd=sd(x)) }
  return(t(apply(x, 1, return.rnorm.one.row)))
}

