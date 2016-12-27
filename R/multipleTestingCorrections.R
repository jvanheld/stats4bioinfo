#' @title Multiple corrections for multiple testing
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description
#' Apply various multiple testing corrections on a list of P-values, report
#' all of them in a table with one row per test and one column per correction,
#' and optionally draw illutrative figures.
#'
#' In particular, we apply the elegant strategy from Storey & Tibshirani (2003)
#' to estimate the respective proportions of true and alternative hypotheses
#' among the tests.
#'
#'
#' @details
#' First version: 2012-02-10
#' Last modification: 2015-02
#'
#' @param p.values             list of p-values
#' @param feature.names=NULL       Names associated to the p.value list (will be used as row.names for to result table)
#' @param lambda=0.5       lambda parameter for estimating pi0 = m0/m
#' @param alpha=0.05   Threshold of significance (alpha).
#' @param main='Multitesting corrections'
#' Prefix for the main title of the plots
#' @param run.qvalue=FALSE
#' For the sake of comparison, compute FDR and generate graphs
#' using Storey's qvalue library
#' (\url{http://genomics.princeton.edu/storeylab/qvalue/}).
#' @param plots=FALSE   If true, generate illutrsative plots
#' @param plot.pch=c(p.value=2,
#'                   fdr=4,
#'                   qval.0=3,
#'                   e.value=1,
#'                   fwer=20)   point type (character) associated to each statistics for the plot
#' @param plot.col=c(p.value='black',
#'                   fdr='blue',
#'                   qval.0='darkviolet',
#'                   e.value='darkbrown',
#'                   fwer='orange')
#' color for the plot
#' @param plot.elements=c("p.value", "fdr", "qval.0","e.value", "fwer")
#' Elements to draw on the plots (can be any subet of the default list).
#' @return
#' A list comprizing the following elements:
#'  \item{lambda}{Value provided for the input parameter lambda.}
#'  \item{alpha}{Value provided for the input parameter threshold (alpha).}
#'  \item{m}{total number of tests (= number of elements in the input list of p-values). m = m0 + m1.}
#'  \item{m0.est}{Estimation of the number of cases under null hypothesis.}
#'  \item{m1.est}{Estimation of the number of cases that do not comply with the null hypothesis (m1.est = m - m0.est).}
#'  \item{pi0}{Estimated proportion of trully null cases among all tests. pi0 = m0.est / m. }
#'  \item{nb.signif}{Number of tests declared significant above the specified alpha threshold}
#'  \item{multitest.table}{
#'  A table with one row per test, and one column per statistics (the fields described below).}
#'  \item{p.value}{P-value = probability to observe an at least as significant result under the null hypothesis. This column contains a replica of the input vector of p-values.}
#'  \item{rank}{rank by increasing P-value}
#'  \item{e.value}{E-value = expected number of false positives}
#'  \item{fdr}{False Discovery Rate = expected proportion of truly null hypotheses among the cases declared positives. The FDR is estimated using Storey and Tibshirani procedure.}
#'  \item{fwer}{Family-Wise Error Rate = probability to observe at least one FP among all the tests, assuming all of them are under null hypothesis.}
#'  \item{qval.0}{FDR assuming that all the tests are under null hypothesis (Benjamini-Hochberg approach).}
#' @references
#' 1. Benjamini, Y. and Hochberg, Y. (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. JOURNAL-ROYAL STATISTICAL SOCIETY SERIES B.
#' 2. Storey, J.D. and Tibshirani, R. (2003) Statistical significance for genomewide studies. Proc Natl Acad Sci USA, 100, 9440–9445.
#' @examples
#' ## Generate a vector of 10,000 fake p-values
#' ## with a given proportion of trully null (m0=9,000) and
#' ## non-null (m1=1,000) cases.
#' my.p.values <- c(runif(n=9000, min=0, max=1),
#'                  runif(n=1000, min=0, max=1e-2))
#'
#' multitest.result <- multipleTestingCorrections(my.p.values)
#' attributes(multitest.result)
#'
#' ## In principle, m0 should be ~9,000
#' print(paste("m0.est =", multitest.result$m0.est))
#'
#' ## In principle, m1 should be ~1,000
#' print(paste("m1.est =", multitest.result$m1.est))
#'
#' ## In principle, pi0 should be ~0.9
#' print(paste("pi0 =", multitest.result$pi0))
#'
#' ## Plot P-value distributions
#' mulitpleTestingCorrections.plotPvalDistrib(multitest.result)
#'
#' ## Compare the different multiple testing corrections (Y axis)
#' ## versus the nominal p-value (X axis).
#' mulitpleTestingCorrections.plotCorrectedVsPval(multitest.result)
#'
#' ## Plot the number of significant features for different
#' ## multiple testing corrections
#' mulitpleTestingCorrections.plotSignifFeatures(multitest.result)
#' 
#' @export
multipleTestingCorrections <- function(p.values,
                                       feature.names=NULL,
                                       lambda=0.5,
                                       alpha=0.05,
                                       main='Multitesting corrections',
#                                       file.prefix=NULL,
                                       run.qvalue=FALSE,
                                       plots=FALSE
                                       ) {

  m <- length(p.values) ## Total number of tests (m = m0 + m1)

  ## Create a table for storing the results
  multitest.table <- data.frame(p.value=p.values)
  if (!is.null(feature.names)) {
    rownames(multitest.table) <- feature.names
  } else if (!is.null(names(p.values))) {
    rownames(multitest.table) <- names(p.values)
  }
  multitest.table$rank <- rank(p.values)

  ## E-value
  multitest.table$e.value <- p.values * m

  ## Estimate FDR using the R function p.adjust
  multitest.table$fdr <- p.adjust(p.values, method="fdr")

  ## Family-wise error rate (fwer)
  multitest.table$fwer <- pbinom(q=0, size=m, prob=multitest.table$p.value, lower.tail=FALSE)
  #  multitest.table$fwer <- 1 - ( 1 - multitest.table$p.value)^m

  ## Compute q-value with Storey's library
  ##
  ## Beware: Storey's qvalue() function assumes that p-values are sorted
  ## by increasing order. However, for multilpeTestingCorrections we don't
  ## want to impose this constraint, because we want to be able to annotate
  ## Any user-provided table without changing the order ot hte rows.
  ## We fix this by sorting p-values before sending them to quvalue(),
  ## and re-ordering them by rank of the original p-value vector before
  ## adding them to the multitesting.result table.
  if (run.qvalue) {
    library(qvalue)
    qobj <- qvalue::qvalue(sort(multitest.table$p.value), lambda=lambda)
    multitest.table$qval.Storey <-  qobj$qvalues[rank(p.values)]
    #     if (plots) {
    #       ## Draw plots of the function qvalue::qplot
    #       qvalue::qplot(qobj)
    #       if (!is.null(file.prefix)) {
    #         dev.copy2pdf(paste(file.prefix, "_qvalue_plots.pdf"),
    #                      width=8,height=6);
    #       }
    #     }
  }
  
  ## q-value_0 (q-value corresponding to the case where all the features
  ## would be trully null).
  multitest.table$qval.0 <- multitest.table$p.value * m/multitest.table$rank
  #   par(mfrow=c(1,2))
  #   plot(multitest.table$qval.0[order(multitest.table$p.value)], log="xy", panel.first=grid(), main="qval.0 before cummax")
  multitest.table$qval.0[order(multitest.table$p.value)] <- cummax(multitest.table$qval.0[order(multitest.table$p.value)])
  #   plot(multitest.table$qval.0[order(multitest.table$p.value)], log="xy", panel.first=grid(), main="qval.0 after cummax")
  #   par(mfrow=c(1,1))
  
  ## Estimate numbers of null and alternative hypotheses
  m0.est <- min(m, sum(p.values >= lambda) / (1-lambda)) ## m0 cannot be larger than m
  m1.est <- m - m0.est
  pi0 <- m0.est/m

  #  print(c(lambda, m0.est))


  ## Prepare the result
  result <- list()
  if (run.qvalue) {
    result$qobj <- qobj
  }
  result$multitest.table <- multitest.table
  result$lambda <- lambda
  result$alpha <- alpha
  result$m <- m
  result$m0.est <- m0.est
  result$m1.est <- m1.est
  result$pi0 <- pi0
  result$nb.signif <- c(
    p.value=sum(multitest.table$p.value <= alpha),
    bonferoni=sum(multitest.table$e.value <= alpha),
    e.value.le.1=sum(multitest.table$e.value <= 1),
    fwer=sum(multitest.table$fwer <= alpha),
    qval.0=sum(multitest.table$qval.0 <= alpha),
    fdr=sum(multitest.table$fdr <= alpha)
  )
  if (run.qvalue) {
    result$nb.signif["qval.Storey"] = sum(multitest.table$qval.Storey <= alpha)
  }


  ################################################################
  ## Draw plots
  if (plots) {
    mulitpleTestingCorrections.plotPvalDistrib(multitest.result)
#    mulitpleTestingCorrections.plotCorrectedVsPval(multitest.result)
#    mulitpleTestingCorrections.plotSignifFeatures(multitest.result)
  }

  return(result)
}


