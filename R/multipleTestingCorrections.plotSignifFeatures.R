#' @title Plot number of significant features as a function of the significance threshold.
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description
#' Plot the number of significant features as a function of the control
#' criterion (nominal p-value, e-value, fdr, ...).
#' @param multitest.result the list returned by the function multipleTestingCorrections().
#' @param main="Significant features"  main title of the plot
#' @param xlab="P-value derived statistics"
#' @param ylab="Significant features"
#' @param alpha=0.05   Threshold of significance (alpha).
#' @param plot.legend=TRUE  Plot a legend indicating the number of features declared significant with the alpha threshold on the selected statistics.
#' @param legend.corner="topleft" corner wher the legend has to be placed.
#' @param legend.cex=1   Font size for the legend.
#' @param plot.pch=c(p.value=2,e.value=1,fwer=20,fdr=4,qval.0=3)  Specific characters to distinguish the plotted statistics.
#' @param plot.col=c(p.value='#000000',e.value='#BBBBBB',fwer='#444444',fdr='#888888',qval.0='#666666') Specific colors or gray levels to distinguish the plotted statistics.
#' @param plot.elements=c("p.value","fdr","qval.0","e.value","fwer") Selection of elements to display on the plot.
#' @param ...        Additional parameters are passed to plot()
#'  
#' @examples
#' ## To obtain the input list (multitest.result), see the documentatio of
#' example(multipleTestingCorrections)
#'
#' ## Plot all the multiple testing corrections at once
#' mulitpleTestingCorrections.plotSignifFeatures(multitest.result)
#'
#' ## Compare e-value and FWER
#' mulitpleTestingCorrections.plotSignifFeatures(multitest.result, plot.elements=c("e.value","fwer"))
#'
#' ## Compare e-value and FDR
#' mulitpleTestingCorrections.plotSignifFeatures(multitest.result, plot.elements=c("fdr","e.value"))
#'
#' ## Compare Benjamini-Hochberg (qval.0) and Storey-Tibshirani (fdr) estimates of FDR
#' mulitpleTestingCorrections.plotSignifFeatures(multitest.result, plot.elements=c("fdr","qval.0"))
#' @export
mulitpleTestingCorrections.plotSignifFeatures <- function (multitest.result,
                                                           main="Significant features",
                                                           xlab="Significance threshold",
                                                           ylab="Significant features",
                                                           alpha=0.05,
                                                           plot.legend = TRUE,
                                                           legend.corner="topleft",
                                                           legend.cex=1,
                                                           xlim=NULL,
                                                           plot.pch=c(p.value=2,
                                                                      fdr=4,
                                                                      qval.0=3,
                                                                      e.value=1,
                                                                      fwer=20),
                                                           plot.col=c(p.value='#000000',
                                                                      fdr='#888888',
                                                                      qval.0='#666666',
                                                                      e.value='#BBBBBB',
                                                                      fwer='#444444'),
                                                           plot.elements=c("p.value",
                                                                           "fdr",
                                                                           "qval.0",
                                                                           "e.value",
                                                                           "fwer"),
                                                           ...
) {
  
  ## Retrieve the multitesting result table
  multitest.table <- multitest.result$multitest.table
  
  ## Define X limits
  p.value.min <- min(multitest.table$p.value) ## Minimal p-value
  m <- multitest.result$m ## Number of tests
  if (is.null(xlim)) {
    xlim <- c(p.value.min, m)
  }
  
  ## Draw the plot frame
  plot(x=multitest.table$p.value,
       y=multitest.table$rank,
       main=main,
       xlab=xlab,
       ylab=ylab,
       log='xy',
       panel.first=c(
         grid(equilogs=F, lty='solid', col='#CCCCCC'),
         abline(h=10^c(-11:3), col='#CCCCCC'),
         abline(h=1,col="#666666")
       ),
       col=plot.col["p.value"], pch=plot.pch["p.value"],
       xlim=xlim, type="n",
       ylim=c(1,nrow(multitest.table)),
       ...)
  
  ## Plot the lines and prepare the legend
  legend.text <- vector() ## Collect legend text
  n.below.alpha.values <- vector() ## Collect number of values below the alpha threshold for each stat
  for (element in plot.elements) {
    lines(x=multitest.table[,element], 
          y=multitest.table$rank, 
          col=plot.col[element], pch=plot.pch[element], type='p')
    n.below.alpha <- sum(as.vector(as.matrix(multitest.table[,element])) <= alpha)   
    n.below.alpha.values <- append(n.below.alpha.values, n.below.alpha)
    legend.text <- append(
      legend.text, 
      paste(sep="", element, " <= ", alpha, ": ", n.below.alpha)
    )
    #    segments(n.below.alpha, alpha, n.below.alpha, p.value.min, lty="dashed")
  }
  
  ## Plot landmark lines
  abline(v=alpha, col='black', lty='dashed')
  abline(h=n.below.alpha.values, lty="dashed")
  
  
  ## Plot the legend
  if (plot.legend) {
    legend(legend.corner, 
           legend=legend.text,
           bg='#FFFFFF', bty="o", 
           cex=legend.cex,
           col=plot.col[plot.elements],
           pch=plot.pch[plot.elements],
    )
  }
}
