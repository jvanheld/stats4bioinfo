#' @title Plot corrected p-values versus nominal p-value
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description
#' Plot the different multi-testing corrected statistics as a
#' function of the nominal P-value.
#' @param multitest.result the list returned by the function multipleTestingCorrections().
#' @param main='Multitesting corrections'  main title of the plot
#' @param alpha=0.05   Threshold of significance (alpha).
#' @param plot.pch=c(p.value=2,fdr=4,qval.0=3,e.value=1,fwer=20)  Specific characters to distinguish the plotted statistics.
#' @param plot.col=c(p.value='#000000',fdr='#888888',qval.0='#666666',e.value='#BBBBBB',fwer='#444444') Specific colors or gray levels to distinguish the plotted statistics.
#' @param plot.elements=c("p.value","e.value","fwer","fdr","qval.0") Selection of elements to display on the plot.
#' @param plot.legend=TRUE  Plot a legend indicating the number of features declared significant with the alpha threshold on the selected statistics.
#' @param legend.corner="topleft" corner wher the legend has to be placed.
#' @param legend.cex=1   Font size for the legend.
#' @param xlab="p-value" Label for the X axis
#' @param ylab="Multi-testing corrected statistics" Label for the Y axis
#' @param ...        Additional parameters are passed to plot()
#' @export
#' @examples
#' ## To obtain the input list (multitest.result), run the example of
#' ## stats4bioinfo::multipleTestingCorrections().
#' example(multipleTestingCorrections)
#' 
#' ## Plot all the multiple testing corrections at once
#' mulitpleTestingCorrections.plotCorrectedVsPval(multitest.result)
#'
#' ## Compare e-value and FWER
#' mulitpleTestingCorrections.plotCorrectedVsPval(multitest.result, plot.elements=c("e.value","fwer"))
#'
#' ## Compare e-value and FDR. 
#' ## This plot highlights the non-linear relationship between FDR and p-value.
#' mulitpleTestingCorrections.plotCorrectedVsPval(multitest.result, plot.elements=c("e.value","fdr"))
#'
#' ## Compare Benjamini-Hochberg (qval.0) and Storey-Tibshirani (fdr) estimates of FDR.
#' mulitpleTestingCorrections.plotCorrectedVsPval(multitest.result, plot.elements=c("qval.0","fdr"))
mulitpleTestingCorrections.plotCorrectedVsPval <- function (multitest.result,
                                                            main="Multitesting corrections",
                                                            xlab="p-value",
                                                            ylab="Multi-testing corrected statistics",
                                                            alpha=0.05,
                                                            legend.corner="topleft",
                                                            legend.cex=1,
                                                            plot.legend=TRUE,
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
  multitest.table <- multitest.result$multitest.table
  m <- multitest.result$m
  m0 <- multitest.result$m0
  Pi0 <- multitest.result$Pi0
  p.value.min = min(multitest.table$p.value)
  plot(multitest.table$p.value,
       multitest.table$p.value,
       main=main,
       xlab=xlab,
       ylab=ylab,
       log='xy',
       panel.first=c(
         abline(h=10^c(-10:3),col="#CCCCCC"),
         abline(v=10^c(-10:0),col="#CCCCCC"),
         abline(h=1,col="#666666")
       ),
       col=plot.col["p.value"], pch=plot.pch["p.value"],
       xlim=c(p.value.min, 1),
       ylim=c(p.value.min, m), type="n",
       ...)
  
  ## Plot the lines and prepare the legend
  legend.text <- vector()
  for (element in plot.elements) {
    lines(multitest.table$p.value, multitest.table[,element], col=plot.col[element], pch=plot.pch[element], type='p')
    n.below.alpha <- sum(as.vector(as.matrix(multitest.table[,element])) <= alpha)
    legend.text <- append(
      legend.text, 
      paste(sep="", element, " <= ", alpha, ": ", n.below.alpha)
    )
  }
  
  ## Show the alpha level
  abline(h=alpha, col='black', lty='dashed')
  
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
