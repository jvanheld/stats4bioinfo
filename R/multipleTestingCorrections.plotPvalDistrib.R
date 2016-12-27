#' @title Multiple corrections for multiple testing
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description
#' Generate a plot that reproduces Fig 1 from Storey and Tibshirani (2003),
#' with some additional details, in order to illustrate the estimation of
#' the parameter PI0 = m0/m1.
#' @param multitest.result the list returned by the function multipleTestingCorrections().
#' @param main='Multitesting corrections'  main title of the plot
#' @param plot.legend=TRUE Plot a legend with some indicative numbers (m0, m1, pi0).
#' @param legend.corner="topright" corner wher the legend has to be placed.
#' @param legend.cex=1   Font size for the legend.
#' @param breaks=seq(from=0,to=1, by=0.05)    Breaks for the histogram
#' @param draw.lambda="arrow" Indicate the level of the lambda parameter. Supported 
#' representations: "arrow", "line", "none". 
#' @param draw.m0.line=TRUE Draw an horizontal line indicating the estimated m0 
#' (number of trully null features) per bin.
#' @param draw.mean.line=FALSE Draw a dashed horizontal line indicating the mean 
#' number of features per bin. The difference between this line and the "m0 per bin" 
#' line reflects the importance of truly alternative features.
#' @param col="#BBBBFF" Histogram background color (passed to hist()).
#' @param overlay=NULL  Boolean vector marking features of a special category 
#' (e.g. truly null features). A colored histogram will be printed on top of the 
#' main histogram, to indicate the number of features belonging to this group.
#' @param overlay.col="#CCCCCC" Color for the overlay histogram.
#' @param mean.line.col="black" Color to draw the line indicating the mean number of feature per bin.
#' @param m0.line.col="black" Color to draw the line indicating the estimated m0 per bin.
#' @param ...   Additional parameters are passed to hist()
#' @export
#' @examples
#' ## To obtain the input list (multitest.result), run the examples of
#' ## stats4bioinfo::multipleTestingCorrections().
#'
#' example(multipleTestingCorrections)
#'
#' ## Plot the p-value distribution + landmarks
#' mulitpleTestingCorrections.plotPvalDistrib(multitest.result, draw.lambda="line")
#'
mulitpleTestingCorrections.plotPvalDistrib <- function(multitest.result,
                                                       main='P-value distribution',
                                                       plot.legend=TRUE,
                                                       legend.cex=1,
                                                       legend.corner="topright",
                                                       breaks=seq(from=0,to=1, by=0.05),
                                                       draw.lambda = "arrow",
                                                       draw.m0.line=TRUE,
                                                       draw.mean.line=TRUE,
                                                       overlay=NULL,
                                                       col="#FFEEDD",
                                                       overlay.col="#CCCCCC",
                                                       mean.line.col="darkred",
                                                       m0.line.col="blue",
                                                       ...
) {
  
  ## Plot the histogram
  distrib <- hist(multitest.result$multitest.table$p.value,
                  breaks=breaks, plot=TRUE,
                  xlab="nominal p-value",
                  main=main, col=col,
                  ...)
  m.per.bin <- sum(distrib$counts)/length(distrib$counts)
  m0.per.bin <- multitest.result$m0.est / length(distrib$count)
  
  ## Plot the histogram of overlay values, if specified
  if (!is.null(overlay)) {
    distrib.overlay <-  hist(multitest.result$multitest.table$p.value[overlay],
                             breaks=breaks, add=TRUE,col=overlay.col)
  }
  
  ## Draw horizontal line to show the expected frequency per histogram bar,
  ## and the frequency based on m0 estimation
  if (draw.mean.line) {
    abline(h=m.per.bin, col=mean.line.col, lty='dashed', lwd=2)
  }
  if (draw.m0.line) {
    abline(h=m0.per.bin, col=m0.line.col, lty='solid', lwd=2)
  }
  
  ## Draw an arrow to highlight the lambda parameter
  if (draw.lambda == "arrow") {
    arrows(multitest.result$lambda, m.per.bin*1.5, 
           multitest.result$lambda, m.per.bin*1.1, 
           angle = 30, length=0.05, code = 2, col=m0.line.col, lwd = 2)
  } else if (draw.lambda == "line") {
    abline(v=multitest.result$lambda, lty="dotted", col=m0.line.col, lwd=2)
  }
  
  ## Add a legend
  if (plot.legend) {
    legend(legend.corner, bty='o', bg='white', border=NA,
           cex=legend.cex,
           legend=(c(paste("lambda =", multitest.result$lambda),
                     paste("m =", multitest.result$m),
                     paste("per bin =", m.per.bin),
                     paste("m0.est =", round(multitest.result$m0.est)),
                     paste("m1.est =", length(multitest.result$multitest.table$p.value) - round(multitest.result$m0.est)),
                     paste("pi0 =", round(multitest.result$pi0, digits=3)))
           ))
  }
}


