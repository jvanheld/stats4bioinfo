#' @title Plot an histogram with a bootstrap test result.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Plot an histogram from the result of tTestPerRow.bootstrap(). 
#' The histogram indicates the number of features (ordinate) with a given support
#' (abcissa), i.e. the number of features declared significant in exactly 
#' X bootstrap iterations. 
#' 
#' @details
#' First version: 2015-04
#' Last modification: 2015-04
#' 
#' @param bootstrap.result    An object returned by the function tTestPerRow.bootstrap().
#' @param support.statistics="fdr"  Statistics for which the support has to be displayed.
#' Supported values: "fdr", "p.value", "e.value".
#' @param support.quantile=0.75 Minimal percent of support required to declare a feature positive.
#' The default is to select features supported in 75\% of bootstrap iterations.
#' @param xlab="Support" Label for the X axis.
#' @param ylab="Number of features" Label for the Y axis.
#' @param main=paste(sep="", "Bootstrap support (", support.statistics,"<=", bootstrap.result$alpha, ")")
#' @param plot.dcdf=FALSE  Plot the decreasing CDF over the histogram, in "step" mode. 
#' This curve indicates the number of features supported by at least X bootstraps.
#' @param plot.legend = TRUE  Only valid if plot.dcdf is TRUE. Plot a legend indicating 
#' the number of features supported by a given number of the iterations. 
#' @param legend.corner="topright"  Corner to plo the legend.
#' @param lwd=1  Line width for the dcdf and support threshold lines.
#' @param ... Additional parameters are passed to hist().
#' 
#' @examples
#' ## First run the examples of tTestPerRow.bootstrap(), in order to get some result to plot
#' example("tTestPerRow.bootstrap")
#' 
#' ## Plot the histogram
#' tTestPerRow.bootstrap.hist(student.bootstrap, col="#BBBBBB")
#' 
#' ## Plot the histogram + the dcdf. 
#' tTestPerRow.bootstrap.hist(student.bootstrap, col="#BBBBBB", 
#'                            plot.dcdf=TRUE, lwd=2)
#'
#' ## Plot the histogram + the dcdf for the e-value. 
#' tTestPerRow.bootstrap.hist(student.bootstrap, support.statistics="e.value",
#'                            col="#BBBBBB", 
#'                            plot.dcdf=TRUE, lwd=2)
#' 
#' @export
tTestPerRow.bootstrap.hist <- function(bootstrap.result, 
                                       xlab = "Support",
                                       ylab = "Number of features",
                                       support.statistics = "fdr",
                                       support.quantile=bootstrap.result$support.quantile,
                                       main = paste(sep="", "Bootstrap support (", support.statistics,"<=", bootstrap.result$alpha, ")"), 
                                       plot.dcdf=FALSE, 
                                       plot.legend = TRUE,
                                       legend.corner="topright",
                                       lwd=1,
                                       ylim=NULL,
                                       ...) {
  support.values <- as.vector(as.matrix(bootstrap.result$stats.per.row[paste(sep=".", support.statistics, "support")]))
  
  ## Compute dcdf if required, and adapt Y axis limits
  if (plot.dcdf) {
    if (is.null(ylim)) {
      ylim <- c(0, nrow(bootstrap.result$stats.per.row))
    }
  }
  
  h <- hist(support.values, 
            breaks=0:bootstrap.result$iterations,
            main=main,
            xlab=xlab, 
            ylab=ylab, 
            ylim=ylim, ...)
  
  if (plot.dcdf) {
    dcdf <- rev(cumsum(rev(h$counts)))
    min.support <- floor(bootstrap.result$iterations * support.quantile)
    supported.features <- dcdf[min.support]
    lines(h$mids-1/2,dcdf, type="s", col="blue", lwd=lwd)
    arrows(min.support, 0, min.support, supported.features, length=0, col="brown", lwd=lwd)
    arrows(min.support, supported.features, 0, supported.features, length=0, col="brown", lwd=lwd)
    if (plot.legend) {
      legend(legend.corner, 
             col=c("blue", "brown"),
             legend=c(
               "dcdf",
               paste(sep="", "N(support >= ",min.support, ") = ", sum(support.values>=min.support))
             ), 
             lwd=lwd)
    }
  }
}
