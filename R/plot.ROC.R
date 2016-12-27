#' @title Plot a ROC curve
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Plot a ROC curve from a table returned by ROCstatistics. 
#' @details
#' First version: 2015-04
#' Last modification: 2015-04
#' @param ROC.table, ## A data fame resulting from the function ROCstatistics()
#' @param main='ROC' main title for the plot
#' @param xlab='negative scores'
#' @param ylab='positive scores'
#' @param add=FALSE If true, the ROC curve is added on a pre-existing plot
#' @param ... Additional parameters are passed to lines().
#' @export
plot.ROC <- function(ROC.table, ## A data fame resulting from the function ROCstatistics()
                               main='ROC',
                               xlab='negative scores',
                               ylab='positive scores',
                               add=FALSE, ## If true, the ROC curve is added on a pre-existing plot
                               ...
) {
  if (!add) { 
    plotROCframe(main=main,xlab=xlab,ylab=ylab) 
  }
  x <- c(0, sort(ROC.table$neg.iCDF), 1)
  y <- c(0, sort(ROC.table$pos.iCDF), 1)
  lines(x,y,...)
}

################################################################
## Old function names maintained for backward compatiility
plot.ROC.curves <- function (...) {
  plotROCfromFiles(...)
}

plot.ROCstatistics <- function (...) {
  plot.ROC(...)
}
