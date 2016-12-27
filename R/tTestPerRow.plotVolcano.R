#' @title Volcano plot from a tTestPerRow result object.
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Draw a volcano plot from the table produced by stats4bioinfo::tTestPerRow().
#' Abcissa represents the mean difference between groups.
#' Ordinate the significance of the Welch test (sig=-log10(e-value)).
#' @param ttpr.result Must be the result obtained with tTestPerRow()
#' @param control.type="fdr"   Type of control (supported: "fdr", "e.value", "p.value")
#' @param alpha=0.05    Alpha threshold for the control of false positives
#' @param plot.ci=FALSE if TRUE, draw confidence intervals for each feature of the volcano plot.
#' @param ... Additional parameters are passed to VolcanoPlot()
#' @return no return object
#' @examples
#' ## Parameters to generate and analyse the simulated dataset
#' n.h1 <- 500 ## Number of rows under H1
#' n.h0 <- 500 ## Number of rows under H0
#' alpha <- 0.05
#' sample.labels <- c(rep("a", 50), rep("b", 50))
#' 
#' ## Generate an artificial dataset with n.h1 truly positive and n.h0 truly negative rows
#' samples.per.groups <- unlist(table(sample.labels))
#' h1.data <- rnormPerGroup(n=samples.per.groups, mean=c(0, 0.5), sd=c(1,1), nrow=n.h1)
#' h0.data <- rnormPerGroup(n=samples.per.groups, mean=c(0, 0), sd=c(1,1), nrow=n.h0)
#' x <- rbind(h1.data$x, h0.data$x)
#' x.status <- rep(c("H1","H0"), times=c(n.h1, n.h0))
#' table(x.status)
#' 
#' ## Run Student t-test on the simulated data set
#' x.student <- tTestPerRow(x, var.equal=TRUE, cl=sample.labels, alpha=alpha, test.group="b")
#' 
#' ## Draw a classical Volcano plot with -log10(p-value) on the Y axis. 
#' tTestPerRow.plotVolcano(x.student, legend.corner="topleft", control.type = "p.value")
#' 
#' ## Add a threshold on the effect size
#' tTestPerRow.plotVolcano(x.student, legend.corner="topleft", control.type = "p.value", effect.threshold = 0.6)
#' 
#' ## Draw volcano plot with status- and density-based colors
#' tTestPerRow.plotVolcano(x.student, legend.corner="topleft", control.type = "p.value", effect.threshold = 0.6, density.colors=TRUE)
#' 
#' ## Draw a Volcano plot with horizontal bars denoting the confidence intervals 
#' ## around the difference between means
#' tTestPerRow.plotVolcano(x.student, legend.corner="topleft", 
#'    plot.ci=TRUE, control.type="p.value", alpha=)
#' 
#' ## Draw a volcano plot where the colors represent the true status (TP, FP, TN, FN) 
#' ## rather than the significance.
#' (confusion.table <- table(x.status, positive=x.student$table$p.value < alpha))
#' pred.status <- rep(NA, times=nrow(x.student$table))
#' pred.status[x.status=="H1" & x.student$table$p.value < alpha] <- "TP" ## True positives
#' pred.status[x.status=="H0" & x.student$table$p.value >= alpha] <- "TN" ## True negatives
#' pred.status[x.status=="H0" & x.student$table$p.value < alpha] <- "FP" ## False positive
#' pred.status[x.status=="H1" & x.student$table$p.value >= alpha] <- "FN" ## False negatives
#' table(pred.status)
#' pred.status.colors <- c("TP"="darkgreen", "FP"="red", "FN"="orange", "TN"="grey")
#' tTestPerRow.plotVolcano(x.student, legend.corner=NULL, 
#'    plot.ci=TRUE, plot.points=TRUE, tick.size=0.08, control.type="p.value", 
#'    col.points=pred.status.colors[pred.status])
#'  legend("topleft", legend=paste(table(pred.status), names(table(pred.status))),
#'         col=pred.status.colors[names(table(pred.status))], lwd=2)
#' 
#' ## Draw two volcano plots separating positive and negative tests, to show that 
#' ## the p-value cutoff is equivalent to a selection of the confidence interval 
#' ## that do not cross the 0 value. Set the p-value threshold to 1/N, which is 
#' ## equivalent to set the control on e-value <= 1.
#' par(mfrow=c(1,2))
#' corrected.alpha <- 1/nrow(x.student$table)
#' tTestPerRow.plotVolcano(x.student, col.points=NULL, col.positive="black",
#'    legend.corner=NULL, plot.ci=TRUE, tick.size=0, 
#'    control.type="p.value", alpha=corrected.alpha, main="Positive tests")
#' tTestPerRow.plotVolcano(x.student, col.points="grey", col.positive=NA,
#'    legend.corner=NULL, plot.ci=TRUE, tick.size=0, 
#'    control.type="p.value", alpha=corrected.alpha, main="Negative tests")
#' par(mfrow=c(1,1))
#' 
#' @export
tTestPerRow.plotVolcano <- function(
  ttpr.result,
  control.type="fdr",
  alpha=ttpr.result$alpha,
  ylab=paste(sep="", "-log10(", control.type, ")"),
  plot.ci=FALSE,
  ... ## additional parameters are passed to the VolcanoPlot function
) {
  
  multi.t <- ttpr.result$table
  
  ## Identify the positive tests
  if (alpha != ttpr.result$alpha) {
    warning(paste("Selecting positive with alpha=", alpha, "different from ttpr alpha=", ttpr.result$alpha))
    if (plot.ci) {
      warning(paste(sep="", "Computing ",100*(1-alpha),"% confidence intervals"))
      ## Recompute the confidence interval width using user-provided alpha
      ## Compute the width of a confidence interval around the difference between means
      if (ttpr.result$alternative == "two.sided") {
        critical.t <- qt(p=1-alpha/2, df=multi.t$df)
      } else {
        critical.t <- qt(p=1-alpha, df=multi.t$df)
      }
      multi.t[, "ci.width"] <- 2*critical.t * multi.t$diff.stder
    }
  }
  
  VolcanoPlot(multitest.table=multi.t, control.type=control.type, alpha=alpha, effect.size.col="means.diff", ylab=ylab, plot.ci=plot.ci, ...)
}


