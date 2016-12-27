#' @title Draw a Volcano plot.
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Draw a volcano plot from a table containing at least one column for the effect size, and another one for the p-value or an equivalent measure of significance (FDR, e-value, FWER, ...).
#' @param result.table A data frame containing one row per feature, and one column per statistics.
#' @param effect.size.col A column number or name, inidicating which column of the result table contains the effect size, which will be 
#' @param control.type="p.value"  A column number or name, indicating which column of the result table contains the p-value or an equivalent indication of the significance of each feature (example: "fdr", "e.value", "p.value")
#' @param nb.tests=nrow(multitest.table) Total number of tests performed (by default, the number of rows in the table). 
#' @param alpha=0.05    Alpha threshold for the control of false positives
#' @param effect.threshold=NULL Threshold on the absolute value of the effect size.
#' @param sort.by.pval=FALSE Sort row by p-value in order to plot significant elements on top of non-significant
#' @param xlab="Effect size" Label for the X axis.
#' @param ylab="sig = -log10(p-value)" Label for the Y axis
#' @param xlim    Range of the X axis.
#' @param ylim    Range of the Y axis.
#' @param density.colors=FALSE Automatically set the colors according to feature status and local density in the Volcano space.
#' @param col.points='#888888' Color(s) for the points. can be either a single value (same color for all points), or a vector of the same length as the numer of genes (rows) of the input table. 
#' @param col.positive='#4455DD' Color to highlight significant points (called positive). When NULL, positive points are not displayed. 
#' @param col.lines='blue'  Color for the line highlighting the effect size thresholds.
#' @param col.alpha="darkred" Color for the line highlighting the significance threshold.
#' @param col.grid='#AAAAAA'   Grid color
#' @param lty.lines="solid" Line type for the lines.
#' @param lty.alpha="dashed" Line type for the horizontal line denoting the alpha theshold.
#' @param lty.grid="dotted" Line type for the grid.
#' @param pch.positive=3   Point shape for the positive tests (default: 3, i.e. + symbol)
#' @param pch.negative=1   Point shape for the negative tests (default: 1, i.e. o symbol)
#' @param legend.corner="bottomleft"  Corner for the legend. When NULL, no legend is displayed.
#' @param full.legend=TRUE Plot additional indications on the legend (number of elements passing the different adjusted p-values)
#' @param legend.cex=1   Font size for the legend.
#' @param cex=0.6     Point size, passed to plot() and lines()
#' @param plot.points=TRUE Plot one point per row with the significance (sig=-log10(p-value)) as a function of the effect size.
#' @param plot.ci=FALSE Plot horizontal bars to denote confidence intervals around the mean difference.
#' @param tick.size=0.05 Height of the vertical lines denoting the boundaries of confidence intervals.
#' @param ... Additional parameters are passed to plot()
#' @return no return object
#' 
#' @export
VolcanoPlot <- function(
  multitest.table,
  effect.size.col,
  control.type="p.value",
  nb.tests=nrow(multitest.table),
  alpha=0.05,
  effect.threshold=NULL,
  sort.by.pval=FALSE,
  plot.points=TRUE,
  plot.ci=FALSE,
  xlab="Effect size",
  ylab=paste(sep="", "-log10(",control.type,")"),
  density.colors=FALSE,
  col.points = '#888888',
  col.positive='#4455DD',
  col.grid = '#AAAAAA',
  col.lines = 'blue',
  col.alpha = "darkred",
  lty.lines = "solid",
  lty.alpha = "dashed",
  lty.grid="dotted",
  pch.positive=3,
  pch.negative=1,
  xlim = NULL,
  ylim = NULL,
  cex=0.6,
  legend.corner="bottomleft",
  full.legend=FALSE,
  legend.cex=1,
  tick.size=0.05,
  #                         Y.score = "sig", ## Score to plot on the Y axis. Supported: sig (default), p.value e.value
  ... ## additional parameters are passed to the plot function
) {
  
  if (is.null(col.positive)) {
    col.positive = NA
  }
  
  ## Identify features declared positive according to the specified alpha
  positive <- multitest.table[,control.type] <= alpha
  # table(positive)
  
  ## Apply threshold on effect size if required
  if (!is.null(effect.threshold)) {
    positive[abs(multitest.table[,effect.size.col]) < effect.threshold] <- FALSE
  }
  
  ## Sort the table before plotting, to ensure that positive points appear visible on top of negative points
  if (sort.by.pval) {
    order <- order(multitest.table[,control.type], decreasing=TRUE)
    multitest.table.sorted <- multitest.table[order,]
    if (length(col.points) == nrow(multitest.table)) {
      col.points <- col.points[order]
    }
  } else {
    multitest.table.sorted <- multitest.table
  }
  
  ## Select Y values depending on the control type (p.value, e.value, fdr)
  y.values.ori <- -log10(multitest.table.sorted[, control.type])
  
  ## Fix a problem with infinite Y values resulting from 0 values for the control (p-value or derived stat)
  y.values <- y.values.ori
  y.value.max <- 320 ## Maximal value for Y corresponds to the precision of floating point computation for the p-values (~ 1e-320)
  y.values[is.infinite(y.values)] <- y.value.max
  
  ## Replace NA values by -1, so they appear outside of the plot
  y.values[is.na(y.values)] <- -1
  
  ## Compute confidence interval limits
  if (plot.ci) {
    ci.min <- multitest.table.sorted[, effect.size.col] - multitest.table.sorted[, "ci.width"]/2
    ci.max <- multitest.table.sorted[, effect.size.col] + multitest.table.sorted[, "ci.width"]/2
  }  
  
  ## Define limits of X and Y scales
  if (is.null(xlim)) {
    if (plot.ci) {
      xmax <- max(abs(c(ci.min, ci.max)), na.rm=TRUE)
    } else {
      xmax <- max(abs(multitest.table[, effect.size.col]), na.rm=TRUE)
    }
    xlim<- c(-xmax, xmax)*1.05 ## Add a 5% margin on X axis
  }
  if (is.null(ylim)) {
    ylim <- c(min(y.values), max(1,y.values))
  }
  
  ## Define point colors and shapes. 
  ## Note: the attributes col.points and pch.points can either be either 
  ## a single value (for all points) or a  vector with one user-specified 
  ## color per point.
  if (density.colors) {
    ## Compute status and density-specific color
    multitest.table.sorted$status <- "negative"
    multitest.table.sorted$status[positive] <- "positive"
    multitest.table.sorted$color <- featureColors(
      multitest.table.sorted$status, 
      positions = data.frame(multitest.table.sorted[, effect.size.col],y.values))
  }  else {
    multitest.table.sorted$color <- col.points
    if (length(col.points) != nrow(multitest.table.sorted)) {
      if ((is.na(col.positive)) | (is.null(col.positive))) {
        multitest.table.sorted[positive, "color"] <- rep(NA, times=sum(positive))
      } else {
        multitest.table.sorted[positive, "color"] <- col.positive
      }
    }
  }
  
  ## Point character for each feature
  multitest.table.sorted$pch <- pch.negative
  if (!is.null(pch.positive)) {
    multitest.table.sorted[positive, "pch"] <- pch.positive
  }
  
  ## Identify the points above Y limits, to denote them by a different symbol and color
  above.ymax <- y.values.ori > ylim[2]
  y.values[above.ymax] <- ylim[2]
  multitest.table.sorted$pch[above.ymax] <- 17
  multitest.table.sorted$color[above.ymax] <- "purple"
  
  ################################################################
  ## Draw the Volcano plot
  plot(multitest.table.sorted[, effect.size.col],y.values,
       xlab=xlab,
       ylab=ylab,
       cex=cex,
       xlim=xlim,
       ylim=ylim,
       type="n",
       panel.first=grid(lty=lty.grid, col=col.grid),
       ...)
  
  ## Plot the extent + boundaries confidence intervals
  if (plot.ci) {
    arrows(x0 = ci.min, 
           y0 = y.values,
           x1 = ci.max,
           y1 = y.values, 
           col=multitest.table.sorted$color, angle=90, length=tick.size/2, code=3)
  } 
  
  ## Plot the points
  if (plot.points) {
    points(x=multitest.table.sorted[, effect.size.col],
           y=y.values, cex=cex,
           col=multitest.table.sorted$color,
           pch=multitest.table.sorted$pch)
  }
  
  
  ## Vertical line to denote the null position (corresponding to no difference between groups)
  abline(v=0,col="black", lty="dashed", lwd=1) 
  
  ## Add vertical lines to denote the threshold on effect size
  if (!is.null(effect.threshold)) {
    abline(v=c(-effect.threshold, effect.threshold),col=col.lines, lwd=1, lty=lty.lines) ## Vertical line to denote the threshold on effect size
  }
  
  ## Add horizontal lines to denote alpha level
  abline(h=-log10(alpha),col=col.alpha, lwd=1, lty=lty.alpha) ## Horizontal line to denote the significance threshold
  
  ## Plot the legend
  if (!is.null(legend.corner)) {
    ## Store legend prameters in a table
    legend.table <- data.frame()
    
    ## Legend for alpha threshold
    if (!is.null(alpha)) {
      legend.table <- rbind(legend.table, 
                            data.frame("name"="alpha", 
                                       "legend"=paste(sep="", control.type," < ", signif(digits=4, alpha)),
                                       "lwd"=2, "col"=col.alpha, "pch"=-1, "lty"=lty.alpha))
    }
    
    
    ## Legend for threshold on effect size
    if(!is.null(effect.threshold)) {
      legend.table <- rbind(legend.table, 
                            data.frame("name"="effect", 
                                       "legend"=paste(sep="", "abs(", effect.size.col,") >= ", signif(digits=4, effect.threshold)),
                                       "lwd"=2, "col"=col.lines, "pch"=-1, "lty"=lty.lines))
    }
    
    ## Legend for the points, depending on their status + display options
    nb.positives <- sum(positive, na.rm = TRUE)
    nb.negatives <- nb.tests - nb.positives
    legend.table <- rbind(legend.table, 
                          data.frame("name"="positive", 
                                     "legend"=paste(sep="", nb.positives, " positives"), 
                                     "lwd"=1, "col"= col.positive[1], "pch"=-1, "lty"="blank"))
    legend.table <- rbind(legend.table, 
                          data.frame("name"="negative", 
                                     "legend"=paste(sep="", nb.negatives, " negatives"), 
                                     "lwd"=1, "col"= col.points[1], "pch"=-1, "lty"="blank"))
    
    row.names(legend.table) <- legend.table[, "name"]
    if (plot.points) {
      legend.table$pch <- -1
      legend.table["positive", "pch"] <- pch.positive
      legend.table["negative", "pch"] <- pch.negative
    }
    
    
    ## Line type for legend elements
    if (plot.ci) {
      legend.table$lty <- "solid"
      legend.table["positive", "lty"] <- "solid"
      legend.table["negative", "lty"] <- "solid"
    }
    
    ## Plot a legend with additional info on the number of positives for different multiple testing corrections
    if (full.legend) {
      legend.pch <- c(as.vector(legend.table$pch), -1,-1,-1,-1)
      legend(legend.corner, 
             c(paste(sep="","N=",nrow(multitest.table)),
               paste(sep="", sum(positive), " positives (",control.type," <= ", alpha, ")"),
               paste(sum(!positive), "negatives"),
               paste(sep="", "P-value <= ", alpha, ": ", sum(multitest.table[,"p.value"] <= alpha)),
               paste(sep="", "FDR <= ", alpha, ": ", sum(multitest.table[,"fdr"] <= alpha)),
               paste(sep="", "E-value <=", alpha, ": ", sum(multitest.table[,"e.value"] <= alpha)),
               paste(sep="", "E-value <= 1: ", sum(multitest.table[,"e.value"] <= 1))
             ),
             pch=legend.pch, 
             cex=legend.cex,
             lwd=as.vector(legend.table$lwd), 
             lty=c(as.vector(legend.table$lty), "blank", "blank", "blank", "blank"),
             col=c(as.vector(legend.table$col), "white","white","white","white"),
             bty="o", bg="white")
    } else {
      legend(legend.corner, 
             #              c(paste(sep="", control.type," = ", alpha),
             #                paste(sep="", sum(positive), " positives"),
             #                paste(sum(!positive), "negatives")
             #              ),
             legend=as.vector(legend.table$legend),
             pch=as.vector(legend.table$pch), 
             cex=legend.cex,
             lwd=as.vector(legend.table$lwd), 
             lty=as.vector(legend.table$lty),
             col=as.vector(legend.table$col),
             bty="o", bg="white")
    }
  }
  # return(multitest.table.sorted) ## For debugging
}
