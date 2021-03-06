% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multipleTestingCorrections.plotSignifFeatures.R
\name{mulitpleTestingCorrections.plotSignifFeatures}
\alias{mulitpleTestingCorrections.plotSignifFeatures}
\title{Plot number of significant features as a function of the significance threshold.}
\usage{
mulitpleTestingCorrections.plotSignifFeatures(multitest.result,
  main = "Significant features", xlab = "Significance threshold",
  ylab = "Significant features", alpha = 0.05, plot.legend = TRUE,
  legend.corner = "topleft", legend.cex = 1, xlim = NULL,
  plot.pch = c(p.value = 2, fdr = 4, qval.0 = 3, e.value = 1, fwer = 20),
  plot.col = c(p.value = "#000000", fdr = "#888888", qval.0 = "#666666",
  e.value = "#BBBBBB", fwer = "#444444"), plot.elements = c("p.value", "fdr",
  "qval.0", "e.value", "fwer"), ...)
}
\arguments{
\item{multitest.result}{the list returned by the function multipleTestingCorrections().}

\item{...}{Additional parameters are passed to plot()}

\item{main="Significant}{features"  main title of the plot}

\item{xlab="P-value}{derived statistics"}

\item{ylab="Significant}{features"}

\item{alpha=0.05}{Threshold of significance (alpha).}

\item{plot.legend=TRUE}{Plot a legend indicating the number of features declared significant with the alpha threshold on the selected statistics.}

\item{legend.corner="topleft"}{corner wher the legend has to be placed.}

\item{legend.cex=1}{Font size for the legend.}

\item{plot.pch=c(p.value=2, e.value=1, fwer=20, fdr=4, qval.0=3)}{Specific characters to distinguish the plotted statistics.}

\item{plot.col=c(p.value='#000000', e.value='#BBBBBB', fwer='#444444', fdr='#888888', qval.0='#666666')}{Specific colors or gray levels to distinguish the plotted statistics.}

\item{plot.elements=c("p.value", "fdr", "qval.0", "e.value", "fwer")}{Selection of elements to display on the plot.}
}
\description{
Plot the number of significant features as a function of the control
criterion (nominal p-value, e-value, fdr, ...).
}
\examples{
## To obtain the input list (multitest.result), see the documentatio of
example(multipleTestingCorrections)

## Plot all the multiple testing corrections at once
mulitpleTestingCorrections.plotSignifFeatures(multitest.result)

## Compare e-value and FWER
mulitpleTestingCorrections.plotSignifFeatures(multitest.result, plot.elements=c("e.value","fwer"))

## Compare e-value and FDR
mulitpleTestingCorrections.plotSignifFeatures(multitest.result, plot.elements=c("fdr","e.value"))

## Compare Benjamini-Hochberg (qval.0) and Storey-Tibshirani (fdr) estimates of FDR
mulitpleTestingCorrections.plotSignifFeatures(multitest.result, plot.elements=c("fdr","qval.0"))
}
\author{
Jacques van Helden (\email{Jacques.van-Helden@univ-amu.fr})
}

