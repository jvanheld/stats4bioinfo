% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meanEqualityTests.R
\name{meanEqualityTests}
\alias{meanEqualityTests}
\title{Multiple test for mean equality}
\usage{
meanEqualityTests(x, g, goi = NULL, selected.tests = NULL, alpha = 0.05,
  verbosity = 0)
}
\arguments{
\item{x}{A matrix or data frame}

\item{g}{A vector describing group assigment
(length should equal the number of columns of the input matrix)}

\item{goi}{Group of interest. Required for the two-groups tests (Student, Welch, Wilcoxon)}

\item{alpha}{Significance threshold (will be applied on corrected p-values)}

\item{verbosity}{level of verbosity (print messages during execition)}

\item{selected.tests=NA}{Selection of a subset of tests to run.}
}
\value{
\describe{
  \item{param}{Parameters of the analysis}
  \item{selected.tests}{Tests ran. If no test was specified, all possible tests are selected.}
  \item{nb.per.group}{A sorted vector indicating the number of individuals (columns) per group.}
  \item{stats.per.row}{data.frame summarizing some descriptive statistics + the significance scores (p-value, e-value, fdr) returned by the
  different tests for each row of the input matrix.}
}
}
\description{
Run different tests of mean equality (parametric or not)
on each row of a data table, and return a table summarizing the results
returned by each test.

The function also runs test to check the conditions of applicability 
for the mean equality tests: 
(1) Test of normality: Shapiro; 
(2) Tests of variance equality: F-test  (parametric), 
    Levene (non-parametric), 
    Brown-Forsythe (non-parametric, robust estimators).
}
\details{
First version: 2015-03
Last modification: 2015-03
}
\examples{
## Load example data set from Den Boer, 2009
library(denboer2009)
data(denboer2009.expr)     ## Load expression table
data(denboer2009.pheno)    ## Load phenotypic data
data(denboer2009.amp)    ## Load absent/marginal/present calls to filter input set

## Select subtypes represented by at least 30 samples
verbose("Selecting subtypes with at least 30 samples")
samples.per.subtype <- table(denboer2009.pheno$sample.labels)
selected.subtypes <- names(samples.per.subtype)[samples.per.subtype >= 30]
selected.samples <- denboer2009.pheno$sample.labels \%in\% selected.subtypes

## Define group labels and group of interest
g <- denboer2009.pheno[selected.samples, "sample.labels"]
goi <- "Bh" ## Select one cancer subtype as group of interest

## Filter out genes called absent in most samples
verbose("Selecting probeset present in at least 30 samples")
selected.probesets <- apply(denboer2009.amp == "P", 1, sum) >= 30
x <- denboer2009.expr[selected.probesets, selected.samples]
verbose(paste("Matrix size:", nrow(x), "x", ncol(x)))

## Run several mean equality tests on each row of the data set.
## We restrict it to the first probesets for the demo.
diff.results <- meanEqualityTests(x=x, g=g, goi=goi, 
    selected.tests=c("welch", "wilcoxon"), verbosity=1)

## Return the parameters of the analysis
print(diff.results$param)
print(diff.results$selected.tests)
print(diff.results$nb.per.group)

## Compare p-values returned by 2 tests
welch.vs.wilcoxon <- meanEqualityTests.compareTwoTests(
   diff.results, test1="welch", test2="wilcoxon")

## Display the names of the result field stat.per.row
names(diff.results$stats.per.row)

## Draw a color Volcano plot for Welch test results
meanEqualityTests.plotVolcano(diff.results, test="welch", legend.cex=0.7, plot.cex=0.5)

## Draw a grayscale Volcano plot for Welch test results
meanEqualityTests.plotVolcano(diff.results, test="welch", legend.cex=0.7, plot.cex=0.5, plot.colors=FALSE)
}
\author{
Jacques van Helden (\email{Jacques.van-Helden@univ-amu.fr})
}

