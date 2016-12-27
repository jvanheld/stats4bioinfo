#' @title Run Wilcoxon rank-sum test (Mann-Whitney) on each row of a data frame
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Apply Wilcoxon rank-sum test (also called Mann-Whitney test) to each row of a 
#' data frame containing multivariate data for two samples. 
#' 
#' @details
#' First version: 2003-09
#' Last modification: 2015-02
#'
#' @param x                 A matrix or data frame
#' @param cl                A vector describing class assignment (length should equal the number of columns of the data table)
#' @param P.threshold=NA       p-value threshold. If specified, the result table only contains rows passing this threshold.
#' @param E.threshold=NA       e-value threshold. If specified, the result table only contains rows passing this threshold.
#' @param FDR.threshold=NA     Threshold on the False Discovery Rate (FDR). If specified, the result table only contains rows passing this threshold.
#' @param robust.est=F      Use robust estimators for central tendency and dispersion
#' @param verbosity=1       Level of verbosity
#' @param volcanoPlot=T     Draw a volcano plot.
#' @param alternative="two-sided"    Alternative hypothesis for the wilcox.test. Supported: "two.sided" (default), "less", "greater".
#' @param test.group=cl[1] Specify which group should be considered as first term for the difference (\eqn{d=m_{test}-m_{others}}). 
#' By default the first label of the class vector (cl) is considered as test group.
#' @param group.col.name   Include group labels in the column name of the output table
#' @param alpha=0.05  threshold to declare a feature positive. The threshold can be applied on any of the following statistics: 
#' @param ...  Additional parameters are passed to the function tTestPerRow.plotVolcano()
#' 
#' @return
#' A data.frame with one row per test result, and one column per statistics.
#' 
#' @examples
#' 
#' ## Load example data set from Den Boer, 2009
#' library(denboer2009)
#' data(denboer2009.expr)     ## Load expression table
#' data(denboer2009.pheno)    ## Load phenotypic data
#' data(denboer2009.group.labels)    ## Load phenotypic data
#'
#' ## Print cancer types and associated group labels
#' print(data.frame(denboer2009.group.labels))
#'
#' ## Compute the number of samples per subtype of cancer (ALL)
#' sort(table(denboer2009.pheno$sample.labels), decreasing=TRUE)
#'
#' ## Create a vector with group labels per sample,
#' ## For the Welch test we compare one group of interest (e.g. Bh)
#' ## to all the other ones (labeled as "other").
#' goi <- "Bh" ## Group of interest
#' sample.groups <- denboer2009.pheno$sample.labels
#' sample.groups[sample.groups != goi] <- "other"
#'
#' ## Check number of samples per group
#' sort(table(sample.groups))
#'
#' ## Run Welch test on each row of the DenBoer dataset
#' system.time(wilcox.result <- wilcoxTestPerRow(x=denboer2009.expr, cl=sample.groups, test.group="Bh"))
#' 
#' ## Draw a volcano plot with Welch result table
#' VolcanoPlot(multitest.table=wilcox.result$table, control.type="e.value", alpha=0.05, 
#'      effect.size.col="U.diff", xlab="U1 - U2",
#'      main=paste(sep="", "Wilcoxon rank-sum test: Den Boer (2009), ", goi, " vs others"),
#'      legend.corner = "topleft")
#' 
#' ## Run Welch test on each row of the DenBoer dataset
#' welch.result <- tTestPerRow(x=denboer2009.expr, cl=sample.groups, test.group="Bh", var.equal=FALSE)
#' 
#' ## Compare e-values from Wilcoxon and Welch tests
#' plotPvalCompa(data.frame(
#'    "Wilcoxon"=wilcox.result$table$e.value,
#'    "Welch"=welch.result$table$e.value), score="e-value", alpha=0.05)
#'    
#' ## Compare FDR from Wilcoxon and Welch tests
#' plotPvalCompa(data.frame(
#'    "Wilcoxon"=wilcox.result$table$fdr,
#'    "Welch"=welch.result$table$fdr), score="FDR", alpha=0.05,
#'    main="Wilcoxon versus Welch (FDR)")
#'    
#' ## Confusion table between Wilcoxon and Welch tests
#' table(wilcox.result$table$fdr < 0.05, welch.result$table$fdr < 0.05) ## Lenient threshold on FDR
#' table(wilcox.result$table$e.value < 0.05, welch.result$table$e.value < 0.05) ## Lenient threshold on E-value
#' table(wilcox.result$table$e.value < 1, welch.result$table$e.value < 1) ## Intermediate threshold on E-value
#' 
#' @export
wilcoxTestPerRow <- function (x,
                         cl,
                         P.threshold=NA,
                         E.threshold=NA,
                         FDR.threshold=NA,
                         robust.est=F, ## use robust estimators for central tendency and dispersion
                         verbosity=1, ## Level of verbosity
                         volcanoPlot = FALSE, ## Draw a volcano plot
                         alternative = "two.sided", ## Supported: "two.sided", "less", "greater"
                         group.col.names=FALSE, ## Include the group label in the column names
                         test.group=cl[1], ## Test group to compute the mean difference $d = m_{test} - m_{others}$
                         alpha=0.05, ## Alpha value, used to define the confidence interval width
                         native.test = TRUE, ## Run the native function stats::wilcox.test()
                         ... ## additional parameters are passed to the function tTestPerRow.plotVolcano()
) {
  
  
  ## Dimensions of the data set
  n <- nrow(x)
  p <- ncol(x)
  
  ## check the dimensions of the class vector
  if (length(cl) != p) {
    stop (paste(sep='', 
                'Invalid column number (',p,')',
                'should equal the length of the class vector (',length(cl),')'))
  }
  
  ## Compute classical or robust estimators of the central tendency (mean) 
  ## and dispersion (variance sd).
  group.sizes <- table(cl)
  groups <- unique(cl)
  if (length(groups) != 2) {
    stop("wilcoxTestPerRow: invalid vector of groups. Should contain exactly two different groups.")
  }
  control.group <- unique(groups[groups != test.group])
  
  ## Report starting time and parameters
  test.name <- "Wilcoxon"
  verbose(paste(sep="", "wilcoxTestPerRow: ",test.name," t-test on ", nrow(x), " rows. ", 
                "Group sizes: ", group.sizes[1], ",", group.sizes[2], ". ",
                test.group, " vs ", control.group), 1)

  ## Compute Wilcoxon statistics  
  n1 <- sum(cl == test.group)
  n2 <- sum(cl == control.group)
  N <- n1 + n2
  ranks <- data.frame(t(apply(x, 1, rank)))
  head(ranks)
  R1 <- apply(ranks[,cl==test.group], 1, sum) ## Sum of ranks for first group
  R2 <- apply(ranks[,cl==control.group], 1, sum) ## Sum of ranks for first group
  U1 <- R1 - n1*(n1+1)/2
  U2 <- R2 - n2*(n2+1)/2
  
  ## Run R wilcox.test() function on each row to compute the p-value
  # if (run.wilcox.t) {
  verbose("Computing wilcoxon p-values with stats::wilcox.test()", 1)
  return.wilcox <- function (x, cl, test.group, control.group, alternative="two.sided") {
    w <- wilcox.test(
      x[cl==test.group], x[cl==control.group], 
      alternative=alternative)
    return(c(pval=w$p.value, w$statistic))
  }
  system.time(wilcox.stats <- data.frame(t(apply(x, 1, return.wilcox, cl,
                                                 test.group=test.group, 
                                                 control.group=control.group, 
                                                 alternative=alternative))))
  W <- wilcox.stats$W
  p.value <- wilcox.stats$pval
  # } else {
  #   ## NOTE: the raw computation of wilcoxon statistics is fine but the p-values obtained with pwilcox 
  #   ## may be incorrect if there are ties in the dataset, which is frequently the case with RNA-seq data for example. 
  #   ## Therefore I prefer to use the native R function stats::wilcox.test(), which handles ties by using a normal approximation. 
  #   ## 
  #   ## NOTE: the function stats::wilcox.test() gives slightly different results from the 
  #   ## apply() computation applied below. This might be due to the fact that the default 
  #   ## wilcox test uses a normal approximation. The option exact=TRUE cannot be used when 
  #   ## there are ties, which is frequently the case with RNA-seq counts for example.
  #   ##
  #   ## I guess it is safer to run the native wilcox.test(), which will treat the ties. On the other hand, 
  #   ## its p-value is restricted to 1e-32, whereas the log.p computation below sets the limit 
  #   ## much lower. Of course the esimtates are too imprecise to trust such precisions, but the 
  #   ## p-value order of magniture is useful to rank genes (even knowing that there might be huge 
  #   ## imprecision). I thus leave it as an option. 
  #   ## 
  #   verbose("Computing wilcoxon p-values with pwilcox(log.p=TRUE)", 1)
  #   
  #   ## Compute Wilcoxon statistics
  #   W <- U1
  #   #W <- pmin(U1,U2)
  #   # hist(W, breaks=50)
  #   
  #   ## Calculate p-value and e-value
  #   #  p.value.normal.approx <- 2*pnorm(abs(W),lower.tail=F)
  #   if (alternative == "greater") {
  #     #    p.value <- pwilcox(W-0.1, m=n1, n=n2, lower.tail=FALSE)
  #     U <- U1
  #     wm <- rep(n1, length.out = n)
  #     wn <- rep(n2, length.out = n)
  #     system.time(p.value <- exp(pwilcox(U, m=wm, n=wm,lower.tail=TRUE, log.p=TRUE)))
  #   } else if (alternative == "less") {
  #     #    p.value <- pwilcox(W, m=n1, n=n2, lower.tail=TRUE)
  #     U <- U2
  #     wm <- rep(n2, length.out = n)
  #     wn <- rep(n1, length.out = n)
  #     system.time(p.value <- exp(pwilcox(U, m=wm, n=wm,lower.tail=TRUE, log.p=TRUE)))
  #   } else if (alternative == "two.sided") {
  #     p.value.lower <- pwilcox(W, m=n1, n=n2,lower.tail=TRUE)
  #     p.value.higher <- pwilcox(W-0.5, m=n1, n=n2,lower.tail=FALSE)
  #     p.value <- 2*pmin(p.value.lower, p.value.higher)
  #     U <- pmin(U1,U2)
  #     wm <- rep(n1, length.out = n); wm[U2<U1] <- n2
  #     wn <- rep(n2, length.out = n); wn[U2<U1] <- n1
  #     system.time(p.value <- 2*exp(pwilcox(U, m=wm, n=wm,lower.tail=TRUE, log.p=TRUE)))
  #   } else {
  #     stop('Invalid alternative option for wilcoxTestPerRow(). Supported: "two.sided", "less", or "greater"')
  #   }
  # } 
  # range(p.value)
  # hist(p.value, breaks=seq(0,1,0.05)) ## For debug
  # compa <- data.frame(cbind(wilcox.stats, U1, U2, W2=W, U, wm, wn, U1.ge.U2 = U1>U2, p.value))
  # head(compa)
  # plot(compa$p.value2, compa$p.value, log="xy")
  # table(compa$pval==compa$p.value)
  # table(compa$p.value2==compa$p.value)
  
  # wt <- wilcox.test(unlist(x[2,cl==test.group]), unlist(x[2,cl==control.group]))
  # wt$p.value
  
  ## Multiple testing corrections
  e.value <- p.value*nrow(x)
  sig <- -log(e.value)/log(10)
  
  multi.corr <- multipleTestingCorrections(p.value, plots=FALSE)
  
  
  ## Build the result
  result.table <- data.frame(
    R1,
    R2,
    R.diff = R1 - R2,
    U1,    
    U2,
    U.diff = U1 - U2,
    W,
    n1,
    n2,
    p.value,
    e.value,
    sig)
  result.table$fwer <- multi.corr$multitest.table$fwer
  result.table$q.value <- multi.corr$multitest.table$qval.Storey
  result.table$fdr <- multi.corr$multitest.table$fdr
  result.table$rank <- multi.corr$multitest.table$rank
  
  ## Adapt column headers to the type of estimators, and to the group names
  if (group.col.names) {
    names(result.table) <- sub(
      x=names(result.table), perl = TRUE,
      pattern="1$", 
      replacement = paste(sep="", ".", test.group)
    )
    names(result.table) <- sub(
      x=names(result.table), perl = TRUE,
      pattern="2$", 
      replacement = paste(sep="", ".", control.group)
    )
  }
  head(result.table)
  
  ## Filtering on p-value and e-value thresholds
  if (!is.na(P.threshold)) {
    result.table <- result.table[result.table$p.value < P.threshold,]
  }
  if (!is.na(E.threshold)) {
    result.table <- result.table[result.table$e.value < E.threshold,]
  }
  if (!is.na(FDR.threshold)) {
    result.table <- result.table[result.table$fdr < FDR.threshold,]
  }
  
  ## Report done time
  verbose("Multiple Wilcoxon test done", 2)
  
  ## Draw a volcano plot
  if (volcanoPlot) {
    tTestPerRow.plotVolcano(result.table, ...)
  }
  
  ## Build the result object with result table + parameters of the analysis
  result <- list()
  result$table <- result.table
  result$alpha <- alpha
  result$groups <- groups
  result$test.group <- test.group
  result$control.group <- control.group
  result$p.threshold <- P.threshold
  result$e.threshold <- E.threshold
  result$fdr.threshold <- FDR.threshold
  result$alternative <- alternative
  
  ##  plot(p.value.normal.approx,p.value,panel.first=grid(col='#0000ff'),log="xy")
  return (result)
}

