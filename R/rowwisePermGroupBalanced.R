#' @title Permute rows in a group-balanced way
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' 
#' @description 
#' #' Random permutation of the values on each row of the input data frame, 
#' with a preservation of the proportions between 2-group labeled columns:
#' group of interest (GOI), and  "others", respecively.
#'
#' Each row thus contains exactly the same values as in the original expression
#' matrix, but there should be no specific distinction between groups.
#'
#' ATTENTION ! This procedure may give rise to surprizing bias. 
#' A priori it seemed to me that a balanced representation of the
#' original groups between the permuted samples would be a good 
#' idea, because I sometimes observed that permutation tests with small 
#' sample sizes would return too many significant results since sometimes 
#' the resampled groups contain different proportions of the original
#' groups. I thus implemented balanced permutation to suppress this 
#' effect. 
#' 
#' However, I noticed the opposite effect: when the effect size
#' is very strong in a given dataset, the balanced permuted set 
#' has an *under-representation* lof low p-values 
#' (e.g. 0 <= pval <= 30%), see example below. The cause of this 
#' surprizing behaviour is that the balanced permutations ensure
#' equality of the resampled group means, but if the original groups
#' have very different means, the resampled distributions are bimodal,
#' and have thus a high variance. The consequence is to artificially 
#' reduce the denominator of the t statistics ($t_{obs}$).
#'
#' 
#' @details
#' First version: 2015-03
#' Last modification: 2015-03
#'
#' @param x       A matrix or data frame
#' @param g       Group labels
#' @param goi     Group of interest. If not specified, the first label is take as group of interest.
#' In case g contains more than two distinct labels, these group labels are converted
#' to "GOI", and "other", respectively.
#'
#' @return
#' A data frame of same dimensions as the input matrix/data frame, with row-wise
#' group-balanced permuted values.
#'
#' @examples
#' ## Run example for rowwiseSample, in order to load 
#' ## the data and parameters
#' example(rowwiseSample)
#'
#' ## Permute the values of denboer2009
#' balanced.perm.profiles <- rowwisePermGroupBalanced(
#'     x=denboer2009.expr[selected.samples],
#'     g=selected.labels,
#'     goi=group1)
#'
#' ## Run Welch test on the row-wise permuted values
#' balanced.perm.welch <- meanEqualityTests(
#'     balanced.perm.profiles, 
#'     g=selected.labels, goi=group1,
#'     selected.tests="welch")
#'                  
#' ## Draw volcano plot of Welch test result with the permuted values, resp.
#' ## NOTE: we already see that the negative control is "too good": 
#' ## the highest significances are at -2 instead of 0.
#' meanEqualityTests.plotVolcano(balanced.perm.welch, 
#'    test="welch", 
#'    legend.corner='topright',
#'    main="Permuted Den Boer 2009, Welch volcano")
#' 
#' ## Plot p-value distribution for the balanced row-wise permuted dataset
#' ## NOTE: this plot clearly shows the bias of balanced permutation:
#' ## the low p-values (<= 30%) are under-represented because when the 
#' ## population means differe, the balanced resampling creates groups
#' ## with same expected mean, but increaseed variance.
#' mulitpleTestingCorrections.plotPvalDistrib(
#'    balanced.perm.welch$welch.multicor, legend.corner="bottomright",
#'    col='#FFBBBB')
#' 
#' @export
rowwisePermGroupBalanced <- function (x,
                                      g,
                                      goi=g[1]) {

  verbose(paste(sep="", "Generating group-balanced row-wise random normal matrix (", 
                nrow(x), " rows x ", ncol(x), " columns)"), 1)
  goi.vs.others.labels <- g
  goi.vs.others.labels[goi.vs.others.labels !=goi] = "other"

  ## Define the number of elements to sample in each group
  N <- ncol(x) ## Number of columns of input table
  n <- sum(g==goi)
  m <- N - n     ## Number of entries from other groups

  verbose(paste(sep="", "Group sizes: ", n, " versus ", m, " (N=",N, ")"), 1)
  
  ## Use a trick to randomly choose between ceiling and floor for the half of odd numbers
  n1.1 <- round((n + runif(1,-0.1,+0.1))^2/N) ## Number to samples from group 1 that remain in group 1
  n1.2 <- n - n1.1 ## Number of samples from group 2 assigned to group 1
  n2.1 <- n - n1.1 ## Number of samples from group 1 assigned to group 2
  n2.2 <- m - n1.2 ## Number of samples from group 2 remaining in group 2

  ## Define the function that will be applied to each row
  permuteOneRowGroupBalanced <- function(x, n) {
    ## We resample the sampled subgroups to break their order
    g1 <- sample(c(sample(1:n, size=n1.1), sample((n+1):N, size=n1.2)))
    g2 <- sample(c(sample(1:n, size=n2.1), sample((n+1):N, size=n2.2)))
    permuted <- x[c(g1, g2)]
    return(permuted)
  }

  return(t(apply(x, 1, permuteOneRowGroupBalanced, n=n)))
}


