#' @title ROC statistics
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Compute the statistics required to draw a ROC curve, 
#' + some additional validation statistics
#' @details
#' First version: 2015-04
#' Last modification: 2015-04
#' @param neg.scores vector of scores obtained with the truly negative cases
#' @param pos.scores vector of scores obtained with the truly positive cases
#' @param neg.total=length(neg.scores) total number of negative cases
#' @param pos.total=length(pos.scores) total number of positive cases
#' @param sort.decreasing=FALSE  Indicate whether the score should be sorted by
#' decreasing order (this option should be used when lower scores are 
#' considered to be more significant, e.g. P-value, E-value, FDR).
#' @export
ROCstatistics <- function(neg.scores,
                          pos.scores,
                          neg.total=length(neg.scores),
                          pos.total=length(pos.scores),
                          sort.decreasing=FALSE) {
  all.scores <- sort(unique(c(neg.scores, pos.scores)),decreasing=sort.decreasing)
  n.scores <- length(all.scores)
  
  ROC.table <- data.frame(row.names=all.scores,
                          score=all.scores,
                          pos.counts=rep(0,n.scores)
  )
  
  ## Distribution of positive scores
  pos.counts <- table(pos.scores)
  ROC.table$pos.counts <- 0
  ROC.table[names(pos.counts),"pos.counts"] <- pos.counts
  ROC.table$pos.counts.icum <-  rev(cumsum(rev(ROC.table$pos.counts)))
  ROC.table$pos.iCDF <-  ROC.table$pos.counts.icum/pos.total
  
  ## Distribution of negative scores
  neg.counts <- table(neg.scores)
  ROC.table$neg.counts <- 0
  ROC.table[names(neg.counts),"neg.counts"] <- neg.counts
  ROC.table$neg.counts.icum <-  rev(cumsum(rev(ROC.table$neg.counts)))
  ROC.table$neg.iCDF <-  ROC.table$neg.counts.icum/neg.total
  
  return(ROC.table)
}


#' @title Compute the area under the ROC
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Compute the area under the curve from a table returned by ROCstatistics
#' @details
#' First version: 2015-04
#' Last modification: 2015-04
#' @param ROC.table A result from ROC.table
#' @param ... Additional parameters are passed to XXXX().
#' @export
AUCfromROCtable <- function(ROC.table) {
  x <- c(0, sort(ROC.table$neg.iCDF), 1)
  y <- c(0, sort(ROC.table$pos.iCDF), 1)
  n <- length(x)
  
  auc <- t(as.matrix(abs(x[2:n] - x[1:n-1])))%*%as.matrix(abs((y[2:n] + y[1:n-1])/2))
  auc <- as.vector(auc)
  return(auc)
}


################################################################
## Add labels on a ROC curve, for a given set of FPR values
labels.for.fpr <- function(ROC.table,
                           neg.values=c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4,0.5)) {  
  label.frame <- data.frame()
  for (neg in neg.values) {
    i <- which(ROC.table$neg.iCDF < neg)
    i <- i[1]
    pos <- ROC.table[i,"pos.iCDF"]
    score <- ROC.table[i,"score"]
    label.frame <- rbind(label.frame,
                         data.frame(i,pos,neg,score))
  }
  segments(label.frame$neg, 0,label.frame$neg,label.frame$pos,col="darkred",lty="dotted")
  text(label.frame$neg,0,labels=round(label.frame$neg*100),pos=3,font=2,col="darkred",srt=90)
  segments(0,label.frame$pos,label.frame$neg,label.frame$pos,col="darkred",lty="dotted")
  text(0,label.frame$pos,labels=round(label.frame$pos*100),adj=c(1,0.5),font=2,col="darkred",bg='white',bty="n")
  text(label.frame$neg,label.frame$pos,labels=label.frame$score,adj=c(0,1),font=2,col="darkred")
  return(label.frame)
}

