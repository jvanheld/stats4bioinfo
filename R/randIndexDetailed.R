#' @title Compute Rand Contingency Index from a contingency table.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#' @description Compute Rand Contingency Index from a contingency table.. 
#'
#' @details
#' First version: 2016-12-26
#' Last modification: 2016-12-26
#' 
#' The Rand Consistency Index is a classical measure of the similarity between 
#' two classifications (e.g. clustering results). It relies on the number of 
#' consistent pairs, where "consistent" means either co-clusttered in both 
#' classifications or separated in both.
#'
#' @param cont.table a contingency table (supposed to contain counts). Can be produced by table() or other means. Can belong to any of the following classes: table, matrix or data.frame.
#' @param permute.cols=FALSE permute columns of the contingency table by decreasing sum.
#' @param permute.rows=FALSE permute rows of the contingency table by decreasing sum.
#' @param ...  Additional parameters are passed to the function image()
#' 
#' @return a list with the values of the computed parameters. 
#' @examples
#' clustering1 <- c(rep(c("a", "b","c"), times=c(4,2,1)))
#' clustering2 <- c(rep(c("X", "Y"), times=c(3,4)))
#' cont.table <- table(clustering1, clustering2)
#' 
#' ## Example from Yeung and Ruzzo (2001) supplementary file
#' cont.table <- data.frame(v1 = c(1, 1, 0), v2 = c(1, 2, 0), v3=c(0, 1, 4), row.names=c("u1", "u2", "u3"))
#' result <- randIndexDetailed(cont.table)
#' print(result$consistency.table)
#' print(signif(unlist(result[c("n", "Npairs", "RI", "ARI")]), digits=3))
#' 
#' ## Check the result with another implementation
#' flexclust::randIndex(as.table(as.matrix(cont.table)), correct=TRUE)
#' 
#' ## Draw heatmaps of the contingency tables
#' heatmap.simple(cont.table)
#' heatmap.simple(result$cont.table)
#' @export
randIndexDetailed <- function(cont.table,
                      permute.cols = FALSE,
                      permute.rows = FALSE) {
  
  # Cast input contingency table if required
  if (class(cont.table) %in% c("table", "matrix'")) {
    cont.table <- as.data.frame.matrix(cont.table)
  }
  if (class(cont.table) != "data.frame") {
    stop("Invalid class for contingency table: ", class(cont.table))
  }
  
  # Permute columns and rows of the contingency table if requested
  if (permute.rows) {
    cont.table <- cont.table[order(rowSums(cont.table), decreasing = TRUE),]
  }
  if (permute.cols) {
    cont.table <- cont.table[,order(colSums(cont.table), decreasing = TRUE)]
  }
  
  # Result object
  result <- list()
  result$cont.table <- cont.table
  result$cont.table$Sums <- rowSums(result$cont.table)
  result$cont.table["Sums",] <- colSums(result$cont.table)
  
  # Total number of elements in the classification
  result$n <- n <- sum(cont.table)
  
  
  # Initialize the consistency table
  result$consistency.table <- data.frame(matrix(nrow=3, ncol=3))
  names(result$consistency.table) <- c("co.clust", "separated", "total")
  row.names(result$consistency.table) <- c("co.clust", "separated", "total")
  
  # Total number of pairs between members of any cluster
  result$consistency.table[3,3] <- result$Npairs <- n * (n-1) /2
  
  # Number of co-clustered pairs per cell of the contingency table, 
  # which correspond to pairs co-clustered in both row cluster and 
  # column cluster.
  result$pairs.per.cell <- cont.table * (cont.table-1) /2
  result$n.per.row <- apply(cont.table, 1, sum)
  result$n.per.col <- apply(cont.table, 2, sum)
  
  # Number of co-clustered pairs for columns and rows, resp.
  result$pairs.per.row <- result$n.per.row * (result$n.per.row-1) / 2
  result$pairs.per.col <- result$n.per.col * (result$n.per.col-1) / 2
  
  # Random expectation for the contingency table
  result$exp.cont.table <- as.vector(result$n.per.row) %*% t(as.vector(result$n.per.col))/result$n
  
  # Consistent co-clustering
  result$consistency.table[1,1] <- sum(result$pairs.per.cell)
  
  # Total number of co-clustered pairs 
  #result$consistency.table[1,3] <- sum(result$n.per.row * (result$n.per.row -1) /2)
  result$consistency.table[1,3] <- sum(result$pairs.per.row)
  #result$consistency.table[3,1] <- sum(result$n.per.col * (result$n.per.col -1) /2)
  result$consistency.table[3,1] <- sum(result$pairs.per.col)
  result$consistency.table[1,2] <- result$consistency.table[1,3] - result$consistency.table[1,1]
  result$consistency.table[2,1] <- result$consistency.table[3,1] - result$consistency.table[1,1]
  result$consistency.table[3,2] <- result$consistency.table[3,3] - result$consistency.table[3,1]
  result$consistency.table[2,3] <- result$consistency.table[3,3] - result$consistency.table[1,3]
  result$consistency.table[2,2] <- result$consistency.table[3,2] - result$consistency.table[1,2]

  result$RI <- (result$consistency.table[1,1] + result$consistency.table[2,2]) / result$consistency.table[3,3]
  
  ## Compute adjusted rand index. 
  ## Details can be found here: http://faculty.washington.edu/kayee/pca/supp.pdf
  expectedIndex <- sum(result$pairs.per.col) * sum(result$pairs.per.row)/result$Npairs
  result$ARI <- 
    (sum(result$pairs.per.cell) - expectedIndex) /
    ((sum(result$pairs.per.col) + sum(result$pairs.per.row))/2 - expectedIndex)
  
  return(result)
}  
