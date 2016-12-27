#' @title Given a vector of feature statuses (and optional positions), return a vector 
#' of associated colors.
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Given a vector of feature statuses (and optional positions), return a vector of associated colors, reflecting their status, and optionally the local density (densCols).
#'
#' @details
#' First version: 2015-11 
#' Last modification: 2015-11 
#'
#' @param feature.status vector of strings, each value indicating a status (e.g. "positive", "negative", "TP", "FP", "TN", "TP").
#' @param positions=NULL data frame indicating the position of each feature, which is 
#' passed to densCols() in order to tune the color according to local density of features. 
#' This requires to have defined, for each status, two colors named status.light and 
#' status.dens (e.g. "positive.light", "positive.dens" if the status is positive).
#' @param status.colors=c("default"="grey","light.dens"="#8888FF","dark.dens"="#FF6666","positive"='#4455DD', "positive.light" = '#88CCFF', "positive.dense" = "#2244BB","negative"="#BBBB00", "negative.light" = "#DDDD00", "negative.dense" = "#666600","TP"="#00DD00", "TP.light" = "#00FF00", "TP.dense" = "#004400","TN"="#888888", "TN.light" = "#CCCCCC", "TN.dense" = "#444444","FP"="#FF0000", "FP.light" = "#FF8888", "FP.dense" = "#880000","FN"="#FF8822", "FN.light"="#FFCC44", "FN.dense" = "#884411","alpha"="#BB2222","points"="#BBBBBB") 
#' Vector of colors associated to typical statuses.
#'
#' @examples
#' ## Create features with status-specific distributions of positions
#' TP <- data.frame(x=rnorm(n=2000, mean=0, sd=1), y=rnorm(n=1000, mean=0, sd=1), status="TP")
#' FP <- data.frame(x=rnorm(n=2000, mean=3, sd=1), y=rnorm(n=1000, mean=0, sd=1), status="FP")
#' FN <- data.frame(x=rnorm(n=2000, mean=0, sd=1), y=rnorm(n=1000, mean=3, sd=1), status="FN")
#' TN <- data.frame(x=rnorm(n=2000, mean=3, sd=1), y=rnorm(n=1000, mean=3, sd=1), status="TN")
#' features <- rbind(TP, TN, FP, FN)
#' table(features$status)
#' 
#' ## Assign colors according to status, irrespective of position
#' features$color <- featureColors(features$status)
#' table(features$status, features$color)
#' plot(features$x, features$y, col=features$color)
#' 
#' features$color <- featureColors(features$status, positions=features[, c("x", "y")])
#' plot(features$x, features$y, col=features$color)
#' 
#' @export
featureColors <- function (
  feature.status,
  positions=NULL,
  status.colors=c("default"="grey", 
                  "light.dens"="#8888FF",
                  "dark.dens"="#FF6666",
                  "positive"='#4455DD', "positive.light" = '#4466FF', "positive.dense" = "#2244BB",
#                  "positive"='#0033FF', "positive.light" = '#0099FF', "positive.dense" = "#0033CC",
                  "negative"="#888888", "negative.light" = "#BBBBBB", "negative.dense" = "#222222",
                  "TP"="#00DD00", "TP.light" = "#00FF00", "TP.dense" = "#004400",
                  "TN"="#00CCFF", "TN.light" = "#00FFFF", "TN.dense" = "#004444",
                  "FP"="#FF0000", "FP.light" = "#FF0000", "FP.dense" = "#440000",
                  "FN"="#FF6600", "FN.light"="#FFCC22", "FN.dense" = "#BB4411",
                  "alpha"="#BB2222",
                  "points"="#BBBBBB")) {
  
  ## Initialisation: assign default color, to make sure that each feature has a color
  if (is.null(positions)) {
    if (is.null(status.colors["default"])) {
      status.colors["default"] <- "grey"
    } 
    feature.colors <- rep(x = status.colors["default"], times=length(feature.status))
  } else {
    feature.colors <- densCols(positions)  
  }
  
  ## Get the list of statuses associated to the input features
  statuses <- sort(unique(feature.status))
  
  ## Assign one color to each status
  s <- 0 ## Initialize status counter
  for (status in statuses) {
    s <- s + 1 ## Increment status counter
    selected.features <- (feature.status == status) ## Select the features for current status
    if ((!is.null(positions)) &&
        (!is.null(status.colors[paste(sep=".", status, "light")])) &&
        (!is.null(status.colors[paste(sep=".", status, "dense")]))) {
      dens.palette <-colorRampPalette(c(status.colors[paste(sep=".", status, "light")], 
                                      status.colors[paste(sep=".", status, "dense")])) 
      feature.densities <- densCols(positions, col=dens.palette)
      feature.colors[selected.features] <- feature.densities[selected.features] 
    } else if (!is.null(status.colors[status])) {
      feature.colors[selected.features] <- status.colors[status]
    } else {
      feature.colors[selected.features] <- s
    }
  }
#   dens.palette <-colorRampPalette(status.colors[c("light.dens","dark.dens")]) 
#   dens.palette.pos <-colorRampPalette(status.colors[c("positive.light","positive.dense")]) 
#   dens.palette.neg <-colorRampPalette(status.colors[c("negative.light","negative.dense")]) 
#   dens.palette.TP <-colorRampPalette(status.colors[c("TP.light","TP.dense")]) 
#   dens.palette.TN <-colorRampPalette(status.colors[c("TN.light","TN.dense")]) 
#   dens.palette.FP <-colorRampPalette(status.colors[c("FP.light","FP.dense")]) 
#   dens.palette.FN <-colorRampPalette(status.colors[c("FN.light","FN.dense")]) 
  return(feature.colors)
}

