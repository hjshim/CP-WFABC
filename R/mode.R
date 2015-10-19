mode <- function(x) {
  d <- density(x,adjust = 1)
  d$x[which.max(d$y)]
}