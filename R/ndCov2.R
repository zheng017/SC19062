ndCov2 <- function(z, ix, dims) {
  p <- dims[1]
  q <- dims[2]
  d <- p + q
  x <- z[ , 1:p]
  y <- z[ix, -(1:p)]
  return(nrow(z) * dCov(x, y)^2)
}