dCov <- function(x, y) {
  x <- as.matrix(x); y<- as.matrix(y)
  n <- nrow(x); m<- nrow(y)
  if (n != m || n < 2)  stop('Sample sizes must agree')
  if (! (all(is.finite(c(x,y)))))
    stop('Data contains missing or infinite values')
  Akl <- function(x) {
    d <- as.matrix(dist(x))
    m <- rowMeans(d); M <- mean(d)
    a <- sweep(d, 1, m); b <- sweep(a, 2, m)
    b + M
  }
  A <- Akl(x); B <- Akl(y)
  sqrt(mean(A * B))
}