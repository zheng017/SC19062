#' @title A random walk Metropolis sampler using R
#' @description A random walk Metropolis sampler using R
#' @param sigma the sd
#' @param x0 the initial value
#' @param N the number of samples
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' N <- 2000
#' sigma <- c(.05, .5, 2, 16)
#' set.seed(1234)
#' x0 <- 25
#' rw1 <- rwMetropolisR(sigma[1], x0, N)
#' rw2 <- rwMetropolisR(sigma[2], x0, N)
#' rw3 <- rwMetropolisR(sigma[3], x0, N)
#' rw4 <- rwMetropolisR(sigma[4], x0, N)
#' rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
#' print(c(rw1$k, rw2$k, rw3$k, rw4$k))
#' print(c(1-rw1$k/N,1-rw2$k/N,1-rw3$k/N,1-rw4$k/N))
#' }
#' @export
rwMetropolisR <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= exp(-((abs(y)) - (abs(x[i-1])))))
      x[i] <- y
    else {
      x[i] <- x[i-1]
      k <- k+1
    }
  }
  return(list(x = x, k = k))
}
