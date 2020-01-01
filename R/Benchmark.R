#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of C functions (\code{rwMetropolisR}) and Cpp functions (\code{rwMetropolisC}).
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm rgamma
#' @useDynLib SC19062
#' @examples
#' \dontrun{
#' ts <- microbenchmark(rwC1 = rwMetropolisC(sigma[1], x0, N), rwC2 = rwMetropolisC(sigma[2], x0, N), rwC3 = rwMetropolisC(sigma[3], x0, N), rwC4 = rwMetropolisC(sigma[4], x0, N), rw1 = rwMetropolisR(sigma[1], x0, N), rw2 = rwMetropolisR(sigma[2], x0, N), rw3 = rwMetropolisR(sigma[3], x0, N), rw4 = rwMetropolisR(sigma[4], x0, N))
#' knitr::kable(summary(ts)[, c(1,3,5,6)])
#' }
NULL
