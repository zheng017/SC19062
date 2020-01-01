#include <Rcpp.h>
using namespace Rcpp;

//' @title A random walk Metropolis sampler using R
//' @description A random walk Metropolis sampler using R
//' @param sigma the sd
//' @param x0 the initial value
//' @param N the number of samples
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' N <- 2000
//' sigma <- c(.05, .5, 2, 16)
//' x0 <- 25
//' rw1 <- rwMetropolisC(sigma[1], x0, N)
//' rw2 <- rwMetropolisC(sigma[2], x0, N)
//' rw3 <- rwMetropolisC(sigma[3], x0, N)
//' rw4 <- rwMetropolisC(sigma[4], x0, N)
//' rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
//' print(c(rw1$k, rw2$k, rw3$k, rw4$k))
//' print(c(1-rw1$k/N,1-rw2$k/N,1-rw3$k/N,1-rw4$k/N))
//' }
//' @export
// [[Rcpp::export]]
NumericVector rwMetropolisC(double sigma, double x0, double N) {
	NumericVector x(N+1);
	x[0] = x0;
	NumericVector u = runif(N);
	double k = 0;
	double y = 0;
	for (int i = 2; i < N+1; i++) {
		y = rnorm(1, x[i-2], sigma)[0];
		if (u[i-2] <= exp(-((abs(y))-(abs(x[i-2])))) ) {
			x[i-1] = y;
		} 
		else {
			x[i-1] = x[i-2];
			k++;
		}
	}
	x[N] = k;
	return(x);
}
















