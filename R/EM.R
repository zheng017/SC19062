#' @title A EM algorithm for Gaussian mixture models
#' @description A EM algorithm for Gaussian mixture models
#' @param x e
#' @param mu initial mean
#' @param sigma1 initial sigma1
#' @param sigma2 initial sigma2
#' @param prop initial proportion
#' @return estimates of 
#' @import mixtools
#' @examples
#' \dontrun{
#' mu<-matrix(c(2,4,50,80),2,2)
#' sigma1<-matrix(c(1,0,0,1),2,2)
#' sigma2<-matrix(c(1,0,0,1),2,2)
#' prop<-c(0.5,2)
#' EM(data,mu,sigma1,sigma2,prop)
#' }
#' @export
EM<-function(x,mu,sigma1,sigma2,prop)
{
  library(mixtools)
  data("faithful")
  data<-as.matrix(faithful)
  n<-nrow(data)
  fi<-function(x,mu,sigma) (2*pi)^(-n/2)/sqrt(det(sigma))*exp(-1/2*t(x-mu)%*%solve(sigma)%*%(x-mu))
  Z<-matrix(0,2,n)
  for (i in 1:20) {
    #compute Z
    for (k in 1:n) {
      Z[1,k]<-prop[1]*fi(data[k,],mu[1,],sigma1)
      Z[2,k]<-prop[2]*fi(data[k,],mu[2,],sigma2)
      Zsum<-Z[1,k]+Z[2,k]
      Z[1,k]<-Z[1,k]/Zsum
      Z[2,k]<-Z[2,k]/Zsum
    }
    #updata prop, mu and sigma
    prop[1]<-mean(Z[1,])
    prop[2]<-mean(Z[2,])
    mu[1,]<-Z[1,]%*%data/sum(Z[1,])
    mu[2,]<-Z[2,]%*%data/sum(Z[2,])
    sigma1<-matrix(0,2,2)
    for (k in 1:n) {
      sigma1star<-(Z[1,k]*(as.matrix(data[k,])-as.matrix(mu[1,]))%*%t(as.matrix(data[k,])-as.matrix(mu[1,])))/sum(Z[1,])
      sigma1<-sigma1+sigma1star
    }
    sigma2<-matrix(0,2,2)
    for (k in 1:n) {
      sigma2star<-(Z[2,k]*(as.matrix(data[k,])-as.matrix(mu[2,]))%*%t(as.matrix(data[k,])-as.matrix(mu[2,])))/sum(Z[2,])
      sigma2<-sigma2+sigma2star
    }
  }
  return (list(mu=mu,sigma1=sigma1,sigma2=sigma2,prop=prop))
}