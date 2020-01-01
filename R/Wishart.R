#generate a random sample from Wishart distribution
x<-function(sigma,n){
  d<-dim(sigma)[1]  #compute dimension d
  L<-t(chol(sigma)) #compute Choleski factorization of sigma
  #generate matrix T based on Bartlett's decomposition
  z<-rnorm(d^2)
  T<-matrix(z,d,d)
  T[upper.tri(T)]=0
  diag<-numeric(d)
  for (i in 1:d) {
    diag[i]<-rchisq(1,n-i+1)
    T[i,i]=diag[i]
  }
  A=T%*%t(T)
  #generate a random sample
  sample<-L%*%A%*%t(L)
  sample
}