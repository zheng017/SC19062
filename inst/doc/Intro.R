## ----message=FALSE, warning=FALSE---------------------------------------------
library(SC19062)
#initial values of parameters
mu<-matrix(c(2,4,50,80),2,2)
sigma1<-matrix(c(1,0,0,1),2,2)
sigma2<-matrix(c(1,0,0,1),2,2)
prop<-c(0.5,2)
#outcomes
EM(data,mu,sigma1,sigma2,prop)

## ----message=FALSE, warning=FALSE---------------------------------------------
n=200
z1 <- rnorm(n)
z2 <- rnorm(n)
eps <- rnorm(n,0,0.2)
eps1 <- rnorm(n,0,0.2)
eps2 <- rnorm(n,0,0.2)
xx <- cbind(z1,z2)
xt <- runif(n,-6,6)
vv <- cbind(z1+eps1,z2+eps2)
y <- 2*xx[,1] + xx[,2] + cos(pi*xt/6) + eps
gg <- cos(pi*xt/6)
## ventana= smoothing parameter
## the kernel is scaled so that their quartiles
## (viewed as probability densities) are at +/- 0.25*ventana
## The smoothing bandwidth = 0.25*ventana/0.675
ventana=2.5
salida=PLMEV(y,vv,xt,palabra="mcd",cw=0.975,k=0.55,ventana)
salida$beta.rob
plot(xt,salida$ghat,xlab='T',ylab='g(t)')
points(xt,gg,col="red")
legend('topright' , legend = c('true' , 'robust orthogonal') , col = c('red', 'black') , lty = 3 , lwd = 5)

## -----------------------------------------------------------------------------
library(bootstrap)     
data(scor)
knitr::kable(scor)   #scan the data in table

## ----message=FALSE, warning=FALSE---------------------------------------------
library(TSA)
#win.graph(width=4.875,height=2.5,pointsize=8)
data(larain)
plot(larain,ylab='Inches',xlab='Year',type='o')

## -----------------------------------------------------------------------------
#win.graph(width=3,height=3,pointsize=8)
plot(y=larain,x=zlag(larain),ylab='Inches',xlab='Previous Year Inches')

## -----------------------------------------------------------------------------
library(SC19062)
data("data1")
cl<-kmeans(data,5)    #using k-means algorithm
plot(data,col=cl$cluster)

## -----------------------------------------------------------------------------
n<-1000
set.seed(1234)
u<-runif(n)   #generate random samples from U(0,1)
sigma<-c(1,2,3)
x<-sigma[1]*sqrt(-2*log(u))  #generate random samples when sigma=1
y<-sigma[2]*sqrt(-2*log(u))  #generate random samples when sigma=2
z<-sigma[3]*sqrt(-2*log(u))  #generate random samples when sigma=3
#plot histograms
hist(x,prob=TRUE,main=expression(sigma==1))
hist(y,prob=TRUE,main=expression(sigma==2))
hist(z,prob=TRUE,main=expression(sigma==3))

## -----------------------------------------------------------------------------
n<-1000
set.seed(1234)
x1<-rnorm(n,0,1)
x2<-rnorm(n,3,1)
p1<-c(0.75,0.8,0.7,0.65,0.6,0.55)
par(mfcol=c(2,3))  #prepare to plot histograms of different values of p1 in one figure
for (i in 1:6) {
  r<-sample(c(0,1),n,replace=TRUE,prob=c(1-p1[i],p1[i]))
  z<-r*x1+(1-r)*x2
  hist(z,prob=TRUE,ylim=c(0,.3),main=p1[i])   #plot the histograms
  lines(density(z))  #lines the density curve
}

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
m<-1e4
set.seed(1234)
# simple Monte Carlo
x<-runif(m,0,pi/3)
theta.hat<-mean(sin(x))*pi/3
print(c(theta.hat,1-cos(pi/3)))

## -----------------------------------------------------------------------------
m<-1e4
set.seed(1234)
x<-runif(m/2)
theta.hat<-mean(c(exp(-x)/(1+x^2),exp(x-1)/(1+(1-x)^2)))
print(theta.hat)

## -----------------------------------------------------------------------------
# compute variances
n<-1000
MC1<-MC2<-numeric(n)
for (i in 1:n) {
  x<-runif(m)
  y<-runif(m/2)
  MC1[i]<-mean(exp(-x)/(1+x^2))
  MC2[i]<-mean(c(exp(-y)/(1+y^2),exp(y-1)/(1+(1-y)^2)))
}
print(c(var(MC2)/var(MC1),(var(MC1)-var(MC2))/var(MC1)))

## -----------------------------------------------------------------------------
m<-1e4    # number of replicates
set.seed(1234) 
k<-5      # number of strata
r<-m/k    
N<-50     # number of times to repeat the estimation
T2<-numeric(k)
estimates<-matrix(0,N,2)
g<-function(x) {
  exp(-x)/(1+x^2)*(x>0)*(x<1)
}
for (i in 1:N) {
  u1<-runif(m)
  x<--log(1-u1*(1-exp(-1)))    # inverse transformated method
  fg1<-g(x)/(exp(-x)/(1-exp(-1)))
  estimates[i,1]<-mean(fg1)
  for (j in 1:k) {
    u<-c(0,0.2,0.4,0.6,0.8,1)
    quintile<--log(1-u*(1-exp(-1)))   # compute quintiles
    u2<-runif(m/5)
    y<--log(exp(-quintile[j])-u2/5*(1-exp(-1)))
    fg2<-g(y)/(5*exp(-y)/(1-exp(-1)))
    T2[j]<-mean(fg2)
  }
  estimates[i,2]<-sum(T2)
}
print(estimates)

apply(estimates,2,mean)
apply(estimates,2,var)
print((var(estimates[,1])-var(estimates[,2]))/var(estimates[,1]))

## -----------------------------------------------------------------------------
n<-20       #sample size
alpha<-.05  #confidence level
m<-1000     #replicates
LCL<-UCL<-numeric(m)
set.seed(1234)
#compute LCL and UCL for each replicate
for (i in 1:1000) {
  x<-rchisq(n,2)
  LCL[i]<-mean(x)-qt(1-alpha/2,n-1)*sd(x)/sqrt(n)
  UCL[i]<-mean(x)+qt(1-alpha/2,n-1)*sd(x)/sqrt(n)
}
#compute the probability that the confidence interval covers the mean
mean((LCL<2)*(UCL>2))

## ----message=FALSE, warning=FALSE---------------------------------------------
library(moments)  #to use function 'skewness' in this package
set.seed(1234)
n<-1000    #sample size
m<-1000
k<-100    #replicates
q<-c(0.025,0.05,0.95,0.975)   #four quantiles
quantile<-matrix(0,k,4)
#MC experiments to generate k replicates of estimate under normality
for (i in 1:k) {
  b1<-replicate(m,expr={
    x<-rnorm(n)
    skewness(x)
  })
  quantile[i,1]<-quantile(b1,q[1])
  quantile[i,2]<-quantile(b1,q[2])
  quantile[i,3]<-quantile(b1,q[3])
  quantile[i,4]<-quantile(b1,q[4])
}
quantile_sample<-c(c(mean(quantile[,1]),mean(quantile[,2]),mean(quantile[,3]),mean(quantile[,4])))
print(quantile_sample)

## -----------------------------------------------------------------------------
hist(b1)
qqnorm(b1)

## -----------------------------------------------------------------------------
se<-numeric(4)
#compute se using (2.14) formula
for(i in 1:4) {
  se[i]<-sqrt(q[i]*(1-q[i])/(m*dnorm(quantile_sample[i],0,sqrt(6/n))))
}
print(se)

## -----------------------------------------------------------------------------
#estimate quantiles using large sample approximation
quantile_large<-numeric(4)
for (i in 1:4) {
  quantile_large[i]<-qnorm(q[i],0,sqrt(6/n))
}
print(quantile_large)

knitr::kable(rbind(quantile_sample,quantile_large),col.names=c(0.025,0.05,0.95,0.975))

## -----------------------------------------------------------------------------
y<-rnorm(n,0,sqrt(6/n))
hist(y)

## -----------------------------------------------------------------------------
alpha<-0.05
set.seed(1234)
n<-30
m<-2500
parameter<-1:10
N<-length(parameter)
pwr<-numeric(N)
# critical value for the skewness test
cv<-qnorm(1-alpha/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))

# define sk
sk<-function(x) {
  #compute the sample skewness coeff
  xbar<-mean(x)
  m3<-mean((x-xbar)^3)
  m2<-mean((x-xbar)^2)
  return(m3/m2^1.5)
}


for (j in 1:N) {
  p<-parameter[j]
  sktests<-numeric(m)
  for (i in 1:m) {
    x<-rbeta(n,p,p)
    sktests[i]<-as.integer(abs(sk(x))>=cv)
  }
  pwr[j]<-mean(sktests)
}

#plot power vs parameter
plot(parameter,pwr,type='b')

## -----------------------------------------------------------------------------
for (j in 1:N) {
  p<-parameter[j]
  sktests<-numeric(m)
  for (i in 1:m) {
    y<-rt(n,p)
    sktests[i]<-as.integer(abs(sk(y))>=cv)
  }
  pwr[j]<-mean(sktests)
}
plot(parameter,pwr,type='b')

## -----------------------------------------------------------------------------
n<-c(10,20,30,50,100,500)
alpha<-.05
p.hat1<-p.hat2<-p.hat3<-numeric(length(n))

m<-10000    #number of replicates
p<-numeric(m)
for (i in 1:length(n)) {
  for (j in 1:m) {
    x<-rchisq(n[i],1)
    ttest<-t.test(x,mu=1)
    p[j]<-ttest$p.value
  }
  p.hat1[i]<-mean(p<alpha)
}
p.hat1


for (i in 1:length(n)) {
  for (j in 1:m) {
    x<-runif(n[i],0,2)
    ttest<-t.test(x,mu=1)
    p[j]<-ttest$p.value
  }
  p.hat2[i]<-mean(p<alpha)
}
p.hat2

for (i in 1:length(n)) {
  for (j in 1:m) {
    x<-rexp(n[i],1)
    ttest<-t.test(x,mu=1)
    p[j]<-ttest$p.value
  }
  p.hat3[i]<-mean(p<alpha)
}
p.hat3

## -----------------------------------------------------------------------------
library(bootstrap)
library(boot)
data<-as.matrix(scor)
plot(scor)   # plot scatter plots
cor(scor)    # compute the sample correlation matrix

## -----------------------------------------------------------------------------
# define function to compute correlation
statistic12<-function(x,i) cor(x[i,1],x[i,2])
statistic34<-function(x,i) cor(x[i,3],x[i,4])
statistic35<-function(x,i) cor(x[i,3],x[i,5])
statistic45<-function(x,i) cor(x[i,4],x[i,5])

boot(data,statistic=statistic12,R=2000)
boot(data,statistic=statistic34,R=2000)
boot(data,statistic=statistic35,R=2000)
boot(data,statistic=statistic45,R=2000)

## -----------------------------------------------------------------------------
# for normal populations
library(boot)
set.seed(1234)
# define a function to compute skewness using function "boot" 
sk<-function(x,i) mean((x[i]-mean(x[i]))^3)/mean((x[i]-mean(x[i]))^2)^1.5
n<-500     # number of replicates in Monte Carlo
m<-200     # sample size
norm_cil<-norm_cir<-basic_cil<-basic_cir<-perc_cil<-perc_cir<-numeric(n)
# Monte Carlo and bootstrap
for (i in 1:n) {
  x<-rnorm(m)
  boot.obj<-boot(x,statistic=sk,R=2000)
  ci<-boot.ci(boot.obj,type=c('norm','basic','perc'))
  norm_cil[i]<-ci$normal[2]
  norm_cir[i]<-ci$normal[3]
  basic_cil[i]<-ci$basic[4]
  basic_cir[i]<-ci$basic[5]
  perc_cil[i]<-ci$percent[4]
  perc_cir[i]<-ci$percent[5]
}
# compute coverage rates 
print(c(mean((norm_cil<0)*(norm_cir>0)),mean((basic_cil<0)*(basic_cir>0)),
        mean((perc_cil<0)*(perc_cir>0))))

## -----------------------------------------------------------------------------
# for chi-square(5) dsitributions
library(boot)
set.seed(1234)
sk<-function(x,i) mean((x[i]-mean(x[i]))^3)/mean((x[i]-mean(x[i]))^2)^1.5
n<-500      # number of replicates in Monte Carlo
m<-200      # sample size
norm_cil<-norm_cir<-basic_cil<-basic_cir<-perc_cil<-perc_cir<-numeric(n)
# Monte Carlo and bootstrap
for (i in 1:n) {
  x<-rchisq(m,5)
  boot.obj<-boot(x,statistic=sk,R=2000)
  ci<-boot.ci(boot.obj,type=c('norm','basic','perc'))
  norm_cil[i]<-ci$normal[2]
  norm_cir[i]<-ci$normal[3]
  basic_cil[i]<-ci$basic[4]
  basic_cir[i]<-ci$basic[5]
  perc_cil[i]<-ci$percent[4]
  perc_cir[i]<-ci$percent[5]
}
# compute coverage rates 
skewness<-sqrt(8/5)   # skewness of chi-square(5)
print(c(mean((norm_cil<skewness)*(norm_cir>skewness)),mean((basic_cil<skewness)*(basic_cir>skewness)),mean((perc_cil<skewness)*(perc_cir>skewness))))

## -----------------------------------------------------------------------------
library(boot)
data(scor,package='bootstrap')
eigen<-eigen(cov(scor))    #compute eigens of covariance matrix Sigma
theta.hat<-eigen$values[1]/sum(eigen$values)    #original theta.hat

n<-nrow(scor)
theta.jack<-numeric(n)
for (i in 1:n) {
  eigen.jack<-eigen(cov(scor[-i,]))
  theta.jack[i]<-eigen.jack$values[1]/sum(eigen.jack$values)
}
bias.jack<-(n-1)*(mean(theta.jack)-theta.hat)     #jackknife estimate of bias
se.jack<-sqrt((n-1)*mean((theta.jack-theta.hat)^2))  #jackkinfe estimate of se

c(original=theta.hat,bias.jack=bias.jack,se.jack=se.jack)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(DAAG)
attach(ironslag)
n<-length(magnetic)
e5<-numeric(n)

#n-fold cross validation
#fit models on leave-one-out samples
for (k in 1:n) {
  y<-magnetic[-k]
  x<-chemical[-k]
  
  J5<-lm(y~x+I(x^2)+I(x^3))
  yhat5<-J5$coef[1]+J5$coef[2]*chemical[k]+J5$coef[3]*chemical[k]^2+J5$coef[4]*chemical[k]^3
  e5[k]<-magnetic[k]-yhat5
}
# the mean of the squared prediction errors
c(Linear=19.55644,Quadratic=17.85248,Exponential=18.44188,Cubic=mean(e5^2))

## -----------------------------------------------------------------------------
L1<-lm(magnetic~chemical)
L2<-lm(magnetic~chemical+I(chemical^2))
L3<-lm(log(magnetic)~chemical)
L5<-lm(magnetic~chemical+I(chemical^2)+I(chemical^3))
#compare four adjusted R squares
c(Linear=summary(L1)$adj.r.squared,Quadradic=summary(L2)$adj.r.squared,Exponential=summary(L3)$adj.r.squared,Cubic=summary(L5)$adj.r.squared)

## -----------------------------------------------------------------------------
n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1

set.seed(1234)
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)

# compute the maximum numebers of extreme points m1 and m2
log(.025) / log(n1 / (n1 + n2))
log(.025) / log(n2 / (n1 + n2))

## -----------------------------------------------------------------------------
m1 <- 4
m2 <- 7

# original statistic
counttest <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer((outx > m1) | (outy > m2)))
}

R <- 9999
z <- c(x,y)
K <- n1 + n2
reps <- numeric(R)
t0 <- counttest(x,y)
for (i in 1:R) {
  k <- sample(K, size = n1, replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k]
  X <- x1 - mean(x1)
  Y <- y1 - mean(y1)
  reps[i] <- counttest(x1, y1)
}

# compute alphahat
alphahat <- mean(c(t0, reps) > t0)
print(alphahat)

## -----------------------------------------------------------------------------
# costruct ndCov2 function for permulation bootstrap
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

ndCov2 <- function(z, ix, dims) {
  p <- dims[1]
  q <- dims[2]
  d <- p + q
  x <- z[ , 1:p]
  y <- z[ix, -(1:p)]
  return(nrow(z) * dCov(x, y)^2)
}

library(boot)
library(Ball)
set.seed(1234)
m <- 100
p <- 2
n1 <- n2 <- 20
R <- 999
n <- n1 + n2
# compute p values using two methods
p.cor <- p.ball <- matrix(0, m ,2)
  
for (i in 1:m) {
  X <- matrix(rnorm(n1 * p), ncol = p)
  e <- matrix(rnorm(n2 * p), ncol = p)
  Y1 <- X/4 + e
  Y2 <- X/4 * e
  z1 <- cbind(X, Y1)
  z2 <- cbind(X, Y2)
  boot.obj1 <- boot(data = z1, statistic = ndCov2, R = R, sim = 'permutation', dims = c(2,2))
  tb1 <- c(boot.obj1$t0, boot.obj1$t)
  boot.obj2 <- boot(data = z2, statistic = ndCov2, R = R, sim = 'permutation', dims = c(2,2))
  tb2 <- c(boot.obj2$t0, boot.obj2$t)
  p.cor[i, 1] <- mean(tb1 >= tb1[1])
  p.cor[i, 2] <- mean(tb2 >= tb2[1])
  p.ball[i, 1] <- bd.test(X, Y1, R = 999, seed = i * 1234)$p.value
  p.ball[i, 2] <- bd.test(X, Y2, R = 999, seed = i * 1234)$p.value 
}

# compute powers in two models using two methods
alpha <- .1
pow1 <- colMeans(p.cor < alpha)
pow2 <- colMeans(p.ball < alpha)
pow1
pow2

## -----------------------------------------------------------------------------
data(data2)
n <- c(10, 20, 30, 50, 80, 100)
plot(n, data[ , 1], type = 'b')
plot(n, data[ , 2], type = 'b')
plot(n, data[ , 3], type = 'b')
plot(n, data[ , 4], type = 'b')

## -----------------------------------------------------------------------------
# write a funciton of random walk Metropolis sampler
rw.Metropolis <- function(sigma, x0, N) {
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

N <- 2000
sigma <- c(.05, .5, 2, 16)
set.seed(1234)

# apply to different variances
x0 <- 25
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)
rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)

print(c(rw1$k, rw2$k, rw3$k, rw4$k))
# compute the acceptance rates of different varainces
print(c(1-rw1$k/N,1-rw2$k/N,1-rw3$k/N,1-rw4$k/N))

## -----------------------------------------------------------------------------
#par(mfcol = c(2, 2)) # display 4 graphs together
refline <- c(log(2 * .025), -log(2 * (1 - .975)))
for (j in 1:4) {
  plot(rw[, j], type = 'l', xlab = bquote(sigma == .(round(sigma[j], 3))), ylab = 'X', ylim = range(rw[, j]))
  abline(h = refline)
}
#par(mfcol = c(1, 1)) # reset to default

a <- c(.05, seq(.1, .9, .1), .95)
# compute true quantiles of standard Laplace distribution
Q1 <- log(2 * a[1:5])
Q2 <- -log(2 * (1 - a[6:11]))
Q <- c(Q1, Q2)

mc <- rw[501:N, ]
Qrw <- apply(mc, 2, function(x) quantile(x, a))
knitr::kable(round(cbind(Q, Qrw), 3), col.names = c('TRUE', 'rw1', 'rw2', 'rw3', 'rw4'))



## -----------------------------------------------------------------------------
x <- seq(1, 100, .1)
isTRUE(log(exp(x)) == exp(log(x)))

## -----------------------------------------------------------------------------
isTRUE(all.equal(log(exp(x)),exp(log(x))))
isTRUE(all.equal(log(exp(x)),x))
isTRUE(all.equal(exp(log(x)),x))

## -----------------------------------------------------------------------------
# find the intervals in which the roots fall  
f <- function(a, k) 1-pt(sqrt(a^2*k/(k+1-a^2)), k)
k <- c(4, 25, 100, 1000)
#par(mfcol =  c(2, 2))
for (i in k) {
  g <- function(x) f(x, i-1)-f(x, i)
    a <- seq(0, sqrt(i), .1)
  plot(a, g(a), type = 'l', main = paste('k=',i))
  abline(h = 0)
}

## -----------------------------------------------------------------------------
# Exercise 11.4
f <- function(a, k) 1-pt(sqrt(a^2*k/(k+1-a^2)), k)
k <- c(4:25, 100, 500, 1000)
Ak <- numeric(length(k))
i <- 1
for (j in k) {
  g <- function(x) f(x, j-1)-f(x, j)
  Ak[i] <- uniroot(g, lower = 1, upper = 2)$root
  i <- i+1
}
knitr::kable(cbind(k,Ak))

## -----------------------------------------------------------------------------
# Exercise 11.5
f <- function(k) 2/sqrt(pi*k)*exp(lgamma((k+1)/2)-lgamma(k/2))
ck <- function(a, k) sqrt(a^2*k/(k+1-a^2))
g <- function(u, k) (1+u^2/k)^(-(k+1)/2)
k <- c(4:25, 100, 500, 1000)
root <- numeric(length(k))
i <- 1
for (j in k) {
  ff <- function(a) f(j)*integrate(function(u) {g(u, j)}, 0, ck(a, j))$value -       f(j-1)*integrate(function(u) {g(u, j-1)}, 0, ck(a, j-1))$value 
  root[i] <- uniroot(ff, lower = 1, upper = 2)$root
  i <- i+1
}
knitr::kable(cbind(k, Ak, root))

## ----echo=FALSE---------------------------------------------------------------
dat <- rbind(Genotype=c('AA','Aa','aa',''),
             Frequency=c('p^2','2pq','q^2',1),
             Count=c('n_AA','n_Aa','n_aa','n'))
knitr::kable(dat)

## -----------------------------------------------------------------------------
nA <- 28
nB <- 24
nOO <- 41
nAB <- 70
# EM algorithm
theta0 <- c(.5, .3)
l <- numeric(1000)
for (j in 1:1000) {
  E <- function(theta) {
    p <- theta[1]
    q <- theta[2]
    r <- 1-p-q
    p0 <- theta0[1]
    q0 <- theta0[2]
    r0 <- 1-p0-q0
    return(2*nA*(log(p)+r0/(p0+2*r0)*log(2*r/p))+2*nB*(log(q)+r0/(q0+2*r0)*log(2*r/q))+2*nOO*log(r)+nAB*log(2*p*q))
  }
  Obj <- function(theta) -E(theta)
  optim <- optim(c(.1, .1), Obj)
  theta0 <- optim$par
  l[j] <- E(theta0)
}
print(theta0)
plot(l[1:10], type = 'l', xlab = 'iterations', ylab = 'log-maximum likelihood values')

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

## -----------------------------------------------------------------------------
attach(mtcars)
# with lapply
lapply(formulas, lm)
# for loop
out1 <- vector('list', length(formulas))
for (i in 1:length(formulas)) {
  out1[[i]] <- lm(formulas[[i]])
}
out1

## -----------------------------------------------------------------------------
set.seed(1234)
bootstrap <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

## ----message=FALSE, warning=FALSE---------------------------------------------
# with lapply
lapply(bootstrap, function(x) lm(mpg ~ disp, data = x))
# for loop
out2 <- vector('list', length(bootstrap))
for (i in 1:length(bootstrap)) {
  out2[[i]] <- lm(mpg ~ disp, data = bootstrap[[i]])
}
out2
# lapply without anonymous function
lapply(bootstrap, lm, formula = mpg ~ disp)

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

## -----------------------------------------------------------------------------
lapply(out1, rsq)
lapply(out2, rsq)

## -----------------------------------------------------------------------------
set.seed(1234)
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

## -----------------------------------------------------------------------------
# with sapply
sapply(trials, function(x) x$p.value)
# using [[
sapply(trials, '[[', 'p.value')

## ----eval=FALSE---------------------------------------------------------------
#  # mcsapply()
#  library(parallel)
#  boot_df <- function(x) x[sample(nrow(x), rep = T), ]
#  rsquared <- function(mod) summary(mod)$r.squared
#  boot_lm <- function(i) {
#    rsquared(lm(mpg ~ wt + disp, data = boot_df(mtcars)))
#  }
#  system.time(sapply(1:1e5, boot_lm))
#  system.time(unlist(mclapply(1:1e5, boot_lm, mc.cores = 4)))

## -----------------------------------------------------------------------------
library(SC19062)
library(microbenchmark)
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
rwC1 <- rwMetropolisC(sigma[1], x0, N)
rwC2 <- rwMetropolisC(sigma[2], x0, N)
rwC3 <- rwMetropolisC(sigma[3], x0, N)
rwC4 <- rwMetropolisC(sigma[4], x0, N)
rwC <- cbind(rwC1[-(N+1)], rwC2[-(N+1)], rwC3[-(N+1)], rwC4[-(N+1)])

print(c(rwC1[N+1], rwC2[N+1], rwC3[N+1], rwC4[N+1]))
# compute the acceptance rates of different varainces
print(c(1-rwC1[N+1]/N,1-rwC2[N+1]/N,1-rwC3[N+1]/N,1-rwC4[N+1]/N))

#par(mfcol = c(2, 2)) # display 4 graphs together
refline <- c(log(2 * .025), -log(2 * (1 - .975)))
for (j in 1:4) {
  plot(rwC[, j], type = 'l', xlab = bquote(sigma == .(round(sigma[j], 3))), ylab = 'X', ylim = range(rwC[, j]))
  abline(h = refline)
}
#par(mfcol = c(1, 1)) # reset to default

## -----------------------------------------------------------------------------
rw1 <- rwMetropolisR(sigma[1], x0, N)
rw2 <- rwMetropolisR(sigma[2], x0, N)
rw3 <- rwMetropolisR(sigma[3], x0, N)
rw4 <- rwMetropolisR(sigma[4], x0, N)

print(c(rw1$k, rw2$k, rw3$k, rw4$k))
# compute the acceptance rates of different varainces
print(c(1-rw1$k/N,1-rw2$k/N,1-rw3$k/N,1-rw4$k/N))

qqplot(rw1$x, rwC1[-(N+1)], xlab = 'rw1', ylab = 'rwC1')
qqplot(rw2$x, rwC2[-(N+1)], xlab = 'rw2', ylab = 'rwC2')
qqplot(rw3$x, rwC3[-(N+1)], xlab = 'rw3', ylab = 'rwC3')
qqplot(rw4$x, rwC4[-(N+1)], xlab = 'rw4', ylab = 'rwC4')



## -----------------------------------------------------------------------------
ts <- microbenchmark(rwC1 = rwMetropolisC(sigma[1], x0, N), rwC2 = rwMetropolisC(sigma[2], x0, N), rwC3 = rwMetropolisC(sigma[3], x0, N), rwC4 = rwMetropolisC(sigma[4], x0, N), rw1 = rwMetropolisR(sigma[1], x0, N), rw2 = rwMetropolisR(sigma[2], x0, N), rw3 = rwMetropolisR(sigma[3], x0, N), rw4 = rwMetropolisR(sigma[4], x0, N))
knitr::kable(summary(ts)[, c(1,3,5,6)])

