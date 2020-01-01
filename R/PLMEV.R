#' @title PLMEV procedure
#' @description a robust estimating method in semiparametric partially linear errors-in-variables model
#' @param y observed response Y
#' @param vv  observed data V
#' @param xt  observed t in nonparametric components
#' @param palabra selects the robust estimators of multivariate location and scatter, "mcd": MCD- (minimum covariance determinant) estimators, "m": M- estimators, "s": S- estimators, default is "mcd" 
#' @param cw quantile of a chi^2, default is 0.975
#' @param k tuning constant of score function, default 0.55 (90% efficiency under normality)
#' @param ventana smoothing parameter
#' @return estimates
#' @examples
#' \dontrun{
#' n=200
#' z1 <- rnorm(n)
#' z2 <- rnorm(n)
#' eps <- rnorm(n,0,0.2)
#' eps1 <- rnorm(n,0,0.2)
#' eps2 <- rnorm(n,0,0.2)
#' xx <- cbind(z1,z2)
#' xt <- runif(n,-6,6)
#' vv <- cbind(z1+eps1,z2+eps2)
#' y <- 2*xx[,1] + xx[,2] + cos(pi*xt/6) + eps
#' gg <- cos(pi*xt/6)
#' ventana=2.5
#' salida=PLMEV.rob.2(y,vv,xt,palabra="mcd",cw=0.975,k=0.55,ventana)
#' }
#' @export
#-----------------------------------------------------------------
# 3-Stepwise Procedure
#-----------------------------------------------------------------
PLMEV.rob.2<- function(y,vv,xt,palabra="mcd",cw=0.975,k=0.55,ventana)
{
  library(robustbase)
  library(isotone)
  library(MASS)
  library(rrcov)
  n=length(y)
  p=dim(vv)[2]
  #-----------------------------------------------------------------
  #Step 1: computation of local M-estimators
  #-----------------------------------------------------------------
  g0.m = Local.M(y,xt,k,ventana)
  gv.m = matrix(rep(0,n*p),n,p)
  for (j in 1:p)
  {
    gv.m[,j] = Local.M(vv[,j],xt,k,ventana)
  }
  #-----------------------------------------------------------------
  #Step 2: estimation of regression parameter
  # The weighted orthogonal estimator can be replaced by any
  # suitable estimator
  #-----------------------------------------------------------------
  v.tilde.m = vv - gv.m
  y.tilde.m = y - g0.m
  Z.m = cbind(v.tilde.m,y.tilde.m)
  estimador.beta= fekgaz(Z.m,palabra,cw)
  #-----------------------------------------------------------------
  #Step 3: estimation of the regression function g
  #-----------------------------------------------------------------
  g.rob<- g0.m - estimador.beta %*% t(gv.m)
  list(beta.rob=estimador.beta,nu0=g0.m,nus= gv.m, ghat=g.rob)
}
#-----------------------------------------------------------------
#Local M-estimation
#-----------------------------------------------------------------
## Local M-estimator with normal kernel
##
## Input:
## y: input argument
## t: tuning constant
## ventana= smoothing parameter
## the kernels scaled so that their quartiles
## (viewed as probability densities) are at +/- 0.25*ventana
## k: tuning constant of score function
## default 0.55 (90% efficiency under normality)
## 0.92 (95% efficiency under normality)
## Output:
## Local M-estimator
##
Local.M<- function(y,t,k=0.55,ventana)
{
  n=length(y)
  g0.m = rep(0,n)
  g0.mediano = 0
  bandwidth= 0.25*ventana/0.675
  for (i in 1: n)
  {
    we=dnorm(t - t[i],0,bandwidth)
    we1=we/sum(we)
    g0.mediano= weighted.median(y, we1)
    saliday = rreg3(y,g0.mediano,we,k)
    g0.m[i] = saliday
  }
  Local.M=g0.m
  Local.M
}
## Weight function related to score function psi_1 arctan for local M-estimator
##
## Input:
## u: argument
## k: tuning constant
## Output:
## psi_1(u)/u
wtan<- function(u,k)
{w=u/k
v=ifelse(abs(w)<0.001,1,atan(w)/w)
v= v
}
## Computation of local M-estimator
##
## Input:
## y: observations
## init: initial value (usually local median)
## k: tuning constant
## default 0.55 (90% efficiency under normality)
## 0.92 (95% efficiency under normality)
## Output:
## psi_1(u)/u
rreg3<-function(y,init,wx,k=0.55)
{
  estfin<-init
  scale <- weighted.median(abs(y-init), wx)
  if(scale==0) stop("bandwidth too small")
  for (i in 1:5)
  {
    estinicial<-estfin
    res<-y-estinicial
    w<- wtan(res/scale, k)*wx
    num<-sum(w*y,na.rm=T)
    den<-sum(w,na.rm=T)
    if (den<10^(-10)) {estfin<-.Machine$double.xmax}
    if (den>=10^(-10)) {estfin<-num/den}
  }
  estfin
}
#-----------------------------------------------------------------
#Functions needed for Fekri-Ruiz Gazen (2004) estimator
#-----------------------------------------------------------------
## Weight Function
##
## Hard rejection function on the interval I:[0,c] c is the cw quantile of a
## chi^2 distribution with n degrees
## Input:
## x: argument
## k: degrees of freedom (dim(beta)+1)
## cw: quantile of a chi^2, default is 0.975
## Output:
## Hard rejection function on the interval I:(0,c) at point x
w1<-function(x,k,cw=0.975)
{
  c<-qchisq(cw,k)
  aux1 <- 1*(x >= 0)
  aux2 <- 1*(x <= c)
  y <- aux1*aux2
  y
}
## Fekri-Ruiz Gazen (2004) estimator
## Model Y = X beta
## Input:
## Z: n x (p+1) matrix (First p columns: x, last column: Y)
## metodo: selects the robust estimators of multivariate location and scatter
## "mcd": MCD- (minimum covariance determinant) estimators
## "m": M- estimators
## "s": S- estimators
## default is "mcd"
## cw: quantile of a chi^2, default is 0.975
## Output:
## weighted orthogonal esstimator of parameter beta
fekgaz<-function(Z,metodo="mcd",cw=0.975)
{
  ##Computation of initial multivariate estimates of mu and Sigma
  #-----------------------------------------------------------------
  if(metodo=="mcd")
  {
    fekgaz <- CovMcd(Z)
    mu.n <- fekgaz@center
    sigma.n <- fekgaz@cov
  }
  if(metodo=="m")
  {
    fekgaz <- CovMest(Z)
    mu.n <- fekgaz@center
    sigma.n <- fekgaz@cov
  }
  if(metodo=="s"){
    fekgaz <- CovSest(Z)
    mu.n <- fekgaz@center
    sigma.n <- fekgaz@cov
  }
  #Computation of one-step reweighted versions of mu and Sigma
  #-----------------------------------------------------------------
  n <- length(Z[,1])
  q <- length(Z[1,])
  numerador.mu <- 0
  denominador.mu <- 0
  numerador.sigma <- 0
  ## Estimator of mu
  auxiliar <- rep(mu.n,n)
  dim(auxiliar) <- c(q,n)
  mu.n.matriz <- t(auxiliar)
  Z.menos.mu.n <- Z - mu.n.matriz
  sigma.n.inv<- solve(sigma.n)
  auxi <- diag(Z.menos.mu.n %*% sigma.n.inv %*% t(Z.menos.mu.n))
  dim(auxi)<-c(n,1)
  pesos <- w1(auxi,q,cw)
  denominador <- sum(pesos)
  dim(pesos)<-NULL
  terminos.mu <- Z * pesos
  numerador.mu <- colSums(terminos.mu)
  mu.r.n <- numerador.mu/denominador
  ## Estimator of Sigma
  numerador.sigma <- 0
  auxiliar1 <- rep(mu.r.n,n)
  dim(auxiliar1) <- c(q,n)
  mu.r.n.matriz <- t(auxiliar1)
  Z.menos.mu.r.n <- Z - mu.r.n.matriz
  auxi2 <- diag(sqrt(pesos)) %*% (Z.menos.mu.r.n)
  numerador.sigma <- t(auxi2) %*% (auxi2)
  sigma.r.n <- numerador.sigma/denominador
  #Final estimators of beta
  ###########################
  p <- length(Z[1,]) - 1
  menor.autovec <- eigen(sigma.r.n)$vectors[,p+1]
  beta.n <- -(1/menor.autovec[p+1]) * menor.autovec[1:p]
  beta.n
}