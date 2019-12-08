
packages <- c("MASS", "mvtnorm", "reticulate", "CVXR", "kernlab")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

library(MASS)
library(mvtnorm)
library(reticulate)
library(CVXR)
library(kernlab)

## This function calculates the kernel
phi = function(theta, x)
{
  d = length(x)
  w = theta[1:d]
  b = theta[d+1]
  return( sqrt(2)*cos(sum(w*x)+2*pi*b) )
}

## This function calculates the U statistic, which is based on the Kernel function phi
psi = function(VX, VY)
{ 
  n = length(VX)
  T1 = (sum(VX)^2 - sum(VX^2))/n/(n-1) 
  T2 = (sum(VY)^2 - sum(VY^2))/n/(n-1) 
  T3 = mean(VX)*mean(VY)
  out = T1 + T2 - 2*T3
  return( out )
}
##computes the test statistics based on on inputs which are some functions of the data
##This function's main job is to perform the optimization that obtains the optimal kernel.
#Dr. Zhang asked me to use the package CVXR for this particular optimization
test = function(theta, Bsi, ps)
{
  N = dim(theta)[1]
  pst = Bsi%*%ps
  p = CVXR::Variable(N)
  obj.p = sum(pst*p)           #pst = the shi-tilde vector in last step
  obj.n = -sum(pst*p)
  
  constr = list(sum(p^2) <= 1, Bsi%*%p >= 0)
  prob.p = CVXR::Problem(CVXR::Maximize(obj.p), constr)
  res.p = CVXR::psolve(prob.p, solver = "ECOS")
  prob.n = CVXR::Problem(CVXR::Maximize(obj.n), constr)
  res.n = CVXR::psolve(prob.n, solver="ECOS")
  val = max(res.n$value, res.p$value)
  return(val)
}

#X- an nxd data matrix containing n vectors from the distribution of X, each of length d
#Y - an nxd data matrix containing n vectors from the distribution of Y, each of length d
#aim - to chekc whether the two random samples are from the same distribution
# N - a parameter to be chosen by the user, ideally close to half of n.
#alfa - significance level, automatically taken as 0.05
#' Title
#'
#' @param X an nxd data matrix containing n vectors from the distribution of X, each of length d
#' @param Y an nxd data matrix containing n vectors from the distribution of Y, each of length d
#' @param N a parameter to be chosen by the user, ideally close to half of n.
#' @param alfa significance level, automatically taken as 0.05
#'
#' @return optimal kernel test statistic, the critical point and the decision whether the test is rejected
#' @export
#'
#' @examples
#' Z = optkern(X,Y, alpha = 0.01).The null hypothesis is the that the distribution of X and Y are equal.
optkern = function(X, Y, N, alfa = 0.05)
{
  #Generating theta, the base distribution for the monte carlo estimation. 
  #Here we take the distribution of theta as uniform. 
  #Further research is ongoing about whether there's any way to fine tune this choice.
  Hi = 0.35
  Lo = -Hi
  w = matrix(runif(N*d, min=Lo, max=Hi), N, d)
  theta = cbind(w, runif(N))
  q = rep(1/2/Hi, N)
  
  ##Performing checks on the data:
  if(dim(X) != dim(Y))
    stop(print("X and Y should have same number of rows and columns"))
  if(nrow(X)<ncol(X))
    stop(print("Number of samples must be larger than the length of each sample"))
  if(alfa >1 || alfa <0)
    stop(print("size of the test must lie between zero and one"))
  
  stat = 0
  n = length(X) 
  ps = rep(NA, N)
  psiX = matrix(NA, n, N) 
  psiY = matrix(NA, n, N)
  #in this loop we apply the kernel function on the data points.
  for(i in 1:N)
  {
    #print("got in loop")
    psiX[,i] = apply(X, 1, phi, theta = theta[i,])
    #print(psiX[,i])
    psiY[,i] = apply(Y, 1, phi, theta = theta[i,])
    ps[i] = psi(psiX[,i], psiY[,i])/q[i]
  }
  #centering the kernelled data
  dpsiX = t(t(psiX) - colMeans(psiX))
  dpsiY = t(t(psiY) - colMeans(psiY))
  #Estimating the Sigma Matrix
  Sig = (crossprod(dpsiX) + crossprod(dpsiY))/(n+n-2)
  Qinv = diag(1/q) 
  
  #Construct the B matrix in the process, perform spectral analysis of B, to be used in later steps
  C = 2*(1/n/(n-1) + 1/n/(n-1) + 2/n/n)
  B = C*Qinv%*%(Sig*Sig)%*%Qinv
  eB = eigen(B)
  ev = eB$values
  ev[ev<0] =  min(ev[ev>0])
  Bsi = eB$vectors%*%diag(1/sqrt(ev))%*%t(eB$vectors)
  #Call the test function with entries 
  #1.theta, which is a random sample from the base distribution 
  #2. Bsi, ps: matrices and vectors computed halfway through the process
  #I'm sorry if this part is difficult to understand, the theory is way too complicated here and the steps don't have any name because the research procedure is still going on.
  stat = test(theta, Bsi, ps)
  
  data = rbind(psiX, psiY)
  #performing permutation test with the available data with 100 replications.

  
  stat.per = rep(NA, 100)
  for(ii in 1:100)
  {
    data.per = data[sample(1:(2*n)), ]
    psiX.per = data.per[1:n,]
    psiY.per = data.per[(n+1):(2*n),]   
    ps.per = rep(NA, N)
    for(jj in 1:N)
    {
      ps.per[jj] = psi(psiX.per[,jj], psiY.per[,jj])/q[jj]
    }
    # tryCatch({stat.per[ii] = test(theta, Bsi, ps.per)}, error=function(e){} )
    
    stat.per[ii] = test(theta, Bsi, ps.per)
  } 
  
  stat.per = stat.per[!is.na(stat.per)]
  
  #obtaining the critical value by selecting the (1-alpha)-th quantile.
  crit = quantile(stat.per, 1 - alfa)
  #if the test statistics is greater than the critical value the null hypothesis is rejected.
  d = stat > crit
  return(list(teststat = stat, criticalvalue = crit, decision = d))
}
