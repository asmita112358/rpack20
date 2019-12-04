library(MASS)
library(mvtnorm)
library(reticulate)
library(CVXR)
library(kernlab)

phi = function(theta, x)
{
  d = length(x)
  w = theta[1:d]
  b = theta[d+1]
  return( sqrt(2)*cos(sum(w*x)+2*pi*b) )
}

psi = function(VX, VY)
{ 
  n = length(VX)
  T1 = (sum(VX)^2 - sum(VX^2))/n/(n-1) 
  T2 = (sum(VY)^2 - sum(VY^2))/n/(n-1) 
  T3 = mean(VX)*mean(VY)
  out = T1 + T2 - 2*T3
  return( out )
}

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
optkern = function(X, Y, N, alfa = 0.05)
{
  #Generating theta, the base distribution for the monte carlo estimation
  Hi = 0.35
  Lo = -Hi
  w = matrix(runif(N*d, min=Lo, max=Hi), N, d)
  theta = cbind(w, runif(N))
  q = rep(1/2/Hi, N)
  
  stat = 0
  n = length(X) 
  ps = rep(NA, N)
  psiX = matrix(NA, n, N) 
  psiY = matrix(NA, n, N)
  for(i in 1:N)
  {
    #print("got in loop")
    psiX[,i] = apply(X, 1, phi, theta = theta[i,])
    #print(psiX[,i])
    psiY[,i] = apply(Y, 1, phi, theta = theta[i,])
    ps[i] = psi(psiX[,i], psiY[,i])/q[i]
  }
  dpsiX = t(t(psiX) - colMeans(psiX))
  dpsiY = t(t(psiY) - colMeans(psiY))
  Sig = (crossprod(dpsiX) + crossprod(dpsiY))/(n+n-2)
  Qinv = diag(1/q) 
  
  C = 2*(1/n/(n-1) + 1/n/(n-1) + 2/n/n)
  B = C*Qinv%*%(Sig*Sig)%*%Qinv
  eB = eigen(B)
  ev = eB$values
  ev[ev<0] =  min(ev[ev>0])
  Bsi = eB$vectors%*%diag(1/sqrt(ev))%*%t(eB$vectors)
  stat = test(theta, Bsi, ps)
  
  data = rbind(psiX, psiY)
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
  crit = quantile(stat.per, 1 - alfa)
  d = stat > crit
  return(list(teststat = stat, criticalvalue = crit, decision = d))
}
