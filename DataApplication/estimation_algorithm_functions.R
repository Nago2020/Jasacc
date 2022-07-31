
#==========================================================
# Copula based Cox proportional hazards
# models for dependent censoring
# [Authors and affiliation blinded]
# July, 2022
#==========================================================



# load libraries 

library(copula)
library(survival)
library(pbivnorm)



#==========================================================
# Estimation algorithm
#==========================================================

# This file contains estimation algorithm for the 
# proposed model in the paper



# Estimation of the cumulative hazard function, Lambda
# SolveL is used to estimate Lambda recursively for
# a semiparametric copula model. 
# The arguments of SolveL:
                        # cop = type of copula
                        # theta = parameter estimates at a given iteration  
                        # dist = distribution of dependent censoring, C
                        # Z = observed survival time
                        # d1 = censoring indicator for T
                        # d2 = censoring indicator for C
                        # X, W =  data matrices

SolveL = function(theta,Z,d1,d2,X,W,cop,dist){
  T1 <- c(0,sort(unique(Z[d1 == 1])))
  m = length(T1)
  k = dim(X)[2]
  beta = theta[1:k]
  
  L = rep(0,m)
  L[1] = 0
  L[2] = sum((Z==T1[2]))/sum((Z>=T1[2])*exp(X%*%beta))
  
  for (i in 3:m){
    csum = sum(L[1:(i-1)])
    psi = CompC(theta,T1[i],X,W,csum,cop,dist)         # used to compute psi value in Section 4 of the paper
    L[i] <- sum(Z == T1[i])/sum((Z>=T1[i])*exp(psi))
  }
  res <- list(lambda = L,cumhaz = cumsum(L), times = T1)
}


# SolveLI - estimates Lambda under independent
# censoring assumption
# The arguments of SolveLI:
                        # theta = parameter estimates at a given iteration  
                        # Z = observed survival time
                        # d1 = censoring indicator for T
                        # d2 = censoring indicator for C
                        # X, W =  data matrices

SolveLI = function(theta,Z,d1,d2,X){
  T1 <- c(sort(unique(Z[d1 == 1])))
  m = length(T1)
  k = dim(X)[2]
  beta = theta[1:k]
  
  L = rep(0,m)
  for (i in 1:m){
    L[i] <- sum(Z == T1[i])/sum((Z>=T1[i])*exp(X%*%beta))
  }
  res <- list(lambda = L,cumhaz = cumsum(L), times = T1)
}


# CompC - support function for SolveL: 
# The arguments of CompC:
                        # theta = parameter estimates at a given iteration  
                        # t =  a fixed time point
                        # ld = cumulative hazard value at previous time point, t-1
                        # cop = type of copula
                        # X,W =  data matrices
                        # dist = distribution of dependent censoring, C

CompC = function(theta,t,X,W,ld,cop,dist){
  k = dim(X)[2]
  l = dim(W)[2]
  beta = theta[1:k]
  eta = theta[(k+1):(l+k)]
  nu = theta[(l+k+1)]
  gm = theta[(l+k+2)]
  lfun = (W%*%eta)
  y = log(t)
  
  # distribution of T
  G1 = 1-exp(-ld*exp(X%*%beta))       # Cox model
  
  # distribution of C
  
  if (dist == "wb"){                  # Weibull margin
    G2 = 1-exp(-exp((y-W%*%eta)/nu))
  }
  if (dist == "lgn"){                 # Log-normal margin
    m = (y-W%*%eta)/nu
    G2 = pnorm(m)
  }
  if (dist == "gomp"){                # Gompertz margin
    Y = Z
    mu2 = exp(W%*%eta)
    G2 = 1-exp(exp(-mu2/nu)*(1-exp(Y/nu)))
  }
  
  # avoid numerical issues
  G1[is.nan(G1)] = 1e-10
  G2[is.nan(G2)] = 1e-10
  
  G1[is.na(G1)] <- 1e-10
  G2[is.na(G2)] <- 1e-10
  G1 = pmax(1e-6,G1)
  G2 = pmax(1e-6,G2)
  G1 = pmin(G1,0.999999)
  G2 = pmin(G2,0.999999)
  
  # Joint distribution 
  
  PD = cbind(G1,G2)
  
  if (cop == 1)                      # Clayton copula
  {  
    clay.cop = claytonCopula(gm, dim = 2)
    cp1 = (G1^(-gm)+G2^(-gm)-1)^(-1/gm-1)*G1^(-gm-1)
    cp1 = pmin(cp1,0.99999)
    z6 =  1-G1-G2+pCopula(PD,clay.cop)                 
    z6 =   pmax(z6,1e-10)
  }
  
  if (cop == 2)                      # Frank copula
  {
    frank.cop = frankCopula(gm, dim = 2)
    cp1 = (exp(- gm*G1)*(exp(- gm*G2)-1))/(exp(-gm)-1 + (exp(-gm*G1)-1)*(exp(-gm*G2)-1))
    cp1 = pmax(1e-10,cp1)
    z6 =  1-G1-G2+pCopula(PD,frank.cop)                    
    z6 =   pmax(z6,1e-10)
  }
  if (cop == 3)                    # Gumbel copula
  {
    gumb.cop = gumbelCopula(gm, dim = 2)
    cp1 = pCopula(PD,gumb.cop)*((-log(G1))^gm+(-log(G2))^gm)^(-1+1/gm)*(-log(G1))^(gm-1)/G1
    cp1[is.na(cp1)] <- 1e-5
    cp1 = pmax(1e-5,cp1)
    z6 =  1-G1-G2+pCopula(PD,gumb.cop)              
    z6 =   pmax(z6,1e-10)
  }
  if (cop == 4)              # normal 
  {
    cp1 = pnorm((qnorm(G2)-gm * qnorm(G1))/sqrt(1 - gm^2))
    p1 = qnorm(G1)
    p2 = qnorm(G2)
    c2 = cbind(p1,p2)
    
    z6 =  1-G1-G2+pbivnorm(c2,rho = gm)       # contribution to the likelihood from A
    z6 =   pmax(z6,1e-20)
  }
  B0 = X%*%beta-ld*exp(X%*%beta)
  B1 = log(pmax(1e-10,1-cp1))-log(z6)
  tot = B0+B1
  return(tot)
}

SearchIndicate = function(t,T1){       # Indicate position of t in T1
  i = 1;
  m = length(T1);
  
  while(T1[i+1]<=t){
    i = i+1;
    if (i == m) break; 
  }
  return(i)
}

# Change hazard and cumulative hazard to 
# long format

Longfun = function(Z,T1,lhat,Lhat){ 
  n = length(Z)
  llong = rep(0,n)
  Llong = rep(0,n)
  
  for(i in 1:n){
    j = SearchIndicate(Z[i],T1)
    llong[i] = lhat[j]
    Llong[i] = Lhat[j]
  }
  long = cbind(llong,Llong)
}

Distance = function(b,a){   # L2-norm of difference
  x = b-a
  n = length(x)
  l2norm = sqrt(sum(x^2)/n)
  return(l2norm)
}

# to obtain quantile from survival function

get_surv_prob <- function(fit, times) {
  stepfun(fit$time, c(1, fit$surv))(times)
}


# Pseudo-likelihood function
# We maximize PseudoL function in order to
# estimate the finite dimensional model parameters,
# including the dependency parameter

# The arguments of PseudoL:
                      # theta = parameter estimates at a given iteration  
                      # Z =  observed survival time
                      # d1 = censoring indicator for T
                      # d2 = censoring indicator for C
                      # lhat = estimated hazard function
                      # cumL = estimated cumulative hazard function
                      # cop = type of copula
                      # X,W =  data matrices
                      # dist = distribution of dependent censoring, C

PseudoL = function(theta,Z,d1,d2,X,W,lhat,cumL,cop,dist){
  
  k = dim(X)[2]
  l = dim(W)[2]
  beta = theta[1:k]
  eta = theta[(k+1):(l+k)]
  nu = theta[(l+k+1)]
  gm = theta[(l+k+2)]
  y = log(Z)
  
  # distribution of T -- Cox PH
  
  g1 = lhat*exp(X%*%beta)*exp(-cumL*exp(X%*%beta))
  G1 = 1-exp(-cumL*exp(X%*%beta))
  
  # dist of C
  
  if (dist == "wb"){    # Weibull
    g2 = 1/nu*exp((y-W%*%eta)/nu)*exp(-exp((y-W%*%eta)/nu))      #log-time scale
    G2 = 1-exp(-exp((y-W%*%eta)/nu))
  }
  if (dist == "lgn"){   # Lognormal
    m = (y-W%*%eta)/nu
    g2 = 1/(sqrt(2*pi)*nu*Z)*exp(-(m^2/2))
    G2 = pnorm(m)
  }
  if (dist == "gomp"){
    Y = Z
    mu2 = exp(W%*%eta)
    g2 = 1/nu*exp((Y-mu2)/nu)*exp(exp(-mu2/nu)*(1-exp(Y/nu)))
    G2 = 1-exp(exp(-mu2/nu)*(1-exp(Y/nu)))
  }
  
  # Joint distribution
  # avoid numerical issues
  
  g1[is.na(g1)] = 1e-20
  g2[is.na(g2)] = 1e-20
  g1[is.nan(g1)] = 1e-20
  g2[is.nan(g2)] = 1e-20
  G1[is.na(G1)] = 1e-10
  G2[is.na(G2)] = 1e-10
  G1[is.nan(G1)] = 1e-10
  G2[is.nan(G2)] = 1e-10
  g1 = pmax(g1,1e-30)
  g2 = pmax(g2,1e-30)
  g1[!is.finite(g1)] = 0
  g2[!is.finite(g2)] = 0
  
  G1 = pmax(G1,1e-10)
  G1 = pmin(G1,0.999999)
  G2 = pmax(G2,1e-10)
  G2 = pmin(G2,0.999999)
  
  PD = cbind(G1,G2)
  
  if (cop == 1)                                          # Clayton copula
  {  
    clay.cop = claytonCopula(gm, dim = 2)
    cp1 = (G1^(-gm)+G2^(-gm)-1)^(-1/gm-1)*G1^(-gm-1)
    cp2 = (G1^(-gm)+G2^(-gm)-1)^(-1/gm-1)*G2^(-gm-1)
    z6 =  1-G1-G2+pCopula(PD,clay.cop)                   # Joint survival function       
    z6 =   pmax(z6,1e-10)
  }
  if (cop == 2)                                          # Frank copula
  {
    frank.cop = frankCopula(gm, dim = 2)
    cp1 = (exp(- gm*G1)*(exp(- gm*G2)-1))/(exp(-gm)-1 + (exp(-gm*G1)-1)*(exp(-gm*G2)-1))
    cp2 = (exp(- gm*G2)*(exp(- gm*G1)-1))/(exp(-gm)-1 + (exp(-gm*G1)-1)*(exp(-gm*G2)-1))
    
    z6 =  1-G1-G2+pCopula(PD,frank.cop)                    
    z6 =   pmax(z6,1e-10)
  }
  if (cop == 3)                                       # Gumbel copula
  {  
    gumb.cop = gumbelCopula(gm, dim = 2)
    cp1 = pCopula(PD,gumb.cop)*((-log(G1))^gm+(-log(G2))^gm)^(-1+1/gm)*(-log(G1))^(gm-1)/G1
    cp2 = pCopula(PD,gumb.cop)*((-log(G1))^gm+(-log(G2))^gm)^(-1+1/gm)*(-log(G2))^(gm-1)/G2
    cp1[is.na(cp1)] <- 0.00001
    cp2[is.na(cp2)] <- 0.00001
    
    z6 =  1-G1-G2+pCopula(PD,gumb.cop)               
    z6 =   pmax(z6,1e-10)
  }
  if (cop == 4)     # normal 
  {
    cp1 = pnorm((qnorm(G2) - gm * qnorm(G1))/sqrt(1 - gm^2))
    cp2 = pnorm((qnorm(G1) - gm * qnorm(G2))/sqrt(1 - gm^2))
    p1 = qnorm(G1)
    p2 = qnorm(G2)
    c2 = cbind(p1,p2)
    z6 =  1-G1-G2+pbivnorm(c2,rho = gm)       # contribution to the likelihood from A
    z6 =   pmax(z6,1e-10)
  }
  cp1[is.na(cp1)] = 1e-20
  cp2[is.na(cp2)] = 1e-20
  cp1[is.nan(cp1)] = 1e-20
  cp2[is.nan(cp2)] = 1e-20
  cp1 = pmin(cp1,0.999999)
  cp2 = pmin(cp2,0.999999)
  cp1 = pmax(cp1,1e-20)
  cp2 = pmax(cp2,1e-20)
  cp1[!is.finite(cp1)] = 0
  cp2[!is.finite(cp2)] = 0
  d3 =  (1-d1-d2)             # Censoring
  term1 <- (1 - cp1)*g1       #  Contribution from T
  term2 <- (1 - cp2)*g2       #  Contribution from C
  
  Logn <- -sum(d1 * log(pmax(term1, 1e-10)) + d2 * log(pmax(term2, 1e-10)) + d3*log(z6))
  return(Logn)
}


# NEIN - main function in implementing the proposed estimation algorithm
# see the Supplementary Material for details on this estimation algorithm 

# The arguments of NEIN:
                        # theta =  initial values for the finite dimensional parameters   
                        # Z =  observed survival time
                        # d1 = censoring indicator for T
                        # d2 = censoring indicator for C
                        # cop = type of copula
                        # X,W =  data matrices
                        # dist = distribution of dependent censoring, C
                        # eps = convergence error

NEIN = function(theta,Z,d1,d2,X,W,cop,dist,eps){
  
  res = SolveL(theta,Z,d1,d2,X,W,cop,dist)       # estimate the cumulative hazard function
  lsml = res$lambda
  Llrge = res$cumhaz
  T1 = res$times
  longfrm = Longfun(Z,T1,lsml,Llrge)         
  lhat = longfrm[,1]
  cumL = longfrm[,2]
  
  k = dim(X)[2]
  l = dim(W)[2]
  
  #lb = lower boundary, ub = upper boundary
  
  if(cop==1|cop==2)                # Clayton or Frank copula
  {
    lb = c(rep(-Inf,k+l),0,1e-5)   
    ub = c(rep(Inf,k+l+2))         
  }
  if(cop==3)                       
  {
    lb = c(rep(-Inf,k+l),0,1.0001) 
    ub = c(rep(Inf,k+l+2))          
  } 
  if(cop==4)
  {
    lb = c(rep(-Inf,k+l),0,-0.99)  
    ub = c(rep(Inf,k+l+1),0.99)
  } 
  
  # estimate theta
  
  parhat = nlminb(start = theta,PseudoL,Z = Z,d1 = d1,d2 = d2,X = X,W = W,lhat = lhat,cumL = cumL,cop = cop,dist = dist,lower = lb ,upper =  ub, control = list(eval.max=300,iter.max=200))$par
  
  a = theta
  b = parhat
  flag = 0
  
  while (Distance(b,a)>eps){    # doing this while loop until the desired convergence criteria are met, eps
    a = b
    res = SolveL(a,Z,d1,d2,X,W,cop,dist)
    lsml = res$lambda
    Llrge = res$cumhaz
    T1 = res$times
    longfrm = Longfun(Z,T1,lsml,Llrge)
    lhat = longfrm[,1]
    cumL = longfrm[,2]
    parhat = nlminb(start = a,PseudoL,Z = Z,d1 = d1,d2 = d2,X = X,W = W,lhat = lhat,cumL = cumL,cop = cop,dist = dist,lower = lb ,upper =  ub, control = list(eval.max=300,iter.max=200))$par
    
    b = parhat
    flag = flag+1;
    if (flag>50)        # stop after iteration 50; this usually gives sufficient convergence results
    {
      flag=0;
      break;
    }
  }
  cumHat = SolveL(b,Z,d1,d2,X,W,cop,dist)
  output = list(parhat = b, hazard = cumHat$lambda, cumhaz = cumHat$cumhaz,T1 = cumHat$times);
  return(output)
}


# LikCopInd- Pseudo-likelihood function under independent censoring
# The arguments of LikCopInd:
                        # theta = parameter estimates at a given iteration
                        # Z =  observed survival time
                        # d1 = censoring indicator for T
                        # d2 = censoring indicator for C
                        # lhat = estimated hazard function
                        # cumL = estimated cumulative hazard function
                        # X,W =  data matrices
                        # dist = distribution of dependent censoring, C

LikCopInd <- function(theta,Z,d1,d2,X,W,lhat,cumL,dist){ # gamma = 0
  
  if (is.vector(X)){
    k = 1
    beta = theta[k]
    mu1 = X*beta
  }
  else if(!is.vector(X)){
    k = dim(X)[2]
    beta = theta[1:k]
    mu1 = X%*%beta
  }
  l = dim(W)[2]
  eta = theta[(k+1):(l+k)]
  mu2 = W%*%eta
  nu = theta[(l+k+1)]
  y = log(Z)
  
  # distribution of T
  
  g1 = lhat*exp(mu1)*exp(-cumL*exp(mu1))
  G1 = 1-exp(-cumL*exp(mu1))
  
  # dist of C
  
  if (dist == "wb"){                                           # Weibull
    g2 = 1/nu*exp((y-mu2)/nu)*exp(-exp((y-mu2)/nu))*(1/Z)      # density of C
    G2 = 1-exp(-exp((y-mu2)/nu))                               # dist. of C
  }
  if (dist == "lgn"){                                          # Log-normal
    m = (y-mu2)/nu
    g2 = 1/(sqrt(2*pi)*nu*Z)*exp(-(m^2/2))                     # density of C
    G2 = pnorm(m)                                              # dist. of C
  } 
  
  # Joint dist,
  
  G1 = pmax(1e-10,G1)
  G2 = pmax(1e-10,G2)
  
  Logn <- -sum(log(pmax(g1[d1==1], 1e-10)))-sum((log(1-G1[(1-d1)==1])))-sum(log(pmax(g2[d2==1], 1e-10)))-sum(log(1-G2[(1-d2)==1]))
  return(Logn)
}  


# NEINInd - main function in implementing the proposed estimation algorithm
# under independent censoring assumption.
# The arguments are theta, Z, d1, d2, X, W, dist and eps,
# which are defined above under NEIN function

NEINInd = function(theta,Z,d1,d2,X,W,dist,eps){   # Likelihood function independent
  
  res = SolveLI(theta,Z,d1,d2,X)             # independent copula
  lsml = res$lambda
  Llrge = res$cumhaz
  T1 = res$times
  longfrm = Longfun(Z,T1,lsml,Llrge)
  lhat = longfrm[,1]
  cumL = longfrm[,2]
  if (is.vector(X)){
    k = 1
  }
  if(!is.vector(X)){
    k = dim(X)[2]
  }
  l = dim(W)[2]
  lb = c(rep(-Inf,k+l),0)
  ub = c(rep(Inf,k+l+1))
  
  parhat  = nlminb(start = theta,LikCopInd, Z = Z,d1 = d1,d2 = d2,X = X,W = W,lhat = lhat,cumL = cumL, dist = dist, lower = lb ,upper =  ub, control = list(eval.max=300,iter.max=200))$par
  
  a = theta
  b = parhat
  flag = 0
  
  while (Distance(b,a)>eps){
    a = b;
    res = SolveLI(a,Z,d1,d2,X)
    lsml = res$lambda
    Llrge = res$cumhaz
    T1 = res$times
    longfrm = Longfun(Z,T1,lsml,Llrge)
    lhat = longfrm[,1]
    cumL = longfrm[,2]
    
    parhat  = nlminb(start = a,LikCopInd, Z = Z,d1 = d1,d2 = d2,X = X,W = W,lhat = lhat,cumL = cumL, dist = dist, lower = lb ,upper =  ub, control = list(eval.max=300,iter.max=200))$par
    
    b = parhat;
    flag = flag+1;
    
    if (flag>50)
    {
      flag=0;
      break;
    }
  }
  cumHat = SolveLI(b,Z,d1,d2,X)
  output = list(parhat1 = b, hazard = cumHat$lambda, cumhaz = cumHat$cumhaz,T1 = cumHat$times);
  return(output)
} 



#========================================================
# Goodness-of-fit test 
#========================================================


# We provide below functions related to the goodness-of-fit
# test proposed in the Supplementary Material


# DistF2 - distribution function of R = min (T,C) at a single 
# time point

# The arguments of DistF:
                    # par = parameter estimates 
                    # r =   a fixed time point
                    # cumL = estimated cumulative hazard function
                    # X,W =  data matrices
                    # cop = type of copula
                    # dist = distribution of dependent censoring time, C


DistF2 = function(par,r,cumL,X,W,cop,dist){
  k1 = length(X[1,])
  k2 = length(W[1,])
  l = k1+k2
  v = k1+1
  beta = par[1:k1]
  eta = par[v:l]
  nu = par[l+1]
  gm = par[l+2]
  y = log(r)
  
  # distribution of T -- Cox model
  
  G1 = 1-exp(-cumL*exp(X%*%beta))
  
  # dist of C
  
  if (dist == "wb"){    # Weibull
    G2 = 1-exp(-exp((y-W%*%eta)/nu))
  }
  if (dist == "lgn"){   # Lognormal
    m = (y-W%*%eta)/nu
    G2 = pnorm(m)
  }
  if(dist == "gomp"){  # Gompertz 
    mu2 = exp(W%*%par[v:l])
    G2 = 1-exp(exp(-mu2/nu)*(1-exp(r/nu)))
  }
  
  # avoid numerical issues
  
  G1[is.na(G1)] = 1e-10
  G2[is.na(G2)] = 1e-10
  G1[is.nan(G1)] = 1e-10
  G2[is.nan(G2)] = 1e-10
  G1 = pmax(G1,1e-10)
  G1 = pmin(G1,0.99999)
  G2 = pmax(G2,1e-10)
  G2 = pmin(G2,0.999999)
  
  JD = cbind(G1,G2)   # Joint distribution function
  
  # copula contribution
  
  if (cop == 1)     # Clayton copula
  {  
    clay.cop = claytonCopula(gm, dim = 2)
    JDr = pCopula(JD,clay.cop)                              
  }
  if (cop == 2)   # Frank copula
  {
    frank.cop = frankCopula(gm, dim = 2)
    JDr = pCopula(JD,frank.cop)    
  }
  if (cop == 3)    # Gumbel copula
  {  
    gumb.cop = gumbelCopula(gm, dim = 2)
    JDr =  pCopula(JD,gumb.cop)
  }
  if (cop == 4)     # normal copula
  {
    p1 = qnorm(G1)
    p2 = qnorm(G2)
    c2 = cbind(p1,p2)
    JDr = pbivnorm(c2,rho = gm)       
  }
  tot = (G1+G2-JDr)
  return(tot)
} 


#==========================================================================

# CMtest - main function to implement the proposed
# goodness-of-fit test

# The arguments of CMtest:
                     # parhat = parameter estimates 
                     # Z =   observed survival time
                     # Lhat = estimated cumulative hazard function
                     # X,W =  data matrices
                     # cop = type of copula
                     # dist = distribution of dependent censoring time, C
                     # d3 = censoring indicator for R, R = min(T,C)
                     # T1 = distinct observed survival time
                     # B = bootstrap size
                     # eps = convergence error

CMtest = function(parhat,Z,d3,Lhat,T1,X,W,cop,dist,B,eps){
  
  # Empirical dist of R using the Kaplan-Meier estimator
  
  surv = survfit(Surv(Z,d3)~1,type = 'kaplan-meier')
  ED  = 1-get_surv_prob(surv,Z)            
  ED[is.na(ED)] = 0
  
  # Compute dist of R under the fitted model
  
  n = length(Z)
  DP = rep(0,n)
  for (i in 1:n){
    Haz <- Lhat[SearchIndicate(Z[i],T1)]
    contr = DistF2(parhat,Z[i],Haz,X,W,cop,dist)   
    DP[i] = sum(contr)/n
  }
  C0 = sum((ED-DP)^2);       # Cramer-Von-Mises statistic
  
  
  # Bootstrapping the distribution of Cramer-Von-Mises statistic
  
  CVM = rep(0,B)
  iseed = 214399
  for (b in 1:B){
    bseed = iseed+b
    CVM[b] = BootGOF(parhat,Z,d3,Lhat,T1,X,W,cop,dist,eps,bseed)  # Bootstrap function
  }
  count1 <- sum(C0 >= quantile(CVM, prob = 0.9))                  # Compute quantile under bootstrap distribution
  count2 <- sum(C0 >= quantile(CVM, prob = 0.95))
  CVb <- mean(CVM)
  pvalue <- sum(CVM>=C0)/B
  temp <- c(C0,CVb,count1,count2)
  return(temp)
}


# BootGOF = used to determine the distribution of our test statistic
# using bootstrap approximation 

# The arguments of BootGOF are parhat,Z,d3,Lhat,T1,X,W,cop,dist, eps 
# and iseed, which are defined above

BootGOF = function(parhat,Z,d3,Lhat,T1,X,W,cop,dist,eps,iseed){
  set.seed(iseed)
  n = length(Z)
  k1 = length(X[1,])
  k2 = length(W[1,])
  l = k1+k2
  v = k1+1
  beta = parhat[1:k1]         # survival model parameters
  eta = parhat[v:l]           # censoring model parameters
  s2 = parhat[l+1]
  gm = parhat[l+2]
  
  # Copulas
  
  if(cop == 1){    # Clayton copula                                  
    clay.cop <- claytonCopula(param = gm,dim = 2)                 
    dat = rCopula(n,clay.cop)
  }
  if(cop == 2){    # Frank copula
    frank.cop = frankCopula(gm, dim = 2)
    dat = rCopula(n,frank.cop)
  }
  if(cop == 3){   # Gumbel copula
    gumb.cop = gumbelCopula(gm, dim = 2)
    dat = rCopula(n,gumb.cop)
  }
  if(cop == 4){   # normal copula
    norm.cop = normalCopula(gm, dim = 2) 
    dat = rCopula(n,norm.cop)
  }
  
  # Survival time of T1 from fitted Cox model
  
  mu1 = exp(X%*%beta)
  mu2 = exp(W%*%eta)
  
  Lm <- (-log(1-dat[,1]))/mu1
  Tb = rep(0,n)
  
  # Find survival time from the quantile of survival function of T
  
  for (i in 1:n){
    if(Lm[i]<=max(Lhat)){
      ind1 = max(which(Lhat<=Lm[i]))
      Tb[i] = T1[ind1]
    }
    else{
      Tb[i] = 10000       
    }
  }
  
  # Dependent censoring time from a parametric model
  
  if(dist == "gomp"){       # Gompertz model
    Cb = s2*(log(1-log(1-dat[,2])*exp(mu2/s2)))
  }
  if(dist == "wb"){         # Weibull model
    Cb = (-log(1-dat[,2]))^s2*mu2
  }
  if (dist == "lgn"){       # Lognormal model 
    Cb = exp(s2*qnorm(dat[,2]))*mu2
  }
  
  Asurv = survfit(Surv(Z,1-d3)~1,type = 'kaplan-meier')             # distribution of A
  Ab = quantile(Asurv, probs = runif(n,1e-5,0.999999))$quantile     # generate bootstrap A from the dist of A
  Ab = as.numeric(Ab)
  Ab[is.na(Ab)] = max(Z)+1
  
  # bootstrap observations
  
  Zb = pmin(Tb,Cb,Ab)           # observed bootstrapped survival time 
  Zb[Zb==0] = 1e-20
  db1 = as.numeric(Zb==Tb)      # censoring indicator for Tb
  db2 = as.numeric(Zb==Cb)      # censoring indicator for Cb
  db3 = db1+db2
  
  # fit model to bootstrap data
  datb = cbind(Zb,db1,db2,X,W)  
  datB = datb[order(datb[,1]),]
  Zb = datB[,1]
  db1 = datB[,2]
  db2 = datB[,3]
  dm = dim(datB)[2]
  Xb = datB[,4:(k1+3)]
  Wb = datB[,(k1+4):dm]
  
  # model fit
  initb = c(0.45,1,0.9,0.35,0.75,2,5)
  outputb = NEIN(initb,Zb,db1,db2,Xb,Wb,cop,dist,eps)       # Fit proposed model
  parhatb = outputb$parhat                    
  Tb1 = outputb$T1
  Lhatb = outputb$cumhaz
  
  # Empirical dist of Rb subject to right censoring
  
  db3 = db1+db2
  survb = survfit(Surv(Zb,db3)~1,type = 'kaplan-meier')
  EDb  = 1-get_surv_prob(survb, Zb)                         # distribution function
  EDb[is.na(EDb)] = 0
  
  # compute dist of Rb = min(Tb,Cb) from fitted model
  
  DPb = rep(0,n)
  for (i in 1:n){
    Hazb <- Lhatb[SearchIndicate(Zb[i],Tb1)]
    contr = DistF2(parhatb,Zb[i],Hazb,Xb,Wb,cop,dist)   
    DPb[i] = sum(contr)/n
  }
  CVM <- sum((EDb-DPb)^2)    # Cramer-Von-Mises statistic 
  return(CVM)
}



#===================================================================
# Real data analyses related functions
#===================================================================

# BSErealdata(): used to compute the bootstrap standard
# errors 

# arguments:
           # data= real data used for resampling
           # bseed = starting value to generate random number
           # cop = type of copula
           # dist = distribution of censoring model
           # eps = convergence error

BSErealdata = function(data,bseed,cop,dist,eps){ # Bootstrap standard error
  
    set.seed(bseed)
    samp1 = sample(length(data[,1]),replace = TRUE)
    datB = data[samp1,]
    datB = datB[order(datB[,1]),]
    Zb = datB[,1];
    db1 = datB[,2];
    db2 = datB[,3];
    Xb = datB[,4:7];
    Wb = datB[,8:12];
    
    fit1 = coxph(Surv(Zb,db1)~Xb)
    fit2 = survreg(Surv(Zb,db2)~Xb, dist = "lognormal")
    initb = c(fit1$coefficients,fit2$coefficients,fit2$scale,2)
    
    outputb = NEIN(initb,Zb,db1,db2,Xb,Wb,cop,dist,eps)
    tb = outputb$parhat
    
    if(cop == 2){
      taub = tau(frankCopula(tb[11]))
      tb[11] = taub}
    if(cop == 3){
      taub = tau(gumbelCopula(tb[11]))
      tb[11] = taub
    }
    if(cop == 4){
      taub = 2*asin(tb[11])/pi
      tb[11] = taub
    }
    beta= tb
    
  return(beta)
}

GOFrealdata = function(parhat,Z,d3,Lhat,T1,X,W,cop,dist,eps,iseed){
  set.seed(iseed)
  n = length(Z)
  k1 = length(X[1,])
  k2 = length(W[1,])
  l = k1+k2
  v = k1+1
  beta = parhat[1:k1]
  eta = parhat[v:l]
  s2 = parhat[l+1]
  gm = parhat[l+2]
  
  # Copulas
  
  if(cop == 1){  # Clayton                                        # Clayton Copula with weibull margins
    clay.cop <- claytonCopula(param = gm,dim = 2)                 # tau = 0.75
    dat = rCopula(n,clay.cop)
  }
  if(cop == 2){  # Frank 
    frank.cop = frankCopula(gm, dim = 2)
    dat = rCopula(n,frank.cop)
  }
  if(cop == 3){  # Gumbel 
    gumb.cop = gumbelCopula(gm, dim = 2)
    dat = rCopula(n,gumb.cop)
  }
  if(cop == 4){ # normal
    norm.cop = normalCopula(gm, dim = 2) 
    dat = rCopula(n,norm.cop)
  }
  
  # Survival time of T1 from fitted Cox model
  
  mu1 = exp(X%*%beta)
  mu2 = exp(W%*%eta)
  
  Lm <- (-log(1-dat[,1]))/mu1
  Tb = rep(0,n)
  
  for (i in 1:n){
    if(Lm[i]<=max(Lhat)){
      ind1 = max(which(Lhat<=Lm[i]))
      Tb[i] = T1[ind1]
    }
    else{
      Tb[i] = 10000
    }
  }
  
  # Dependent censoring time from a parametric model
  
  if(dist == "gomp"){
    Cb = s2*(log(1-log(1-dat[,2])*exp(mu2/s2)))
  }
  if(dist == "wb"){
    Cb = (-log(1-dat[,2]))^s2*mu2
  }
  if (dist == "lgn"){
    Cb = exp(s2*qnorm(dat[,2]))*mu2
  }
  
  Asurv = survfit(Surv(Z,1-d3)~1,type = 'kaplan-meier')                 # distribution of A
  Ab = quantile(Asurv, probs = runif(n,1e-5,0.999999))$quantile         # quantile of KM estimator
  Ab = as.numeric(Ab)
  Ab[is.na(Ab)] = max(Z)+1
  
  # bootstrap observations
  
  Zb = pmin(Tb,Cb,Ab)
  Zb[Zb==0] = 1e-20
  db1 = as.numeric(Zb==Tb)
  db2 = as.numeric(Zb==Cb)
  db3 = db1+db2
  
  # fit model to bootstrap data
  
  datb = cbind(Zb,db1,db2,X,W)
  datB = datb[order(datb[,1]),]
  Zb = datB[,1]
  db1 = datB[,2]
  db2 = datB[,3]
  dm = dim(datB)[2]
  Xb = datB[,4:(k1+3)]
  Wb = datB[,(k1+4):dm]
  
  # model fit
  
  fit1 = coxph(Surv(Zb,db1)~Xb)
  fit2 = survreg(Surv(Zb,db2)~Wb-1, dist = "lognormal")
  initb = c(fit1$coefficients,fit2$coefficients,fit2$scale,gm)
  
  outputb = NEIN(initb,Zb,db1,db2,Xb,Wb,cop,dist,eps)
  parhatb = outputb$parhat                    
  Tb1 = outputb$T1
  Lhatb = outputb$cumhaz
  
  # Empirical dist of Rb subject to right censoring
  
  db3 = db1+db2
  survb = survfit(Surv(Zb,db3)~1,type = 'kaplan-meier')
  EDb  = 1-get_surv_prob(survb, Zb)                             # distribution function
  EDb[is.na(EDb)] = 0
  
  # compute dist of Rb = min(Tb,Cb) from fitted model
  
  DPb = rep(0,n)
  for (i in 1:n){
    Hazb <- Lhatb[SearchIndicate(Zb[i],Tb1)]
    contr = DistF2(parhatb,Zb[i],Hazb,Xb,Wb,cop,dist)   
    DPb[i] = sum(contr)/n
  }
  CVM <- sum((EDb-DPb)^2)
  return(CVM)
}

BSErealdataI = function(data,bseed,dist,eps){
  
    set.seed(bseed)
    samp1 = sample(length(data[,1]),replace = TRUE)
    datB = data[samp1,]
    datB = datB[order(datB[,1]),]
    Zb = datB[,1];
    db1 = datB[,2];
    db2 = datB[,3];
    Xb = datB[,4:7];
    Wb = datB[,8:12];
    
    fit1 = coxph(Surv(Zb,db1)~Xb)
    fit2 = survreg(Surv(Zb,db2)~Xb,dist = "lognormal")
    initb = c(fit1$coefficients,fit2$coefficients,fit2$scale)
    
    outputb = NEINInd(initb,Zb,db1,db2,Xb,Wb,dist,eps)
    beta = outputb$parhat
    
  return(beta)
}

GOFrealdataI = function(parhat,Z,d3,Lhat,T1,X,W,dist,eps,iseed){
  set.seed(iseed)
  n = length(Z)
  k1 = length(X[1,])
  k2 = length(W[1,])
  l = k1+k2
  v = k1+1
  beta = parhat[1:k1]
  eta = parhat[v:l]
  s2 = parhat[l+1]
  
  # Survival time of T1 from fitted Cox model
  
  mu1 = exp(X%*%beta)
  mu2 = exp(W%*%eta)
  
  # generate two indep unif random vars
  
  norm.cop = normalCopula(0, dim = 2) 
  dat = rCopula(n,norm.cop)
  
  Lm <- (-log(1-dat[,1]))/mu1
  Tb = rep(0,n)
  
  for (i in 1:n){
    if(Lm[i]<=max(Lhat)){
      ind1 = max(which(Lhat<=Lm[i]))
      Tb[i] = T1[ind1]
    }
    else{
      Tb[i] = 10000
    }
  }
  
  # Dependent censoring time from a parametric model
  
  if(dist == "gomp"){
    Cb = s2*(log(1-log(1-dat[,2])*exp(mu2/s2)))
  }
  if(dist == "wb"){
    Cb = (-log(1-dat[,2]))^s2*mu2
  }
  if (dist == "lgn"){
    Cb = exp(s2*qnorm(dat[,2]))*mu2
  }
  
  Asurv = survfit(Surv(Z,1-d3)~1,type = 'kaplan-meier')                 # distribution of A
  Ab = quantile(Asurv, probs = runif(n,1e-5,0.999999))$quantile         # quantile of KM estimator
  Ab = as.numeric(Ab)
  Ab[is.na(Ab)] = max(Z)+1
  
  # bootstrap observations
  
  Zb = pmin(Tb,Cb,Ab)
  Zb[Zb==0] = 1e-20
  db1 = as.numeric(Zb==Tb)
  db2 = as.numeric(Zb==Cb)
  db3 = db1+db2
  
  # fit model to bootstrap data
  
  datb = cbind(Zb,db1,db2,X,W)
  datB = datb[order(datb[,1]),]
  Zb = datB[,1]
  db1 = datB[,2]
  db2 = datB[,3]
  dm = dim(datB)[2]
  Xb = datB[,4:(k1+3)]
  Wb = datB[,(k1+4):dm]
  
  # model fit
  
  fit1 = coxph(Surv(Zb,db1)~Xb)
  fit2 = survreg(Surv(Zb,db2)~Wb-1, dist = "lognormal")
  initb = c(fit1$coefficients,fit2$coefficients,fit2$scale)
  
  outputb = NEINInd(initb,Zb,db1,db2,Xb,Wb,dist,eps)
  parhatb = outputb$parhat                   
  Tb1 = outputb$T1
  Lhatb = outputb$cumhaz
  
  # Empirical dist of Rb subject to right censoring
  
  db3 = db1+db2
  survb = survfit(Surv(Zb,db3)~1,type = 'kaplan-meier')
  EDb  = 1-get_surv_prob(survb, Zb)                             # distribution function
  EDb[is.na(EDb)] = 0
  
  # compute dist of Rb = min(Tb,Cb) from fitted model
  
  DPb = rep(0,n)
  for (i in 1:n){
    Hazb <- Lhatb[SearchIndicate(Zb[i],Tb1)]
    contr = DistFI(parhatb,Zb[i],Hazb,Xb,Wb,dist)   
    DPb[i] = sum(contr)/n
  }
  CVM <- sum((EDb-DPb)^2)
  return(CVM)
}




