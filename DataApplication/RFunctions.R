

# data generating model

DataGenerator = function(n,par,iseed,cop,dist){      # cop- type of copula
  
  set.seed(iseed)
  
  coef1 = par[[1]]
  beta = coef1[1:2]
  ld = coef1[3]
  nu1 = coef1[4]
  
  eta = par[[2]]
  sd = par[[3]]
  nu2 = sd[1]
  gm = sd[2]                  # cop para
  
  x0 = rep(1,n)
  x1 = rbinom(n,1,0.5)
  x2 = rnorm(n,0,1)
  M = cbind(x1,x2)
  W = cbind(x0,x1,x2)
  
  
  if(cop == 1 & dist == "wb"){                                    # Clayton Copula with weibull margins
    clay.cop <- claytonCopula(param = gm,dim = 2)                 # tau = 0.75
    U = rCopula(n,clay.cop)
    u1 = U[,1]
    u2 = U[,2]
    
    # Weibull marginals
    mu1 = exp(M%*%beta)
    mu2 = exp(W%*%eta)
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    T2 = (-log(1-u2))^(nu2)*mu2
    A = runif(n,0,20)
  }
  else if (cop == 2 & dist == "wb"){                                           # frank.cop <- frankCopula(param = rho,dim=2)   # tau = 0.75
    frank.cop <- frankCopula(param = gm,dim = 2)  
    U = rCopula(n,frank.cop)
    u1 = U[,1]
    u2 = U[,2]
    
    # Censoring time C using Weibull 
    mu1 = exp(M%*%beta)
    mu2 = exp(W%*%eta) 
    
    T2 = (-log(1-u2))^(nu2)*mu2
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    A = runif(n,0,15)                             
  }
  else if (cop == 3 & dist == "wb"){                                           # Gumbel with Weibull margins
    gumbel.cop <- gumbelCopula(param = gm,dim = 2)  
    U = rCopula(n,gumbel.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
    
    #  C using Weibull  
    mu1 = exp(M%*%beta)
    mu2 = exp(W%*%eta)
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    T2 = (-log(1-u2))^(nu2)*mu2
    A = runif(n,0,25)                             
  }
  else if (cop == 2 & dist == "lgn"){                     # Gumbel with log-normal margin
    frank.cop <- frankCopula(param = gm,dim = 2)  
    U = rCopula(n,frank.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
    
    # survival time of T and C 
    mu1 = exp(M%*%beta)
    mu2 = W%*%eta  
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    T2 = exp(nu2*qnorm(u2)+mu2)
    A = runif(n,0,25)                                     # 25 = 45% and 40% , 13 = 35% and 40%                      
  }  
  else if (cop == 4 & dist == "lgn"){                     # Gumbel with log-normal margin
    normal.cop <- normalCopula(param = gm,dim = 2)  
    U = rCopula(n,normal.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
    
    # survival time of T and C 
    mu1 = exp(M%*%beta)
    mu2 = W%*%eta  
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    T2 = exp(nu2*qnorm(u2)+mu2)
    A = runif(n,0,25)                                 # 25 = 45% and 40% , 13 = 35% and 40%                      
  }  
  Z = pmin(T1,T2,A)
  d1 = as.numeric(Z==T1)
  d2 = as.numeric(Z==T2)
  data = cbind(Z,d1,d2,M,W)
  return(data)
}

DatasimGof = function(n,par,iseed,cop,dist){      # cop- type of copula
  
  set.seed(iseed)
  
  coef1 = par[[1]]
  beta = coef1[1:2]
  ld = coef1[3]
  nu1 = coef1[4]
  
  eta = par[[2]]
  sd = par[[3]]
  nu2 = sd[1]
  gm = sd[2]                  # copula
  
  x0 = rep(1,n)
  x1 = rbinom(n,1,0.5)
  x2 = rnorm(n,0,1)
  M = cbind(x1,x2)
  W = cbind(x0,x1,x2)
  
  if(cop == 1 & dist == "wb"){                                                  # Clayton Copula with weibull margins
    clay.cop <- claytonCopula(param = gm,dim = 2)                # tau = 0.75
    U = rCopula(n,clay.cop)
    u1 = U[,1]
    u2 = U[,2]
    
    # Weibull marginals
    mu1 = exp(M%*%beta)
    mu2 = exp(W%*%eta)
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    
    T2 = (-log(1-u2))^(nu2)*mu2
    A = runif(n,0,20)
  }
  else if (cop == 2 & dist == "wb"){                                 # frank.cop <- frankCopula(param = rho,dim=2)   # tau = 0.75
    frank.cop <- frankCopula(param = gm,dim = 2)  
    U = rCopula(n,frank.cop)
    u1 = U[,1]
    u2 = U[,2]
    
    # Censoring time C using Weibull 
    mu1 = exp(M%*%beta)
    mu2 = exp(W%*%eta) 
    
    T2 = (-log(1-u2))^(nu2)*mu2
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    A = runif(n,0,15)                             
  }
  else if (cop == 2 & dist == "lgn"){                                # lognormal + Frank copula                                         
    frank.cop <- frankCopula(param = gm,dim = 2)  
    U = rCopula(n,frank.cop)
    u1 = U[,1]
    u2 = U[,2]
    
    # generate T and C 
    mu1 = exp(M%*%beta)
    mu2 = W%*%eta  
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    T2 = exp(nu2*qnorm(u2)+mu2)
    A = runif(n,0,25)  
  }
  else if (cop==3 & dist == "wb"){           # Gumbel with Weibull margins
    gumbel.cop <- gumbelCopula(param = gm,dim = 2)  
    U = rCopula(n,gumbel.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
    
    # C using Weibull  
    mu1 = exp(M%*%beta)
    mu2 = exp(W%*%eta)
    T2 = (-log(1-u2))^(nu2)*mu2
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    A = runif(n,0,25)                             
  }
  else if (cop==4 & dist == "lgn"){                                        # Gumbel with log-normal margin
    normal.cop <- normalCopula(param = gm,dim = 2)  
    U = rCopula(n,normal.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
    
    # survival time of T and C 
    mu1 = exp(M%*%beta)
    mu2 = W%*%eta  
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    T2 = exp(nu2*qnorm(u2)+mu2)
    A = runif(n,0,25)                                 # 25 = 45% and 40% , 13 = 35% and 40%                      
  }  
  else if (cop == 3 & dist == "lgn")
  {
    gumbel.cop <- gumbelCopula(param = gm,dim = 2)  
    U = rCopula(n,gumbel.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
    
    # survival time of T and C 
    mu1 = exp(M%*%beta)
    mu2 = W%*%eta  
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    T2 = exp(nu2*qnorm(u2)+mu2)
    A = runif(n,0,25)     
  }
  
  Z = pmin(T1,T2,A)
  d1 = as.numeric(Z==T1)
  d2 = as.numeric(Z==T2)
  data = cbind(Z,d1,d2,M,W)
  return(data)
}

dat.covmis = function(n,par,iseed,cop,dist){
  set.seed(iseed)
  
  coef1 = par[[1]]
  beta = coef1[1:3]
  ld = coef1[4]
  nu1 = coef1[5]
  
  eta = par[[2]]
  sd = par[[3]]
  nu2 = sd[1]
  gm = sd[2]                  # copula
  
  x0 = rep(1,n)
  x1 = rbinom(n,1,0.5)
  x2 = rnorm(n,0,1)
  x22 = x2^2
  M = cbind(x1,x2,x22)
  W = cbind(x0,x1,x2,x22)
  
  if (cop == 2 & dist == "wb"){                                           # frank.cop <- frankCopula(param = rho,dim=2)   # tau = 0.75
    frank.cop <- frankCopula(param = gm,dim = 2)  
    U = rCopula(n,frank.cop)
    u1 = U[,1]
    u2 = U[,2]
    
    # Censoring time C using Weibull 
    mu1 = exp(M%*%beta)
    mu2 = exp(W%*%eta) 
    
    T2 = (-log(1-u2))^(nu2)*mu2
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    A = runif(n,0,15)                             
  }
  else if (cop == 2 & dist == "lgn"){   # lognormal + Frank copula                                         
    frank.cop <- frankCopula(param = gm,dim = 2)  
    U = rCopula(n,frank.cop)
    u1 = U[,1]
    u2 = U[,2]
    
    # generate T and C 
    mu1 = exp(M%*%beta)
    mu2 = W%*%eta  
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    T2 = exp(nu2*qnorm(u2)+mu2)
    A = runif(n,0,25)  
  }
  else if (cop==3 & dist == "wb"){                            # Gumbel with Weibull margins
    gumbel.cop <- gumbelCopula(param = gm,dim = 2)  
    U = rCopula(n,gumbel.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
    
    # C using Weibull  
    mu1 = exp(M%*%beta)
    mu2 = exp(W%*%eta)
    T2 = (-log(1-u2))^(nu2)*mu2
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    A = runif(n,0,25)                             
  }
  else if (cop==4 & dist == "lgn"){                                 # Gumbel with log-normal margin
    normal.cop <- normalCopula(param = gm,dim = 2)  
    U = rCopula(n,normal.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
    
    # survival time of T and C 
    mu1 = exp(M%*%beta)
    mu2 = W%*%eta  
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    T2 = exp(nu2*qnorm(u2)+mu2)
    A = runif(n,0,25)                                 # 25 = 45% and 40% , 13 = 35% and 40%                      
  }  
  else if (cop == 3 & dist == "lgn")
  {
    gumbel.cop <- gumbelCopula(param = gm,dim = 2)  
    U = rCopula(n,gumbel.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
    
    # survival time of T and C 
    mu1 = exp(M%*%beta)
    mu2 = W%*%eta  
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)
    T2 = exp(nu2*qnorm(u2)+mu2)
    A = runif(n,0,25)     
  }
  
  Z = pmin(T1,T2,A)
  d1 = as.numeric(Z==T1)
  d2 = as.numeric(Z==T2)
  data = cbind(Z,d1,d2,M,W)
  return(data)
}


# Estimating cumulative hazard function

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
    psi = CompC(theta,T1[i],X,W,csum,cop,dist)
    L[i] <- sum(Z == T1[i])/sum((Z>=T1[i])*exp(psi))
  }
  res <- list(lambda = L,cumhaz = cumsum(L), times = T1)
}

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
  G1 = 1-exp(-ld*exp(X%*%beta))
  
  # dist of C
  
  if (dist == "wb"){        # Weibull
    G2 = 1-exp(-exp((y-W%*%eta)/nu))
  }
  if (dist == "lgn"){      # Log-normal
    m = (y-W%*%eta)/nu
    G2 = pnorm(m)
  }
  if (dist == "gomp"){
    Y = Z
    mu2 = exp(W%*%eta)
    G2 = 1-exp(exp(-mu2/nu)*(1-exp(Y/nu)))
  }
  
  G1[is.nan(G1)] = 1e-10
  G2[is.nan(G2)] = 1e-10

  G1[is.na(G1)] <- 1e-10
  G2[is.na(G2)] <- 1e-10
  G1 = pmax(1e-6,G1)
  G2 = pmax(1e-6,G2)
  G1 = pmin(G1,0.999999)
  G2 = pmin(G2,0.999999)
  
  # Joint dist 
  
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

# Likelihood function

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
  # avoid NA
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
  G1 = pmin(G1,0.99999)
  G2 = pmax(G2,1e-10)
  G2 = pmin(G2,0.999999)
  
  PD = cbind(G1,G2)
  
  if (cop == 1)                                         # Clayton copula
  {  
    clay.cop = claytonCopula(gm, dim = 2)
    cp1 = (G1^(-gm)+G2^(-gm)-1)^(-1/gm-1)*G1^(-gm-1)
    cp2 = (G1^(-gm)+G2^(-gm)-1)^(-1/gm-1)*G2^(-gm-1)
    z6 =  1-G1-G2+pCopula(PD,clay.cop)                  # joint survival function       
    z6 =   pmax(z6,1e-10)
  }
  if (cop == 2)                                         # Frank copula
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

NEIN = function(theta,Z,d1,d2,X,W,cop,dist,eps){
  
  res = SolveL(theta,Z,d1,d2,X,W,cop,dist)
  lsml = res$lambda
  Llrge = res$cumhaz
  T1 = res$times
  longfrm = Longfun(Z,T1,lsml,Llrge)
  lhat = longfrm[,1]
  cumL = longfrm[,2]
  
  k = dim(X)[2]
  l = dim(W)[2]
  
  if(cop==1|cop==2)
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
  
  parhat = nlminb(start = theta,PseudoL,Z = Z,d1 = d1,d2 = d2,X = X,W = W,lhat = lhat,cumL = cumL,cop = cop,dist = dist,lower = lb ,upper =  ub, control = list(eval.max=300,iter.max=200))$par
  
  a = theta
  b = parhat
  flag = 0
  
  while (Distance(b,a)>eps){
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
    if(flag==30){
      print(flag)
    }
    
    if (flag>50)
    {
      flag=0;
      break;
    }
  }
  cumHat = SolveL(b,Z,d1,d2,X,W,cop,dist)
  output = list(parhat = b, hazard = cumHat$lambda, cumhaz = cumHat$cumhaz,T1 = cumHat$times);
  return(output)
}


NEIN2 = function(theta,Z,d1,d2,X,W,cop,dist,eps){
  
  res = SolveL(theta,Z,d1,d2,X,W,cop,dist)
  lsml = res$lambda
  Llrge = res$cumhaz
  T1 = res$times
  longfrm = Longfun(Z,T1,lsml,Llrge)
  lhat = longfrm[,1]
  cumL = longfrm[,2]
  
  k = dim(X)[2]
  l = dim(W)[2]
  
  if(cop==1|cop==2)
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
  
  parhat = nlminb(start = theta,PseudoL,Z = Z,d1 = d1,d2 = d2,X = X,W = W,lhat = lhat,cumL = cumL,cop = cop,dist = dist,lower = lb ,upper =  ub, control = list(eval.max=300,iter.max=200))$par
  
  a = theta
  b = parhat
  flag = 0
  
  while (Distance(b,a)>eps){     # catch an error to let the algorithm run
    k = 0
    a = b
    res = SolveL(a,Z,d1,d2,X,W,cop,dist)
    lsml = res$lambda
    Llrge = res$cumhaz
    T1 = res$times
    longfrm = Longfun(Z,T1,lsml,Llrge)
    lhat = longfrm[,1]
    cumL = longfrm[,2]
    
    parhat = try(nlminb(start = a,PseudoL,Z = Z,d1 = d1,d2 = d2,X = X,W = W,lhat = lhat,cumL = cumL,cop = cop,dist = dist,lower = lb ,upper =  ub, control = list(eval.max=300,iter.max=200))$par
                 ,TRUE)
    if(isTRUE(class(parhat)=="try-error"))
    { k = 1
    res = SolveL(theta,Z,d1,d2,X,W,cop,dist)
    lsml = res$lambda
    Llrge = res$cumhaz
    T1 = res$times
    longfrm = Longfun(Z,T1,lsml,Llrge)
    lhat = longfrm[,1]
    cumL = longfrm[,2]
    
    parhat = nloptr(x0 = theta,eval_f = PseudoL,Z = Z,d1 = d1,d2 = d2,X = X,W = W,lhat = lhat,cumL = cumL,cop = cop,dist = dist,lb = lb,ub = ub, eval_g_ineq=NULL,
                    opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_rel"=1.0e-8,"maxeval" = 500 ,"xtol_rel"=rep(1.0e-6)))$solution
    a = theta
    b = parhat
    
    while (Distance(b,a)>eps){
      a = b
      res = SolveL(a,Z,d1,d2,X,W,cop,dist)
      lsml = res$lambda
      Llrge = res$cumhaz
      T1 = res$times
      longfrm = Longfun(Z,T1,lsml,Llrge)
      lhat = longfrm[,1]
      cumL = longfrm[,2]
      
      parhat = nloptr(x0 = a,eval_f = PseudoL,Z = Z,d1 = d1,d2 = d2,X = X,W = W,lhat = lhat,cumL = cumL,cop = cop,dist = dist,lb = lb,ub = ub, eval_g_ineq=NULL,
                      opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_rel"=1.0e-8,"maxeval" = 500 ,"xtol_rel"=rep(1.0e-6)))$solution
      b = parhat
      flag = flag+1;
      if(flag==30){
        print(flag)
      }
      if (flag>50)
      {
        flag=0;
        break;
      }
    }
    }
    b = parhat
    flag = flag+1;
    if(flag==30){
      print(flag)
    }
    if (flag>50)
    {
      flag=0;
      break;
    }
    if (k==1){
      break;
    }
  }
  cumHat = SolveL(b,Z,d1,d2,X,W,cop,dist)
  output = list(parhat = b, hazard = cumHat$lambda, cumhaz = cumHat$cumhaz,T1 = cumHat$times);
  return(output)
}


get_surv_prob <- function(fit, times) {
  stepfun(fit$time, c(1, fit$surv))(times)
}

ff = function(x) {
  dd = x/(exp(x)-1)
  return(dd)
}

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
  
  if (dist == "wb"){                                                   # Weibull
    g2 = 1/nu*exp((y-mu2)/nu)*exp(-exp((y-mu2)/nu))*(1/Z)      # log-time scale
    G2 = 1-exp(-exp((y-mu2)/nu))
  }
  if (dist == "lgn"){                                          # Log-normal
    m = (y-mu2)/nu
    g2 = 1/(sqrt(2*pi)*nu*Z)*exp(-(m^2/2))
    G2 = pnorm(m)
  }
  
  # Joint dist,
  
  G1 = pmax(1e-10,G1)
  G2 = pmax(1e-10,G2)
  
  Logn <- -sum(log(pmax(g1[d1==1], 1e-10)))-sum((log(1-G1[(1-d1)==1])))-sum(log(pmax(g2[d2==1], 1e-10)))-sum(log(1-G2[(1-d2)==1]))
  return(Logn)
}  


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
    if(flag==30){
      print(flag)
    }
    
    #print(b)
    
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



################################
# Goodness-of-fit test functions
################################


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
  
  # distribution of T -- Cox PH

  G1 = 1-exp(-cumL*exp(X%*%beta))
  
  # dist of C
  
  if (dist == "wb"){    # Weibull
    G2 = 1-exp(-exp((y-W%*%eta)/nu))
  }
  if (dist == "lgn"){   # Lognormal
    m = (y-W%*%eta)/nu
    G2 = pnorm(m)
  }
  if(dist == "gomp"){
    mu2 = exp(W%*%par[v:l])
    G2 = 1-exp(exp(-mu2/nu)*(1-exp(r/nu)))
  }

  G1[is.na(G1)] = 1e-10
  G2[is.na(G2)] = 1e-10
  G1[is.nan(G1)] = 1e-10
  G2[is.nan(G2)] = 1e-10
  G1 = pmax(G1,1e-10)
  G1 = pmin(G1,0.99999)
  G2 = pmax(G2,1e-10)
  G2 = pmin(G2,0.999999)
  
  JD = cbind(G1,G2)

  # copula contribution
  
  if (cop == 1)     # Clayton copula
  {  
    clay.cop = claytonCopula(gm, dim = 2)
    JDr = pCopula(JD,clay.cop)                              
  }
  if (cop == 2)  # Frank copula
  {
    frank.cop = frankCopula(gm, dim = 2)
    JDr = pCopula(JD,frank.cop)    
  }
  if (cop == 3)    # Gumbel copula
  {  
    gumb.cop = gumbelCopula(gm, dim = 2)
    JDr =  pCopula(JD,gumb.cop)
  }
  if (cop == 4)     # normal 
  {
    p1 = qnorm(G1)
    p2 = qnorm(G2)
    c2 = cbind(p1,p2)
    JDr = pbivnorm(c2,rho = gm)       
  }
  tot = (G1+G2-JDr)
  return(tot)
} 

CMtest = function(parhat,Z,d3,Lhat,T1,X,W,cop,dist,B,eps){
  
  # Empirical dist of R = min(T,C)...

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
    CVM[b] = BootGOF(parhat,Z,d3,Lhat,T1,X,W,cop,dist,eps,bseed)
  }
  count1 <- sum(C0 >= quantile(CVM, prob = 0.9))        # Compute quantile under bootstrap distribution
  count2 <- sum(C0 >= quantile(CVM, prob = 0.95))
  CVb <- mean(CVM)
  pvalue <- sum(CVM>=C0)/B
  temp <- c(C0,CVb,count1,count2)
  return(temp)
}


BootGOF = function(parhat,Z,d3,Lhat,T1,X,W,cop,dist,eps,iseed){
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
  initb = c(0.45,1,0.9,0.35,0.75,2,5)
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


# Real data analysis functions

BSErealdata = function(data,B,cop,dist,eps){ # Boot sandard error
  
  beta.star = matrix(NA,nrow = B,ncol = 11)
  for(b in 1:B)
  {
    set.seed(2143+b)
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
    beta.star[b,] = tb
  }
  return(beta.star)
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

BSErealdataI = function(data,B,dist,eps){
  
  beta.star = matrix(NA,nrow = B,ncol = 10)
  for(b in 1:B)
  {
    set.seed(2143+b)
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
    beta.star[b,] = outputb$parhat
  }
  return(beta.star)
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



