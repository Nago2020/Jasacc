
#==========================================================
# Copula based Cox proportional hazards
# models for dependent censoring
# [Authors and affiliation blinded]
# July, 2022
#==========================================================


#==========================================================
# Data generating functions
#==========================================================

# This file contains functions for simulating data used in 
# Section 6 of the paper


# To simulate data used in Scenario 1--4 and 
# Scenario 5, Case 1 and 2 (see Section 6)


# load libraries 

library(copula)
library(survival)
library(pbivnorm)


DataGenerator = function(n,par,iseed,cop,dist){ # cop= type of copula
                                                # n =  sample size
  set.seed(iseed)                               # par = true parameter values
                                                # dist = distribution of dependent censoring, C
  coef1 = par[[1]]                              # iseed = starting value to generate random number 
  beta = coef1[1:2]             # Survival model true parameter value
  ld = coef1[3]
  nu1 = coef1[4]
  
  eta = par[[2]]                # dependent censoring model true parameter value
  sd = par[[3]]
  nu2 = sd[1]
  gm = sd[2]                    # copula parameter
  
  x0 = rep(1,n)
  x1 = rbinom(n,1,0.5)         # x1  = generated from Bernoulli distribution
  x2 = rnorm(n,0,1)            # x2 =  generated from standard normal distribution
  M = cbind(x1,x2)
  W = cbind(x0,x1,x2)

  if(cop == 1 & dist == "wb"){                           
    clay.cop <- claytonCopula(param = gm,dim = 2)    # Clayton Copula               
    U = rCopula(n,clay.cop)
    u1 = U[,1]
    u2 = U[,2]
    
    mu1 = exp(M%*%beta)
    mu2 = exp(W%*%eta)
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)   # Cox model for T
    T2 = (-log(1-u2))^(nu2)*mu2          # Weibull marginal for C
    A = runif(n,0,20)     # A = generated from uniform distribution
  }
  else if (cop == 2 & dist == "wb"){                
    frank.cop <- frankCopula(param = gm,dim = 2)   # Frank copula
    U = rCopula(n,frank.cop)
    u1 = U[,1]
    u2 = U[,2]
    
    mu1 = exp(M%*%beta)
    mu2 = exp(W%*%eta) 
    T2 = (-log(1-u2))^(nu2)*mu2         # Weibull model for C
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)  # Cox model for T
    A = runif(n,0,15)       # A = generated from uniform distribution                       
  }
  else if (cop == 3 & dist == "wb"){                                    
    gumbel.cop <- gumbelCopula(param = gm,dim = 2)    # Gumbel copula
    U = rCopula(n,gumbel.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
    
    mu1 = exp(M%*%beta)
    mu2 = exp(W%*%eta)
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)   # Cox model for T
    T2 = (-log(1-u2))^(nu2)*mu2          # Weibull model for C
    A = runif(n,0,25)       # A = generated from uniform distribution        
  }
  else if (cop == 2 & dist == "lgn"){                     
    frank.cop <- frankCopula(param = gm,dim = 2)   # Frank copula
    U = rCopula(n,frank.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
  
    mu1 = exp(M%*%beta)
    mu2 = W%*%eta  
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)  # Cox model for T
    T2 = exp(nu2*qnorm(u2)+mu2)         # Lognormal margin for C
    A = runif(n,0,25)     # A = generated from uniform distribution                                                    
  }  
  else if (cop == 4 & dist == "lgn"){                     
    normal.cop <- normalCopula(param = gm,dim = 2)  # Gaussian copula
    U = rCopula(n,normal.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
    
    # survival time of T and C 
    mu1 = exp(M%*%beta)
    mu2 = W%*%eta  
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)    # Cox model for T
    T2 = exp(nu2*qnorm(u2)+mu2)           # Lognormal margin for C
    A = runif(n,0,25)                     # A = generated from uniform distribution                                       
  }  
  Z = pmin(T1,T2,A)             # Observed survival time
  d1 = as.numeric(Z==T1)        # censoring indicator for C
  d2 = as.numeric(Z==T2)        # censoring indicator for A
  data = cbind(Z,d1,d2,M,W)
  return(data)
}


# To simulate data used in Scenario 5, Case 3 (see the Supplemetary Material)
# The arguments of dat.covmis:
                              # cop = type of copula
                              # par = true parameter values  
                              # n = sample size
                              # iseed = starting value to generate random number 
                              # dist = distribution of dependent censoring, C

dat.covmis = function(n,par,iseed,cop,dist){
  set.seed(iseed)
  
  coef1 = par[[1]]
  beta = coef1[1:3]             # Survival model true parameter value
  ld = coef1[4]
  nu1 = coef1[5]
  
  eta = par[[2]]               # dependent censoring model true parameter value
  sd = par[[3]]
  nu2 = sd[1]
  gm = sd[2]                   # copula parameter
  
  x0 = rep(1,n)
  x1 = rbinom(n,1,0.5)         # x1  = generated from Bernoulli distribution
  x2 = rnorm(n,0,1)            # x2  = generated from the standard normal distribution
  x22 = x2^2
  M = cbind(x1,x2,x22)
  W = cbind(x0,x1,x2,x22)
  
  if (cop == 2 & dist == "wb"){                                           
    frank.cop <- frankCopula(param = gm,dim = 2)    # Frank copula
    U = rCopula(n,frank.cop)
    u1 = U[,1]
    u2 = U[,2]
    
  
    mu1 = exp(M%*%beta)
    mu2 = exp(W%*%eta) 
    
    T2 = (-log(1-u2))^(nu2)*mu2                   # Weibull model for C
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)            # Cox model for T
    A = runif(n,0,15)                             
  }
  else if (cop == 2 & dist == "lgn"){                                         
    frank.cop <- frankCopula(param = gm,dim = 2)    # Frank copula   
    U = rCopula(n,frank.cop)
    u1 = U[,1]
    u2 = U[,2]
    
    # generate T and C 
    mu1 = exp(M%*%beta)
    mu2 = W%*%eta  
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)              # Cox model for T
    T2 = exp(nu2*qnorm(u2)+mu2)                     # Lognormal model for C
    A = runif(n,0,25)  
  }
  else if (cop==3 & dist == "wb"){                            
    gumbel.cop <- gumbelCopula(param = gm,dim = 2)   # Gumbel copula
    U = rCopula(n,gumbel.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
    
     
    mu1 = exp(M%*%beta)
    mu2 = exp(W%*%eta)
    T2 = (-log(1-u2))^(nu2)*mu2                  # Weibull model for C
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)           # Cox model for T
    A = runif(n,0,25)                             
  }
  else if (cop == 4 & dist == "lgn"){                                 
    normal.cop <- normalCopula(param = gm,dim = 2)   # Gaussian copula
    U = rCopula(n,normal.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
    
    
    mu1 = exp(M%*%beta)
    mu2 = W%*%eta  
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)   # Cox Model for T
    T2 = exp(nu2*qnorm(u2)+mu2)          # Lognormal model for C
    A = runif(n,0,25)                    # A = simulated from uniform distribution                                                
  }  
  else if (cop == 3 & dist == "lgn")
  {
    gumbel.cop <- gumbelCopula(param = gm,dim = 2)  # Gumbel copula
    U = rCopula(n,gumbel.cop)
    U = pmin(U,0.999)
    u1 = U[,1]
    u2 = U[,2]
    
    
    mu1 = exp(M%*%beta)
    mu2 = W%*%eta  
    T1 = (-log(1-u1)/(ld*mu1))^(1/nu1)     # Cox Model for T
    T2 = exp(nu2*qnorm(u2)+mu2)            # Lognormal model for C
    A = runif(n,0,25)                      # A = simulated from uniform distribution 
  }
  
  Z = pmin(T1,T2,A)             # Observed survival time
  d1 = as.numeric(Z==T1)        # censoring indicator for C
  d2 = as.numeric(Z==T2)        # censoring indicator for A
  data = cbind(Z,d1,d2,M,W)
  return(data)
}