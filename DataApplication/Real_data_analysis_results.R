

#=====================================================================
# Copula based Cox proportional hazards
# models for dependent censoring
# [Authors and affiliation blinded]
# July, 202
#=====================================================================


#=====================================================================
# Reproducibility: MASTER script real data application
#=====================================================================


# load library
 
library(copula)
library(survival)
library(pbivnorm)


# load also parallel packages

library(parallel)
library(foreach)
library(doParallel)

numCores = detectCores()
numCores = numCores-1
clust = makeCluster(numCores)
registerDoParallel(clust)


#===========================================================================
# load estimation algorithm
#===========================================================================


# The functions used to implement the proposed model on the real
# data example are given in estimation_algorithm_functions.R

source("estimation_algorithm_functions.R")     # load estimation algorithm


#============================================================================
# Follicular cell lymphoma data set
# This data set is analyzed in Section 7 of the paper
#============================================================================

# The data are available in the R package 'randomForestSRC' publicly for
# free


library(randomForestSRC)
data(follic)                          # load the data 




#================================================================================
# Change the data to the form suitable for analysis
#================================================================================

# obtain variables

follic = follic[order(follic$time),]
Z = round(follic$time,digits = 3)
d1 = as.numeric(follic$status==1)
d2 = as.numeric(follic$status==2)
x1 = as.numeric(follic$ch=="Y")                     # treatment
x2 = (follic$age-mean(follic$age))/sd(follic$age)   # standardized age 
x3 = (follic$hgb-mean(follic$hgb))/sd(follic$hgb)   # standardized hemoglobin 
x4 = as.numeric(follic$clinstg==1)                  # clinical stage 
X = cbind(x1,x2,x3,x4)                    # data matrix for survival time T
W = cbind(rep(1,length(Z)),X)             # data matrix for dependent censoring C
data = cbind(Z,d1,d2,X,W)
d3 = d1+d2                                # censoring indictor for min(T,C)
B = 200             # Bootstrap size to compute standard errors





#=================================================================================
# Reproduce Table 4
# Case 1: Cox proportional hazards model and Weibull model 
# coupled by Gumbel copula
#=================================================================================


# Get initial values by fitting a Cox model and Weibull model
# under the assumption of the independent censoring



fit1 = coxph(Surv(Z,d1)~X)                         
fit2 = survreg(Surv(Z,d2)~X, dist = "lognormal")
cop = 3         # Gumbel copula
eps = 1e-4
init = c(fit1$coefficients,fit2$coefficients,fit2$scale,2) # initial values
dist = "wb"                   # Weibull distribution 
outputW = NEIN(init,Z,d1,d2,X,W,cop,dist,eps)         # Fit the proposed model
parhat = outputW$parhat                  # obtain parameter estimates


# compute bootstrap standard error

beta.star <- foreach(b = 1:B, .combine = rbind) %dopar% {   # parallel computing
      source("estimation_algorithm_functions.R")
      bseed = 2143+b
      bootsd = BSErealdata(data,bseed,cop,dist,eps)
      bootsd
}
Bootse = apply(beta.star,2,sd)



# Summary results 

taub = tau(gumbelCopula(parhat[11]))
parhat[11] = taub
zratio = parhat/Bootse
pvalue = 2*( 1-pnorm(abs(zratio)))       # p-value

out.res = cbind(parhat,Bootse,pvalue)    # combine results

colnames(out.res) <- c("Est.", "BSE", "Pvalue")
rownames(out.res) <- c("Treatment","Age","Hgb","Clinstg","Intercept","Treatment"
                       ,"Age","Hgb","Clinstg","sigma","tau")

print("Table 4: Cox-Weibull model")

print(out.res)    # print the output




#=====================================================================
# Goodness-of-fit test for Cox-Weibull model
# see Table 8 (in the supplementary material)
#=====================================================================


# The test statistic is based on the distribution
# of R = min(T,C)

parhat = outputW$parhat
Lhat = outputW$cumhaz
T1 = outputW$T1

# Compute nonparametric part

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
C0 = sum((ED-DP)^2);          # Cramer-Von-Mises

# 0.05032167

# Bootstrapping the distribution of C0


CVM <- foreach(b = 1:B, .combine = rbind) %dopar% {   # parallel computing 
  bseed = 214399+b
  source("estimation_algorithm_functions.R")
  Tb = GOFrealdata(parhat,Z,d3,Lhat,T1,X,W,cop,dist,eps,bseed)
  Tb
}
CVb <- mean(CVM)
pvalueW <- sum(CVM>=C0)/B              # goodness-of-fit pvalue
pvalueW

### 0.365



#=================================================================================
# Reproduce Table 4
# Case 2: Cox proportional hazards model and Lognormal model 
# coupled by Gumbel copula
#=================================================================================

# obtain initial values

fit1 = coxph(Surv(Z,d1)~X)
fit2 = survreg(Surv(Z,d2)~X, dist = "lognormal")
cop = 3   # Gumbel copula
init = c(fit1$coefficients,fit2$coefficients,fit2$scale,2)
dist = "lgn"
outputL = NEIN(init,Z,d1,d2,X,W,cop,dist,eps)  # Fit the proposed model
parhat = outputL$parhat
parhat


# Bootstrap data to obtain standard error

beta.star = foreach(b = 1:B, .combine = rbind) %dopar% {   # parallel computing
  source("estimation_algorithm_functions.R")
  bseed = 2143+b
  bootsd = BSErealdata(data,bseed,cop,dist,eps)
  bootsd
}

Bootse = apply(beta.star,2,sd)           # bootstrap standard error

# Summary results 

taub = tau(gumbelCopula(parhat[11]))
parhat[11] = taub
zratio = parhat/Bootse
pvalue = 2*( 1-pnorm(abs(zratio)))       # p-value

out.res = cbind(parhat,Bootse,pvalue)    # combine results

colnames(out.res) <- c("Est.", "BSE", "Pvalue")
rownames(out.res) <- c("Treatment","Age","Hgb","Clinstg","Intercept","Treatment"
                       ,"Age","Hgb","Clinstg","sigma","tau")

print("Table 4: Cox-Lognormal model")
print(out.res)                           # print the output




#==============================================================================
# Goodness-of-fit test for the Cox-Lognormal model
# see Table 8 (in the supplementary material)
#==============================================================================

parhat = outputL$parhat
Lhat = outputL$cumhaz
T1 = outputL$T1

# Compute nonparametric part

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
C0 = sum((ED-DP)^2);                                  # Cramer-Von-Mises
#0.07584506

# Bootstrapping the distribution of C0

CVM <- foreach(b = 1:B, .combine = rbind) %dopar% {    # parallel computing 
  bseed = 214399+b
  source("RFunctions.R")
  Tb = GOFrealdata(parhat,Z,d3,Lhat,T1,X,W,cop,dist,eps,bseed)
  Tb
}
CVb <- mean(CVM)
pvalueL <- sum(CVM>=C0)/B                               # goodness-of-fit pvalue
pvalueL

### 0.56 



#=================================================================================
# Reproduce Table 5
# Case 1: Cox proportional hazards model and Weibull model 
# under the independence copula
#=================================================================================


# obtain initial values

fit1 = coxph(Surv(Z,d1)~X)
fit2 = survreg(Surv(Z,d2)~X, dist = "lognormal")
init = c(fit1$coefficients,fit2$coefficients,fit2$scale)
dist = "wb"   
outputI = NEINInd(init,Z,d1,d2,X,W,dist,eps)  # Fit the independence copula
parhat = outputI$parhat

# compute bootstrap standard error

beta.star = foreach(b = 1:B, .combine = rbind) %dopar% {    # parallel computing
  source("estimation_algorithm_functions.R")
  bseed = 2143+b
  bootsd = BSErealdataI(data,bseed,cop,dist,eps)
  bootsd
}
Bootse = apply(beta.star,2,sd)



# Summary results 

zratio = parhat/Bootse
pvalue = 2*( 1-pnorm(abs(zratio)))        # pvalue

out.res = cbind(parhat,Bootse,pvalue)     # combine results

colnames(out.res) <- c("Est.", "BSE", "Pvalue")
rownames(out.res) <- c("Treatment","Age","Hgb","Clinstg","Intercept","Treatment"
                       ,"Age","Hgb","Clinstg","sigma")
print("Table 5: Cox-Weibul model")
print(out.res)                             # print the output




#==============================================================================
# Goodness-of-fit test for the Cox-Weibull model
# see Table 8 (in the supplementary material)
#==============================================================================


parhat = c(outputI$parhat)
Lhat = outputI$cumhaz
T1 = outputI$T1

# Compute nonparametric part

surv = survfit(Surv(Z,d3)~1,type = 'kaplan-meier')
ED  = 1-get_surv_prob(surv,Z)            
ED[is.na(ED)] = 0

# Compute dist of R under the fitted model

n = length(Z)
DP = rep(0,n)
for (i in 1:n){
  Haz <- Lhat[SearchIndicate(Z[i],T1)]
  contr = DistFI(parhat,Z[i],Haz,X,W,dist)   
  DP[i] = sum(contr)/n
}
C0 = sum((ED-DP)^2);                      # Cramer-Von-Mises
# 0.07198023

# Bootstrapping the distribution of C0

CVM = rep(0,B)
for (b in 1:B){
  bseed = 214399+b
  CVM[b] = GOFrealdataI(parhat,Z,d3,Lhat,T1,X,W,dist,eps,bseed)
}
CVb <- mean(CVM)
pvalueWI <- sum(CVM>=C0)/B                 # goodness-of-fit p-value 
pvalueWI

## 0.08


#=================================================================================
# Reproduce Table 5
# Case 2: Cox proportional hazards model and Lognormal model 
# under the independence copula
#=================================================================================

# obtain initial values

fit1 = coxph(Surv(Z,d1)~X)
fit2 = survreg(Surv(Z,d2)~X, dist = "lognormal")
init = c(fit1$coefficients,fit2$coefficients,fit2$scale)
dist = "lgn"   
outputI = NEINInd(init,Z,d1,d2,X,W,dist,eps)
parhat = output$parhat


# compute bootstrap standard error

beta.star = foreach(b = 1:B, .combine = rbind) %dopar% {   # parallel computing
  source("estimation_algorithm_functions.R")
  bseed = 2143+b
  bootsd = BSErealdataI(data,bseed,cop,dist,eps)
  bootsd
}

Bootse = apply(beta.star,2,sd)


# Summary results 

zratio = parhat/Bootse
pvalue = 2*( 1-pnorm(abs(zratio)))        # pvalue
out.res = cbind(parhat,Bootse,pvalue)     # combine results

colnames(out.res) <- c("Est.", "BSE", "Pvalue")
rownames(out.res) <- c("Treatment","Age","Hgb","Clinstg","Intercept","Treatment"
                       ,"Age","Hgb","Clinstg","sigma")
print("Table 5: Cox-Lognormal model")
print(out.res) 





#=====================================================================
# Goodness-of-fit test for the Cox-Lognormal model under
# the independence copula, see Table 8 (in the supplementary material)
#=====================================================================



parhat = outputI$parhat
Lhat = outputI$cumhaz
T1 = outputI$T1

# Compute nonparametric part

surv = survfit(Surv(Z,d3)~1,type = 'kaplan-meier')
ED  = 1-get_surv_prob(surv,Z)            
ED[is.na(ED)] = 0

# Compute dist of R under the fitted model

n = length(Z)
DP = rep(0,n)
for (i in 1:n){
  Haz <- Lhat[SearchIndicate(Z[i],T1)]
  contr = DistFI(parhat,Z[i],Haz,X,W,dist)   
  DP[i] = sum(contr)/n
}
C0 = sum((ED-DP)^2);                      # Cramer-Von-Mises


# Bootstrapping the distribution of C0

CVM = rep(0,B)
for (b in 1:B){
  bseed = 214399+b
  CVM[b] = GOFrealdataI(parhat,Z,d3,Lhat,T1,X,W,dist,eps,bseed)
}

CVb <- mean(CVM)
pvalueLI <- sum(CVM>=C0)/B                  # goodness-of-fit p-value 
pvalueLI
# 0




#==============================================================================
# Cumulative hazards plot, produce Figure 5
#==============================================================================

# Cox-Weibull

l1 = outputW$hazard
L1 = outputW$cumhaz
T1 = outputW$T1
time = seq(0,31,0.01)
longfrm = Longfun(time,T1,l1,L1)
cumW = longfrm[,2]


plot(time,cumW, type = "l", cex = 1.5, ylim = c(0,2.3), xlab = "Time (in years)",
          ylab = "Cumulative hazard")


# Cox lognormal

l1 = outputL$hazard
L1 = outputL$cumhaz
T1 = outputL$T1
longfrm = Longfun(time,T1,l1,L1)
cumL = longfrm[,2]

lines(Z,cumL, lty = 2, type = "l")
lines(time,cumL, lty = 2, cex = 1.5)

# Independent copula

l1 = outputI$hazard
L1 = outputI$cumhaz
T1 = outputI$T1
longfrm = Longfun(time,T1,l1,L1)
cumLI = longfrm[,2]

lines(time,cumLI, lty = 3, cex = 1.5,type = "l")


stopCluster(clust)




