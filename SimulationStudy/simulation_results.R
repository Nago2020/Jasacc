
#==========================================================
# Copula based Cox proportional hazards
# models for dependent censoring
# [Authors and affiliation blinded]
# July, 2022
#==========================================================




# Function to rbind in parallel

comb<- function(...){
  mapply('rbind', ..., SIMPLIFY = F)
}


#=============================================================================
# To reproduce simulation results in Section 6
# Tables 1-3 (in the paper) and Tables 6-7 (in the supplementary Material)
# Figures 1-2 (in the paper) and Figures 3-4 (in the supplementary Material)
#=============================================================================




#=============================================================================
# Scenario 1: To reproduce Table 1 and Figure 1
#=============================================================================

# In this scenario we assess the performance of the
# proposed model with the model that assumes
# that T and C are independent


# The arguments for Scenario1:

                        # n = sample size
                        # nsim = number of simulations
                        # iseed = starting value to generate random number
                        # eps = convergence error



Scenario1 = function(n,nsim,iseed,eps){
  
  # Specify the setting here
  
  gamma = c(1.861,4.2,18.2)               # we consider three copula parameters
  cop = 2                                 # 2 = Frank, 3 = Gumbel, 4 = normal
  dist = "wb"                             # wb = Weibull distribution for dependent Censoring C  
                  
 summary = c()
 summaryI = c()
 
  for (j in gamma)
  {
    results <- foreach(i = 1:nsim, .combine = 'comb', .multicombine = T, .maxcombine = 500, .inorder = F)%dopar%{ 
    
      source("data_simulation_functions.R")                   # load the functions with required libraries 
      source("estimation_algorithm_functions.R")
      
  
      par = list(beta = c(0.45,1,0.25,0.75),eta = c(1.35,0.3,1), sd = c(1,j))     # true regression parameters 
      
      data = DataGenerator(n, par, iseed+i,cop,dist)                                # 2 =  generate data from Frank copula  
      data = data[order(data[,1]),] 
      Z = data[,1];
      Z[Z==0] = 1e-10
      d1 = data[,2];
      d2 = data[,3];
      X = data[,4:5];
      W = data[,6:8];
      
      # Fit a Frank  copula model
        
        init = c(0.3,0.75,1.1,0.25,0.8,0.9,j*0.75)     # initial values
        output = NEIN(init,Z,d1,d2,X,W,cop,dist,eps)   
        parhat = output$parhat
        tau = tau(frankCopula(parhat[7]))
        parhat[7] = tau
        
        
        # Cumulative hazard functions
        
        # Frank copula
        T1 = output$T1
        Lhat = output$cumhaz
        time = seq(0.1,10,by = 0.1)
        LD1 = rep(0,length(time))
        for (l in 1:length(time)){                        # obtain estimated cumulative hazards
          LD1[l] = Lhat[SearchIndicate(time[l],T1)]
        }
        
        # Fit Independent copula
        
        initI = c(0.3,0.75,1.1,0.25,0.8,0.9)               # initial values
        outputI = NEINInd(initI,Z,d1,d2,X,W,dist,eps)
        parhatI = outputI$parhat1
          
        # Cumulative hazard functions
          
        # Independence copula
          T2 = outputI$T1
          Lhat1 = outputI$cumhaz
          time = seq(0.1,10,by = 0.1)
          LI = rep(0,length(time))
          for (l in 1:length(time)){
            LI[l] = Lhat1[SearchIndicate(time[l],T2)]
          }
          list(parhat,parhatI,LD1,LI)
      }
    
    # Summary measures
    
    parhat = results[[1]]
    parhatI = results[[2]]
    if(j == 18.2){
      Lhat = results[[3]]
      LhatI = results[[4]]
    }
    
    par1 = c(0.45,1,1.35,0.3,1,1)    # True model parameters
    
    # Frank copula model
    
    par0 = c(par1,tau(frankCopula(j)))
    
    par0m = matrix(par0,nsim,7,byrow = TRUE)
    Bias = apply(parhat-par0m,2,mean)              # Bias
    ESD = apply(parhat,2,sd)                       # Empirical standard deviation
    RMSE = sqrt(apply((parhat-par0m)^2,2,mean))    # Root mean square error
    sum = cbind(Bias,ESD,RMSE)
    summary <- cbind(summary,sum)                  # Merge the tables, for tau = 0.2, 0.4 and 0.8
  

    # Independent copula model
    
    par0I = matrix(par1,nsim,6,byrow=TRUE)
    Bias = apply(parhatI-par0I,2,mean)              # Bias
    ESD = apply(parhatI,2,sd)                       # Empirical standard deviation
    RMSE = sqrt(apply((parhatI-par0I)^2,2,mean))    # Root mean square error
    sumI = cbind(Bias,ESD,RMSE)
    summaryI <- cbind(summaryI,sumI)                # Merge the tables, for tau = 0.2, 0.4 and 0.8
  }
 
 
#-------------------
# produce Table 1
#-------------------
 
# Frank copula
 
 summary = round(summary,3)
 summaryI = round(summaryI,3)
 
 colnames(summary) <- c("Bias", "ESD", "RMSE", "Bias", "ESD", "RMSE", "Bias", "ESD", "RMSE")
 rownames(summary) <- c("beta1", "beta2", "eta0", "eta1", "eta2", "sigma","tau")
 
 print("This is Table 1")
 
 print("For tau = 0.2, 0.4, 0.8 from left to right")
 print("Frank Copula")
 print(summary)
 
# Independent copula

 colnames(summaryI) <- c("Bias", "ESD", "RMSE", "Bias", "ESD", "RMSE", "Bias", "ESD", "RMSE")
 rownames(summaryI) <- c("beta1", "beta2", "eta0", "eta1", "eta2", "sigma")
 
 print("Independence Copula")
 print(summaryI)
 
 
 # Make a cumulative hazard plot
 
 # Frank
 
 time = seq(0.1,10,by = 0.1)
 meanF = apply(Lhat,2,mean)           # average of the estimated cumulative hazard functions under the Frank copula
 
 # independent
 
 meanI = apply(LhatI,2,mean)          # average of the estimated cumulative hazard functions under the independence copula
 

 # True curve
 
 Lbd = (0.25*time^(3/4))
 
#--------------------
# produce Figure 1
#--------------------
 
 plot(time,meanF, type = "l",ylim = c(0,1.6),lty = 2,lwd = 2,col = "grey", xlab = "Time", ylab = "Cumulative hazard", cex=3, cex.lab=1.5,cex.axis=1.5)
 title(main = expression(paste(tau, " = Scenario 1, Figure 1")), line = 0.75)
 lines(time,Lbd,lty = 1,lwd = 2,col = "black", cex = 3)
 lines(time,meanI,lty = 2,lwd = 2,col = "black",cex = 2)
}




#===========================================================================
# Scenario 2: To reproduce Table 2 and Figure 3 (left panel)
#===========================================================================

# Under this scenario we would like to know whether the 
# estimated model parameters are sensitive to the 
# misspecification of the copula structure


# The arguments for Scenario2: we use the same arguments as 
# Scenario1, and estimate the proposed model based on
# the Gumbel and Gaussian copulas 
# when the data are generated from a Frank copula 



Scenario2 = function(n,nsim,iseed,eps){
  
  # Specify the setting here
  
  gamma = c(1.861,4.2,18.2)               # we consider three copula parameters
  cop = 2                                 # 2 = Frank, 3 = Gumbel, 4 = Gaussian 
  dist = "wb"                             # wb = Weibull distribution for dependent Censoring C  
  summary = c()
  summaryI = c()
  
  for (j in gamma)
  {
    results <- foreach(i = 1:nsim, .combine = 'comb', .multicombine = T, .maxcombine = 500, .inorder = F)%dopar%{ 
      
      source("data_simulation_functions.R")               # load the functions with required libraries 
      source("estimation_algorithm_functions.R")
      
      
      par = list(beta = c(0.45,1,0.25,0.75),eta = c(1.35,0.3,1), sd = c(1,j))     # true regression parameters 
      
      data = DataGenerator(n, par, iseed+i,cop,dist)                                # 2 =  generate data from Frank copula  
      data = data[order(data[,1]),] 
      Z = data[,1];
      Z[Z==0] = 1e-10
      d1 = data[,2];
      d2 = data[,3];
      X = data[,4:5];
      W = data[,6:8];
      
      #  Fit a misspecified Gumbel copula model
      
      init = c(0.3,0.75,1.1,0.25,0.8,0.9,1.5)        # initial values
      outputg = NEIN(init,Z,d1,d2,X,W,3,dist,eps)    # 3 =  Gumbel copula
      parhatg = outputg$parhat
      taug  = tau(gumbelCopula(parhatg[7]))
      parhatg[7] = taug
      
      
      # Cumulative hazard functions
      
      # Gumbel copula
      T1 = outputg$T1
      Lhatg = outputg$cumhaz
      time = seq(0.1,10,by = 0.1)
      LDg = rep(0,length(time))
      for (l in 1:length(time)){                        # obtain estimated cumulative hazards
        LDg[l] = Lhatg[SearchIndicate(time[l],T1)]
      }
      
      #  Fit a misspecified Gaussian copula model
      
      init = c(0.3,0.75,1.1,0.25,0.8,0.9,0.3)          # initial values
      outputn = NEIN(init,Z,d1,d2,X,W,4,dist,eps)      # 4 = Gaussian copula
      parhatn = outputn$parhat
      taun = tau(normalCopula(parhatn[7]))
      parhatn[7] = taun
      
      
      # Cumulative hazard functions
      
      # Gumbel copula
      T1 = outputn$T1
      Lhatn = outputn$cumhaz
      time = seq(0.1,10,by = 0.1)
      LDn = rep(0,length(time))
      for (l in 1:length(time)){                        # obtain estimated cumulative hazards
        LDn[l] = Lhatn[SearchIndicate(time[l],T1)]
      }
      list(parhatg,parhatn,LDg,LDn)
    }
    
    # Summary measures
    
    parhat = results[[1]]        # parameter estimates for Gumbel copula
    parhatn = results[[2]]       # parameter estimates for Gaussian copula
    if(j == 18.2){               # cumulative hazards estimates
      Lhat = results[[3]]        
      Lhatn = results[[4]]
    }
    
    par1 = c(0.45,1,1.35,0.3,1,1)    # True model parameters
    
    # Gumbel copula model
    
    par0 = c(par1,tau(frankCopula(j)))
    
    par0m = matrix(par0,nsim,7,byrow = TRUE)
    Bias = apply(parhat-par0m,2,mean)              # Bias
    ESD = apply(parhat,2,sd)                       # Empirical standard deviation
    RMSE = sqrt(apply((parhat-par0m)^2,2,mean))    # Root mean square error
    sum = cbind(Bias,ESD,RMSE)
    summary <- cbind(summary,sum)                  # Merge the tables, for tau = 0.2, 0.4 and 0.8
    
    
    # Gaussian copula model
  
    Bias = apply(parhatn-par0m,2,mean)              # Bias
    ESD = apply(parhatn,2,sd)                       # Empirical standard deviation
    RMSE = sqrt(apply((parhatn-par0m)^2,2,mean))    # Root mean square error
    sumI = cbind(Bias,ESD,RMSE)
    summaryI <- cbind(summaryI,sumI)                # Merge the tables, for tau = 0.2, 0.4 and 0.8
  }
  
  
  #-------------------
  # produce Table 2
  #-------------------
  
  # Gumbel copula
  
  summary = round(summary,3)
  colnames(summary) <- c("Bias", "ESD", "RMSE", "Bias", "ESD", "RMSE", "Bias", "ESD", "RMSE")
  rownames(summary) <- c("beta1", "beta2", "eta0", "eta1", "eta2", "sigma","tau")
  
  print("This is Table 2")
  
  print("For tau = 0.2, 0.4, 0.8 from left to right")
  print("Gumbel Copula")
  print(summary)
  
  # Gaussian copula
  
  summaryI = round(summaryI,3)
  colnames(summaryI) <- c("Bias", "ESD", "RMSE", "Bias", "ESD", "RMSE", "Bias", "ESD", "RMSE")
  rownames(summaryI) <- c("beta1", "beta2", "eta0", "eta1", "eta2", "sigma","tau")
  
  print("Gaussian Copula")
  print(summaryI)
  
  # Make a cumulative hazard plot
  
  # Gumbel
  time = seq(0.1,10,by = 0.1)
  meanG = apply(Lhat,2,mean)           # average of the estimated cumulative hazard functions under the Gumbel copula
  
  # Gaussian copula 
  meanN = apply(Lhatn,2,mean)
  
  # True curve
  
  Lbd = (0.25*time^(3/4))

  #------------------------------
  # produce Figure 3, left panel
  #-----------------------------
  
  plot(time,meanG, type = "l",ylim = c(0,1.6),lty = 3,lwd = 2,col = "black", xlab = "Time", ylab = "Cumulative hazard", cex=3, cex.lab=1.5,cex.axis=1.5)
  title(main = "Scenario 2, Figure 3", line = 0.75)
  lines(time,Lbd,lty = 1,lwd = 2,col = "black", cex = 1.5)
  lines(time,meanN,lty = 2,lwd = 2,col = "grey",cex = 1.5)
  
  legend(3, 0.6, legend=c("Gumbel copula", "Gaussian copula","True function"), col=c("black","grey", "black"), lty=c(3,2,1),lwd = 2, cex=0.8,bty = "n",seg.len=1)
}


#============================================================================
# Scenario 3: To reproduce Table 3, Figure 2 and Figure 3 (right panel)
#============================================================================

# In this scenario, we examine how the sampling 
# distribution of our parameter estimators
# are close to a normal limit


Scenario3 = function(n,nsim,B,iseed,eps){
  
  # Specify the setting here
  
  cop = 4                                   # 4 = Gaussian copula
  dist = "lgn"                              # lgn = Lognormal distribution for dependent Censoring C  
  
  par = list(beta = c(0.45,1,0.25,0.5),eta = c(1,0.5,0.75), sd = c(1.5,0.75))  # True model parameters
  per = c(0,0) 
  
  results <- foreach(i = 1:nsim, .combine = 'comb', .multicombine = T, .maxcombine = 500, .inorder = F)%dopar%{ 
    
    source("data_simulation_functions.R")                 # load the functions with required libraries 
    source("estimation_algorithm_functions.R")
    
  data = DataGenerator(n, par, iseed+i,4,dist)            # generate data from Gaussian copula  
  
  data = data[order(data[,1]),] 
  Z = data[,1];
  Z[Z==0] = 1e-10
  d1 = data[,2];
  d2 = data[,3];
  X = data[,4:5];
  W = data[,6:8];
  per = per+c(table(d1)[2],table(d2)[2])
  init = c(0.4,0.9,0.85,0.35,0.6,1.2,0.5)
  output = NEIN(init,Z,d1,d2,X,W,4,dist,eps)
  temp = output$parhat
  tau = 2*asin(temp[7])/pi
  temp[7] = tau
  
  # Obtain cumulative hazards
  
  T1 = output$T1
  Lhat = output$cumhaz
  time = seq(0.1,10,by = 0.1)
  LD1 = rep(0,length(time))
  for (l in 1:length(time)){
    LD1[l] = Lhat[SearchIndicate(time[l],T1)]
  }
  
  # Bootstrap the data to estimate the standard errors 
  
  beta.star = matrix(NA,nrow = B,ncol = 7)
  for(b in 1:B)
  {
    set.seed(2143+b)
    samp1 = sample(length(data[,1]),replace = TRUE)    # resample the data with replacement
    datB = data[samp1,]
    datB = datB[order(datB[,1]),]
    Zb = datB[,1];
    db1 = datB[,2];
    db2 = datB[,3];
    Xb = datB[,4:5];
    Wb = datB[,6:8];
    initb = c(0.4,0.9,0.85,0.35,0.6,1.2,0.5)
    outputb = NEIN(initb,Zb,db1,db2,Xb,Wb,4,dist,eps)   # estimate the proposed model
    tb = outputb$parhat
    taub = 2*asin(tb[7])/pi
    tb[7] = taub
    beta.star[b,] = tb
  }
  Bootse = apply(beta.star,2,sd)     # bootstrap standard error
  
  # Confidence interval for tau based on Fisher's Z transformation 
  
  zt = 0.5*(log((1+temp[7])/(1-temp[7])))             
  se_z = (1/(1-temp[7]^2))*Bootse[7]
  zt_l = zt-1.96*(se_z)
  zt_u = zt+1.96*(se_z)
  
  # back transform
  
  r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)                # back transform to tau scale
  r_u = (exp(2*zt_u)-1)/(exp(2*zt_u)+1)
  EC1 = cbind(matrix(c(temp[1:6]-1.96*(Bootse[1:6]),r_l),ncol=1),matrix(c(temp[1:6]+1.96*(Bootse[1:6]),r_u), ncol=1))
  out = t(c(temp,Bootse,c(t(EC1))))

  list(out,LD1)      # save the output
  }
  
  #------------------
  # produce Table 3
  #------------------
  
  lognorm1 = results[[1]]
  par = list(beta = c(0.45,1,0.25,0.5),eta = c(1,0.5,0.75), sd = c(1.5,0.5398))  # true parameter values, tau = 0.54
  par0 = c(par[[1]][1:2],par[[2]],par[[3]])
  par0m = matrix(par0,nsim,7,byrow = TRUE)
  Bias = apply((lognorm1[,1:7]-par0m),2,mean)
  ESD = apply(lognorm1[,1:7],2,sd)
  BSE =  apply(lognorm1[,8:14],2,mean)
  RMSE = sqrt(apply((lognorm1[,1:7]-par0m)^2,2,mean))
  
  CP = rep(0,7)
  datacp = lognorm1[,15:28]
  for(i in 1:7){
    index = c(2*i-1,2*i)
    CP[i] = sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  summary = cbind(Bias,ESD,BSE,RMSE,CP)
  
  summary = round(summary,3)
  colnames(summary) <- c("Bias", "ESD", "BSE","RMSE", "CR")
  rownames(summary) <- c("beta1", "beta2", "eta0", "eta1", "eta2", "sigma","tau")
  
  print("This is Table 3")
  print("Gaussain copula")
  print(summary)
  
  
  #-------------------
  # Produce Figure 2 
  #-------------------
  
  # QQ plot for parameter estimates 
  
  par(mfrow = c(3,3))        
  
  print("This is Figure 2")
  
  qq = lognorm1[,1:7]
  qqnorm(qq[,1],ylim = c(0,1),cex.lab=1,cex.axis=1, main = expression(hat(beta)[1]))     # for beta1
  qqline(qq[,1],col = 1)
  qqnorm(qq[,2],ylim = c(0.6,1.3),cex.lab=1,cex.axis=1, main = expression(hat(beta)[2]))  # for beta2
  qqline(qq[,2],col = 1)
  qqnorm(qq[,3],cex.lab=1,cex.axis=1, main = expression(hat(eta)[0]))                      # for eta0
  qqline(qq[,3],col = 1)
  
  qqnorm(qq[,4],ylim = c(0,1),cex.lab=1,cex.axis=1, main = expression(hat(eta)[1]))         # for eta1
  qqline(qq[,4],col = 1)
  qqnorm(qq[,5],ylim = c(0.2,1.2), cex.lab=1,cex.axis=1, main = expression(hat(eta)[2]))     # for eta2
  qqline(qq[,5],col = 1)
  qqnorm(qq[,6],cex.lab=1,cex.axis=1, main = expression(hat(sigma)))      # for sigma
  qqline(qq[,6],col = 1)
  qqnorm(qq[,7],cex.lab=1,cex.axis=1, main = expression(hat(tau)))        # for tau
  qqline(qq[,7],col = 1)
  
  # Fisher's Z transform
  
  zt <- 0.5 * (log((1 + qq[,7])/(1 - qq[,7])))      # Fisher's z transform for rho
  
  qqnorm(zt, cex.lab=1,cex.axis=1, main = expression(hat(omega)))     
  qqline(zt,col=1)
  
  
  #------------------------------
  # produce Figure 3, right panel
  #--------------------------------
  
  # Estimated cumulative hazards
  
  cum3 = results[[2]]
  mean3 = apply(cum3,2,mean)
  mean3 = apply(cum3,2,mean)
  
  # True curve
  time = seq(0.1,10,by = 0.1)
  Lbd = (0.25*time^(0.5))
  
  # Make a plot
  
  plot(time,mean3, type = "l",lty = 3,lwd = 3,col = "black", xlab = "Time", ylab = "Cumulative hazard", cex=3, cex.lab=1.5,cex.axis=1.5)
  lines(time,Lbd,lty = 1,lwd = 3,col = "black", cex = 1.5)
  title(main = "Scenario 3, right panel", line = 0.75)
  legend(3, 0.4, legend=c("Gaussian copula","True function"), col=c("black","black"), lwd = 2, lty=c(3,1), cex=0.8,bty = "n",seg.len=1)
}



#===========================================================================
# Scenario 4, To reproduce Table 6, see the Supplementary Material
#===========================================================================

# Under this scenario, we like to know whether the estimated 
# parameters for the survival model are sensitive to 
# model misspecification. 


# We consider two cases:
# Case 1: Marginal distribution of C is mispecified
# Case 2: Both the margin of C and the copula structure
# are misspecified. 



# The arguments for Scenario4: 
                          # n = sample size
                          # nsim = number of simulations
                          # B = bootstrap size
                          # iseed = starting value to generate random number
                          # eps = convergence error


Scenario4 = function(n,nsim,B,iseed,eps){
  
  cop = 2                                   # 2 = Frank, 3 = Gumbel, 4 = normal
  dist = "wb"                               # wb = Weibull, lgn = Lognormal
  dist1 = "lgn"
  mod = c(2,3)                              # 2 = Frank, 3 = Gumbel

  par = list(beta = c(0.45,1,0.25,0.5),eta = c(1,0.5,0.75), sd = c(1.5,7))   # True parameters
  
  # parallel 
  
  results <- foreach(i = 1:nsim, .combine = 'comb', .multicombine = T, .maxcombine = 500, .inorder = F)%dopar%{ 
    
    source("data_simulation_functions.R")                 # load the functions with required libraries 
    source("estimation_algorithm_functions.R")
    
  data = DataGenerator(n, par,iseed,cop,dist1)             # generate data from a Frank copula with lognormal margin 
  data = data[order(data[,1]),] 
  Z = data[,1];
  Z[Z==0] = 1e-10
  d1 = data[,2];
  d2 = data[,3];
  X = data[,4:5];
  W = data[,6:8];
  
  #  mis-specified margin or mis-specified margin + copula

  for (j in mod){  # j = 2, mispecified margin of C and j = 3, misspecified margin of C and copula
    
    init = c(0.4,0.9,0.85,0.35,0.6,1.2,3)
    output = NEIN(init,Z,d1,d2,X,W,j,dist,eps)
    temp = output$parhat
    
    # Obtain cumulative hazards
    
    T1 = output$T1
    Lhat = output$cumhaz
    time = seq(0.1,10,by = 0.1)
    LD1 = rep(0,length(time))
    for (l in 1:length(time)){
      LD1[l] = Lhat[SearchIndicate(time[l],T1)]
    }
    
    # Bootstrap standard error and obtain confidence interval
    
  
    beta.star = matrix(NA,nrow = B,ncol = 7)
    for(b in 1:B)
    {
      set.seed(2143+b)
      samp1 = sample(length(data[,1]),replace = TRUE)
      datB = data[samp1,]
      datB = datB[order(datB[,1]),]
      Zb = datB[,1];
      db1 = datB[,2];
      db2 = datB[,3];
      Xb = datB[,4:5];
      Wb = datB[,6:8];
      initb = c(0.4,0.9,0.85,0.35,0.6,1.2,3)
      outputb = NEIN(initb,Zb,db1,db2,Xb,Wb,j,dist,eps)
      tb = outputb$parhat
      
      if(j == 2){
        taub = tau(frankCopula(tb[7]))
        tb[7] = taub
      }
      if(j == 3){
        taub = tau(gumbelCopula(tb[7]))
        tb[7] = taub
      }
      beta.star[b,] = tb
    }
    Bootse = apply(beta.star,2,sd)
    
    # Confidence interval for tau
    
    if (j == 2){
      tau = tau(frankCopula(temp[7]))
      temp[7] = tau 
      zt = 0.5*(log((1+temp[7])/(1-temp[7])))             # Fisher's z transform
      se_z = (1/(1-temp[7]^2))*Bootse[7]
    }
    if (j == 3){
      tau = tau(gumbelCopula(temp[7]))
      temp[7] = tau
      zt = 0.5*(log((1+temp[7])/(1-temp[7])))             # Fisher's z transform
      se_z = (1/(1-temp[7]^2))*Bootse[7]
    }
    zt_l = zt-1.96*(se_z)
    zt_u = zt+1.96*(se_z)
    
    # back transform
    r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)                  # back transform to tau scale
    r_u = (exp(2*zt_u)-1)/(exp(2*zt_u)+1)
    EC1 = cbind(matrix(c(temp[1:6]-1.96*(Bootse[1:6]),r_l),ncol=1),matrix(c(temp[1:6]+1.96*(Bootse[1:6]),r_u), ncol=1))
    
    if (j == 2){
      parhatF = t(c(temp,Bootse,c(t(EC1))))
      LambdaF = LD1
     } 
    if (j == 3){
      parhatG = t(c(temp,Bootse,c(t(EC1))))
      LambdaG = LD1
    }
  }
  list(parhatF,parhatG,LambdaF,LambdaG)      # save the output
  }
  
  # Case 1, misspecified margin as Weibull
  
  weib = results[[1]]        # retrieve parameter estimates
  par = list(beta = c(0.45,1,0.25,0.5),eta = c(1,0.5,0.75), sd = c(1.5,0.56226))   # True parameters
  par0 = c(par[[1]][1:2],par[[2]],par[[3]])
  par0m = matrix(par0,nsim,7,byrow=TRUE)
  Bias = apply((weib[,1:7]-par0m),2,mean)
  ESE = apply(weib[,1:7],2,sd)
  BSE =  apply(weib[,8:14],2,mean)
  RMSE = sqrt(apply((weib[,1:7]-par0m)^2,2,mean))
  
  CP = rep(0,7)
  datacp = weib[,15:28]
  for(i in 1:7){
    index = c(2*i-1,2*i)
    CP[i] = sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  sum = cbind(Bias,ESE,BSE,RMSE,CP)
  
  #retain estimates corresponding to survival model
  summary = sum[c(1,2,7),]
  summary = round(summary,3)
  colnames(summary) <- c("Bias", "ESD", "BSE","RMSE", "CR")
  rownames(summary) <- c("beta1", "beta2", "tau")

  #------------------------------------------------
  # produce Table 6, Case 1, Supplementary Material
  #------------------------------------------------
  
  print("This is Table 6, Case 1")
  print("Misspecified Margin for C")
  print(summary)
  
  
  # Case 2, misspecified margin  and copula
  
  weibcop = results[[2]]                       # retrieve parameter estimates
  par0 = c(par[[1]][1:2],par[[2]],par[[3]])
  par0m = matrix(par0,nsim,7,byrow=TRUE)
  Bias = apply((weibcop[,1:7]-par0m),2,mean)
  ESE = apply(weibcop[,1:7],2,sd)
  BSE =  apply(weibcop[,8:14],2,mean)
  RMSE = sqrt(apply((weibcop[,1:7]-par0m)^2,2,mean))
  
  CP = rep(0,7)
  datacp = weibcop[,15:28]
  for(i in 1:7){
    index = c(2*i-1,2*i)
    CP[i] = sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  sum = cbind(Bias,ESE,BSE,RMSE,CP)
  
  sum = cbind(Bias,ESE,BSE,RMSE,CP)
  
  #retain estimates corresponding to survival model
  
  summary = sum[c(1,2,7),]
  summary = round(summary,3)
  colnames(summary) <- c("Bias", "ESD", "BSE","RMSE", "CR")
  rownames(summary) <- c("beta1", "beta2", "tau")
  
  #------------------------------------------------
  # produce Table 6, Case 2, Supplementary Material
  #------------------------------------------------
  
  print("This is Table 6 , Case 2")
  print("Misspecified Margin for C and Copula")
  print(summary)
  
  
  # Cumulative hazard plots
  
  # Scenario 4, Case 1
 
  cum4 = results[[3]]           # estimated cumulative hazard functions under case 1
  mean4 = apply(cum4,2,mean)    # average of the cumulative hazard
   
  # Sceenario 4, Case 2
  
  cum42 = results[[4]]            # estimated cumulative hazard functions under case 2
  mean42 = apply(cum42,2,mean)    # average of the cumulative hazard
  
  # True curve
  time = seq(0.1,10,by = 0.1)
  Lbd = (0.25*time^(0.5))
  
  #--------------------------------------------------
  # produce Figure 4, see the Supplementary Material
  #--------------------------------------------------
  
  par(mfrow=c(1,1)) 
  plot(time,mean42, type = "l",lty = 3,lwd = 3,ylim = c(0,1),col = "black", xlab = "Time", ylab = "Cumulative hazard", cex=3, cex.lab=1.5,cex.axis=1.5)
  lines(time,Lbd,lty = 1,lwd = 3,col = "black", cex = 1.5)
  lines(time,mean4,lty = 2,lwd = 3,col = "black", cex = 1.5)
  
  legend(6, 0.4, legend=c("Case 1","Case 2", "True function"), lty=c(2,3,1), cex= 0.8,lwd = 2,bty = "n",seg.len=1.5)
}




#==========================================================================
# Scenario 5,Table 7, see the Supplementary Material
#==========================================================================

# Under this scenario, we assess the performance 
# of the goodness-of-fit test using
# simulations


# We test a null hypothesis under the following  three cases:
#### Case 1: Correctly specified model
#### Case 2: Misspecified margin of C
#### Case 3: Misspecified regression functions

# This scenario is very time consuming
# since it uses bootstrapping to determine
# the null distribution of our test 
# statistic


Scenario5 = function(n,nsim,B,iseed,eps){          # Compute GOF
  
  # Specify the setting here
  
  gm = 14.14                                       # copula parameter
  cop = 2                                          # 2 = Frank, 3 = Gumbel, 4 = Gaussian 
  dist = "lgn"                                     # lgn = lognormal
  case1 = "lgn"; case2 = "wb"; case3 = "misrf";    # lgn = lognormal, wb = Weibull, misrf = misspecified regression function
        
   results <- foreach(i = 1:nsim, .combine = 'comb', .multicombine = T, .maxcombine = 500, .inorder = F)%dopar%{ 
    
    source("data_simulation_functions.R")                           # load the functions with required libraries 
    source("estimation_algorithm_functions.R")
  
        par = list(beta = c(0.5,1,0.25,0.5),eta = c(1,0.35,0.75), sd = c(2,gm))
        data = DataGenerator(n, par, iseed+i,cop,dist)        # generate data from frank copula  
        data = data[order(data[,1]),] 
        Z = data[,1];
        Z[Z==0] = 1e-10
        d1 = data[,2];
        d2 = data[,3];
        d3 = d1+d2
        X = data[,4:5];
        W = data[,6:8];
        
        # Case 1, fit the correct model
        
        init = c(0.4,0.7,0.7,0.25,0.8,0.5,10)
        output1 = NEIN(init,Z,d1,d2,X,W,cop,case1,eps)
        parhat1 = output1$parhat          
        T1 = output1$T1
        Lhat1 = output1$cumhaz
        
        # conduct the test under case 1
        
        results1 = CMtest(parhat1,Z,d3,Lhat1,T1,X,W,cop,case1,B,eps)  
        
        # Case 2, fit the misspecified censoring model
        
        output2 = NEIN(init,Z,d1,d2,X,W,cop,case2,eps)
        parhat2 = output2$parhat          
        T1 = output2$T1
        Lhat2 = output2$cumhaz
        
        # conduct the test
        
        results2 = CMtest(parhat2,Z,d3,Lhat2,T1,X,W,cop,case2,B,eps)  
        
        # GOF under case 3
        par = list(beta = c(0.5,1,0.5,0.25,0.5),eta = c(1,0.35,0.75,-0.6), sd = c(2,gm))  # True parameters
        data = dat.covmis(n, par, iseed+i,cop,dist)                                         # generate data from misspecified regression functions 
        data = data[order(data[,1]),] 
        Z = data[,1];
        Z[Z==0] = 1e-10;
        d1 = data[,2];
        d2 = data[,3];
        d3 = d1+d2;
        X = data[,4:5];
        W = data[,7:9];
        
        # Fit the misspecified regression functions
        
        init = c(0.3,0.7,0.7,0.25,0.8,0.5,10)
        output3 = NEIN(init,Z,d1,d2,X,W,cop,dist,eps)
        parhat3 = output3$parhat          
        T1 = output3$T1
        Lhat3 = output3$cumhaz
        
        # Conduct the test under case 3
        results3 = CMtest(parhat3,Z,d3,Lhat3,T1,X,W,cop,dist,B,eps)
        
        list(results1,results2,results3)    # save the output
   }
  
  
  # Obtain results
  case1 = results[[1]] 
  case2 = results[[2]]
  case3 = results[[3]]
  
  sum1 = apply(case1,2,mean)       # average the number of counts
  sum2 = apply(case2,2,mean)
  sum3 = apply(case3,2,mean)
  
  summary = rbind(sum1,sum2,sum3)
  summary = round(summary,3)
  colnames(summary) <- c("Expected T_CM", "Bootstrap T_CM", "10%","5%")
  rownames(summary) <- c("Case 1", "Case 2", "Case 3")
  
  print("This is Table 7")
  print(paste0("The sample size is ", n))
  print(summary)                                    # output the GOF results
} 





###  End


    