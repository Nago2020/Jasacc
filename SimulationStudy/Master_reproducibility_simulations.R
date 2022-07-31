
#===================================================================
# Copula based Cox proportional hazards
# models for dependent censoring
# [Authors and affiliation blinded]
# July, 202
#===================================================================


#===================================================================
# Reproducibility: MASTER script simulation study
#===================================================================

# This file contains instructions for reproducing all the simulated 
# data analyses, including tables and figures contained in the paper. 

# Download the repository from Github [url blinded].
# The following script assumes the working directory has
# been set to the simulations sub-folder in this folder. 
# All of the necessary code to generate data and to analyze
# the simulated data are now available.



# As denoted below, Scenarios 3, 4 and 5 are computationally 
# intensive and take a *very* long time to run on personal computer, 
# even after applying parallelization in R.  



#======================================================================
# Install necessary packages as well as the R version used for
# all analyses
#======================================================================

# > sessionInfo()
# R version 3.6.0 (2019-04-26)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)



# Install packages

library(copula)
library(survival)
library(pbivnorm)

# parallelization  

library(parallel)
library(foreach)
library(doParallel)

numCores = detectCores()
numCores = numCores-1
clust = makeCluster(numCores)
registerDoParallel(clust)



#========================================================================
# simulation_results.R 
# Master script to conduct simulations in the paper.
# The functions in this script used to reproduce all the tables 
# and figures. 
#========================================================================


source("simulation_results.R")       # load functions


# We follow the following five steps to reproduce the simulation
# study section


## !!! NOTE: Running the following five functions in the simulations_results.R
# will take EXTREMELY LONG time !!!

# We ran all of our simulations using parallel computing on 
# the laboratory supercomputer (the necessary code
# for doing so is omitted for brevity).



#=============================================================================
# Step 1: Scenario1()
# This function is used to reproduce Table 1 and Figure 1
# in the paper.
#=============================================================================


# specify the setting here

n = 500                       # sample size
nsim = 1000                   # number of data created/ number of simulations
iseed = 684132                # starting value to generate random number
eps = 1e-4                    # convergence error



# Run simulations; this may take very long time to see results



Scenario1(n,nsim,iseed,eps)





#=============================================================================
# Step 2: Scenario2()
# This function is used to reproduce Table 2 (in the paper) and Figure 3,
# left panel, (in the Supplementary Material)
#=============================================================================


# specify the setting here

n = 500                       # sample size
nsim = 1000                   # number of data created/ number of simulations
iseed = 684132                # starting value to generate random number
eps = 1e-4                    # convergence error



# Run simulations; this may take very long time on 
# personal laptop


Scenario2(n,nsim,iseed,eps)





#=============================================================================
# Step 3: Scenario3()
# This function is used to reproduce Table 3 (in the paper) and Figure 3,
# right panel, (in the Supplementary Material)
#=============================================================================


# specify the setting here

n = 500                       # sample size
nsim = 500                    # number of data created/ number of simulations
B =  150                      # Bootstrap size
iseed = 684132                # starting value to generate random number
eps = 1e-4                    # convergence error



# Run simulations; this is computationally very time consuming
# since it uses bootstrapping to obtain bootstrap 
# standard errors



Scenario3(n,nsim,B,iseed,eps)






#=============================================================================
# Step 4: Scenario4()
# This function is used to reproduce Table 6 and Figure 4 
# (in the Supplementary Material)
#=============================================================================


# specify the setting here

n = 500                       # sample size
nsim = 500                    # number of data created/ number of simulations
B =  150                      # Bootstrap size
iseed = 684132                # starting value to generate random number
eps = 1e-4                    # convergence error



# Run simulations; this is computationally very time consuming
# since it uses bootstrapping to obtain bootstrap 
# standard errors


Scenario4(n,nsim,B,iseed,eps)




#=============================================================================
# Step 5: Scenario5()
# This function is used to reproduce Table 7 
# (in the Supplementary Material)
#=============================================================================


# specify the setting here

nn = c(500,1000)              # sample size, n = 500 and 1000
nsim = 500                    # number of data created/ number of simulations
B =  250                      # Bootstrap size
iseed = 684132                # starting value to generate random number
eps = 1e-4                    # convergence error



# Run simulations; this is extremely time consuming
# since it uses bootstrapping to determine the null  
# distribution of the goodness-of-fit test statistic

for (n in nn){
  Scenario5(n,nsim,B,iseed,eps)       # print Table 7
}



stopCluster(clust)

#### End 

