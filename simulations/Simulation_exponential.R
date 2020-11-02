##################################################################################
#                    Power simulation composite endpoints
##################################################################################
# 
# 
# 
#
################################################################################## 
#
# Last update: 10/03/2020
#
# R version: R 3.6.1
#
# Authors: Jordi Cortes (jordi.cortes-martinez@upc.edu)
# 
#
################################################################################## 
#
# References:
#   - Gomez G. and Lagakos S. Statistical considerations when using a composite 
#     endpoint for comparing treatment groups. Statistics in Medicine, 2013.
#
##################################################################################

rm(list=ls())

################################################################################## 
# Load packages
##################################################################################
library(copula)
library(numDeriv)
library(reshape)
library(rootSolve)
library(data.table)
library(survival)
library(snowfall)
library(snow)
library(BivarP)
library(parallel)

################################################################################## 
# Load functions and parameters
##################################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('../functions/functions.R')
source('../functions/functions_simulations.R')
source('../functions/functions_ARE.R')
source('../parameters/parameters.R')


################################################################################## 
# Load data
##################################################################################
load('../outputs/HRs_2020-11-02.RData')

############################################################
# Sample sizes with different values
############################################################
##-- Events gAHR and AHR
d$e_gAHR <- schoendfeld.formula(alpha=0.05,power=0.8,HR = d$gAHR_tau_100)
d$e_AHR  <- schoendfeld.formula(alpha=0.05,power=0.8,HR = d$AHR_tau_100)

##-- Sample sizes gAHR and AHR
d$N_gAHR<-  2*ceiling(d$e_gAHR/(d$pstar_0 + d$pstar_1)) 
d$N_AHR <-  2*ceiling(d$e_AHR /(d$pstar_0 + d$pstar_1))


############################################################
# Simulations to know the power
############################################################
## Subset
dd <- subset(d,copula=='Frank' & beta1==1 & beta2==1 & HR1>=0.6 & HR2>=0.6)
dim(dd)

# Parameters
n.SETTINGS <- nrow(dd)     # Number of settings
nsim  <- 10000             # number of simulations
alpha <- 0.05              # Type I error --> https://www.reddit.com/r/statistics/comments/9xx5ni/onesided_log_rank_test/  
power <- 0.8               # Power

# Parallel computation
sfInit(parallel=TRUE, cpus=detectCores()-1)
sfLibrary(survival)


# Store information
dd$power_gAHR <- dd$Events_mean <- NA

# Simulations
t0 <- Sys.time()
ite.ini <- 1
set.seed(12345)

for(j in ite.ini:n.SETTINGS){
  
  n_gAHR <- ceiling(dd$N_gAHR[j]/2)
  Nmax <- n_gAHR

  if(Nmax<10000){
    
    p10 <- dd$p1[j]        # Probability in reference group for relevant endpoint
    HR1 <- dd$HR1[j]       # HR for relevant endpoint
    p20 <- dd$p2[j]        # Probability in reference group for aditional endpoint
    HR2 <- dd$HR2[j]       # HR for relevant endpoint
    
    beta1 <- dd$beta1[j]   # Shape parameter for relevant endpoint
    beta2 <- dd$beta2[j]   # Shape parameter for relevant endpoint
    
    copula <- dd$copula[j]
    rho <- dd$rho[j]
    case <- dd$case[j]
    
    ###################################################
    ##-- Estimate parameters for distributions
    ###################################################
    ##-- Find parameters
    theta <- CopulaSelection(copula=copula,rho=rho)[[2]]
    
    par.shape <- as.numeric(dd[j,c("beta1","beta2","beta1","beta2")])            # Weibull shape parameters
    par.scale <- as.numeric(dd[j,paste0('scale_',c("b10","b20","b11","b21"))])   # Weibull scale parameters
 
    ##-- Select copula
    if(copula=='Frank')  cop <- frankCopula(param=theta, dim = 2)
    if(copula=='Clayton')cop <- claytonCopula(param=theta, dim = 2)
    if(copula=='Gumbel') cop <- gumbelCopula(param=theta, dim = 2)
    
    ##-- Controls
    MVDC0 <- mvdc(copula = cop, 
                  margins = c('weibull','weibull'), 
                  paramMargins =  list(list(shape = par.shape[1], scale = par.scale[1]), 
                                       list(shape = par.shape[2], scale = par.scale[2])), 
                  marginsIdentical = FALSE,
                  check = TRUE, fixupNames = TRUE)
    BI0 <- rMvdc(Nmax*nsim, MVDC0)
    T10 <- matrix(BI0[,1],ncol=nsim,byrow=FALSE)
    T20 <- matrix(BI0[,2],ncol=nsim,byrow=FALSE)
    
    ##-- Treateds
    MVDC1 <- mvdc(copula = cop, 
                  margins = c('weibull','weibull'), 
                  paramMargins =  list(list(shape = par.shape[3], scale = par.scale[3]), 
                                       list(shape = par.shape[4], scale = par.scale[4])), 
                  marginsIdentical = FALSE,
                  check = TRUE, fixupNames = TRUE)
    BI1 <- rMvdc(Nmax*nsim, MVDC1)
    T11 <- matrix(BI1[,1],ncol=nsim,byrow=FALSE)
    T21 <- matrix(BI1[,2],ncol=nsim,byrow=FALSE)

    ##-- Time for composite
    TC0 <- pmin(T10,T20)
    TC1 <- pmin(T11,T21)
    TC <- rbind(TC0,TC1)
    
    ##-- n_gAHR
    Tx_gAHR <- c(rep(0,n_gAHR),rep(1,n_gAHR))
    TC_gAHR <- rbind(TC0[1:n_gAHR,],TC1[1:n_gAHR,])
    
    #-- Power log-rank
    res_gAHR <- sfApply(TC_gAHR,2,f.log.rank,tau=1,Treatment=Tx_gAHR)
    dd$power_gAHR[j] <- sum(res_gAHR<0.05,na.rm=TRUE)/nsim
    
    dd$Events_mean[j] <- mean(apply(TC_gAHR,2,sum_events,tau=1))
    
    ##-- Time ----------------------------------------------------
    print_time(t0,j,n.SETTINGS)
    # cat('power_log-rank:',fc(dd$power_gAHR[j]), #'power_ARE:',fc(power_ARE[j]),
    #     'Events needed:',dd$e_gAHR[j],'Events mean:',dd$Events_mean[j],
    #     'parameters:',paste0(as.character(dd[j,1:9]),collapse = ','),'\n')
    gc()
  }
  # Store results
  if(j %% 100==0 | j==n.SETTINGS){
    write.table(dd,file='data/simulations_exponential_case3_and_4.txt',
                sep='\t',append=FALSE,row.names = FALSE,col.names=TRUE,quote = FALSE)
  }
}
sfStop()

##-- Remove big objects
sort(sapply(sapply(ls(),get),object.size))
rm(TC_gAHR,TC,BI1,BI0,TC0,TC1,T10,T11,T20,T21)

##-- Save image
save.image(file='../outputs/simulations_exponential_case3_and_4.RData')









