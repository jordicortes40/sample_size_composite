##################################################################################
#                    Calculate HRs
##################################################################################
# Computation of the Combined HR (t) for several copulas and marginal distributions.
# Two different ways of input are allowed: by means of probabilities (Taking into 
# consideration Cases 1,2,3,4); and inputing the scale parameters b10, b20.
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
library(cowplot)
library(ggplot2)
library(gridExtra)
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
source('../functions/simulations cause-specific/functions_marginal_selection_modified.R')
source('../functions/functions_simulations.R')
source('../functions/functions_ARE.R')
source('../parameters/parameters.R')


################################################################################## 
# Scenarios parameters
##################################################################################
HR1   <- 0.95 
HR2   <- 0.35 
P1    <- 0.14 
P2    <- 0.05 
BETA1 <- c(0.70,0.5,1,2)
BETA2 <- c(0.91,0.5,1,2)
RHO   <- c(0.1,0.5)
CASE  <- 4
COPULA<- "Frank"

############################################################
# Settings
############################################################
# d0: all combinations
# d1: to store results
d0 <- expand.grid(copula=COPULA,HR2=HR2,HR1=HR1,rho=RHO,
                  p2=P2,p1=P1,beta2=BETA2,beta1=BETA1,case=CASE,stringsAsFactors = FALSE)[,9:1]
d0 <- d0[d0$beta1==0.70 & d0$beta2==0.91   & (d0$rho==0.1 | d0$rho==0.5),]
d1 <- matrix(as.numeric(NA),ncol=length(VAR.NAMES),nrow=nrow(d0))


############################################################
# Estimate HR composite
############################################################
t_values <- c(0.0001,seq(0.001,1,0.001))  # times to assess
t0 <- Sys.time()
for (i in 1:nrow(d0)) {
  ##-- Estimate summary measures of HR*
  COMBINED_HR_values <- COMBINED_HR(t=t_values,
                                    rho0=d0$rho[i],
                                    beta1=d0$beta1[i],beta2=d0$beta2[i],
                                    HR1=d0$HR1[i],HR2=d0$HR2[i],
                                    p1=d0$p1[i],p2=d0$p2[i],
                                    case=d0$case[i],copula=d0$copula[i],
                                    tau=seq(0.25,1,0.25))
  ##-- Store information
  d1[i,] <- as.numeric(unlist(COMBINED_HR_values)) 
  
  ##-- Print execution time
  print_time(t0,i,nrow(d0))
}
d <- cbind(d0,d1)
names(d) <- c(names(d0),VAR.NAMES)
View(cbind(d[,1:9],data.frame(gAHR=d$gAHR_tau_100)),'d')

##################################################################################
#                    Power simulation composite endpoints
##################################################################################

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
dd <- subset(d,!is.na(d$gAHR_tau_100) & d$gAHR_tau_100!=0)
dim(dd)

# Parameters
n.SETTINGS <- nrow(dd)     # Number of settings
nsim  <- 10000             # number of simulations
alpha <- 0.05              # Type I error  
power <- 0.8               # Power

# Parallel computation
sfInit(parallel=TRUE, cpus=detectCores()-1)
sfLibrary(survival)


# Store information
power_gAHR <- Events_mean <- c()

# Simulations
t0 <- Sys.time()
ite.ini <- 1
set.seed(12345)

for(j in ite.ini:n.SETTINGS){
  
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
  
  ##-- Nmax is the maximum sample size in each group from the 3 options
  ##-- It is need in order to simulate the data 
  n_gAHR <- ceiling(dd$N_gAHR[j]/2)
  Nmax <- n_gAHR
  if(Nmax>10000){nsim_aux <- nsim; nsim <- 100} # 
  
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
  power_gAHR[j] <- sum(res_gAHR<0.05,na.rm=TRUE)/nsim

  
  Events_mean[j] <- mean(apply(TC_gAHR,2,sum_events,tau=1))
  ##-- Time ----------------------------------------------------
  print_time(t0,j,n.SETTINGS)

  if(Nmax>10000){nsim <- nsim_aux}
  
  cat('power_log-rank:',fc(power_gAHR[j]), #'power_ARE:',fc(power_ARE[j]),
      'Events needed:',dd$e_gAHR[j],'Events mean:',Events_mean[j],'\n')

}
sfStop()

##-- Add results to data.frame
dd$power_gAHR <- power_gAHR


##-- Remove big objects
sort(sapply(sapply(ls(),get),object.size))
rm(TC_gAHR,TC,BI1,BI0,TC0,TC1,T10,T11,T20,T21)

##-- Save image
try(setwd('C:/Users/jcortes/Dropbox/Sample Size T2E/Simulations/data'))
try(setwd('C:/Users/jordi/Dropbox/Sample Size T2E/Simulations/data'))
save(dd,file=paste0('../outputs/ARREST_4.RData'))




