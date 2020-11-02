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


################################################################################## 
# Load functions and parameters
##################################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('../functions/functions.R')
source('../functions/functions_simulations.R')
source('../functions/functions_ARE.R')
source('../parameters/parameters.R')

################################################################################## 
# Scenarios parameters
##################################################################################
HR1   <- c(0.2,0.4,0.6,0.8,0.9)
HR2   <- c(0.2,0.4,0.6,0.8,0.9)
P1    <- c(0.05,0.1,0.3,0.5)
P2    <- c(0.05,0.1,0.3,0.5)
BETA1 <- c(0.5,1,2)
BETA2 <- c(0.5,1,2)
RHO   <- c(0.1,0.3,0.5,0.8)
CASE  <- c(3,4)
COPULA<- c("Frank","Clayton","Gumbel")

############################################################
# Settings
############################################################
# d0: all combinations
# d1: to store results
d0 <- expand.grid(case=CASE,copula=COPULA,HR2=HR2,HR1=HR1,rho=RHO,
                  p2=P2,p1=P1,beta2=BETA2,beta1=BETA1,stringsAsFactors = FALSE)[,c(1,9:2)]
d0 <- subset(d0,!(case==4 & p1==0.5 & p2==0.5))
d1 <- matrix(as.numeric(NA),ncol=length(VAR.NAMES),nrow=nrow(d0))


############################################################
# Estimate HR composite
############################################################
t_values <- c(0.0001,seq(0.001,1,0.001))  # times to assess
t0 <- Sys.time()
for (i in 1:nrow(d0)) {
  ##-- Estimate summary measures of HR*
  COMBINED_HR_values <- COMBINED_HR(t=t_values,                          # time values
                                    rho0=d0$rho[i],                      # correlation
                                    beta1=d0$beta1[i],beta2=d0$beta2[i], # shape parameters of weibull distribution
                                    HR1=d0$HR1[i],HR2=d0$HR2[i],         # Hazard ratios
                                    p1=d0$p1[i],p2=d0$p2[i],             # Probabilities of observing the event
                                    case=d0$case[i],copula=d0$copula[i], # CASE (1,2,3,4) according to the paper og Gomez and Lagakos (2013) and copula name
                                    tau=seq(0.25,1,0.25))                # Tiem of follow-up to assess tha average hazard ratio measures
  ##-- Store information
  d1[i,] <- as.numeric(unlist(COMBINED_HR_values)) 
  
  ##-- Print execution time
  print_time(t0,i,nrow(d0))
}
d <- cbind(d0,d1)
names(d) <- c(names(d0),VAR.NAMES)
View(head(d),'d')

save(d,file=paste0('../outputs/HRs_',as.character(Sys.Date()),'.RData'))






