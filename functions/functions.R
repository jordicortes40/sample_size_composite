##-- Sample sizes using SHOENDFELD formula
schoendfeld.formula <- function(alpha,power,HR) E <- 4*(qnorm(1-alpha/2) +  qnorm(power))^2 / (log(HR))^2

#######################################################################################
# Functions: several copula densities (Table 3.3 of Trivedi, 2007: An Introduction for Practitioners) 
#
#######################################################################################
# Description: It computes the Combined HR(t) (HR_Star) funtion 
# 
# u	      value between 0 and 1 for endpoint 1 u=F(t_1)
# v	      value between 0 and 1 for endpoint 2 v=F(t_2)
# theta	  Copula parameter
#######################################################################################

dFrank <- function(u,v,theta){(theta*(1-exp(-theta))*exp(-theta*(u+v)))/ (exp(-theta) +  exp(-theta*(u+v))- exp(-theta*u)-exp(-theta*v))^2}
dClayton <- function(u,v,theta){ (u*v)^(-theta-1) * (theta+1) * (u^(-theta) + v^(-theta) - 1)^(-2 - 1/theta)}
dGumbel <- function(u,v,theta){
  u1 <- -log(u)
  u2 <- -log(v)
  num1 <- exp(-(u1^theta + u2^theta)^(1/theta))
  num2 <- (u*v)^(-1)
  num3 <- (u1*u2)^(theta-1)
  num4 <- (u1^theta + u2^theta)^(1/theta) + theta - 1
  num <- num1 * num2 * num3 * num4
  den <- (u1^theta + u2^theta)^(2-1/theta)
  num/den
}
dFGM <- function(u,v,theta) 1 + theta * (1- 2*u) * (1 - 2*v)
dNormal <- function(u,v,theta){
  x <- qnorm(u)
  y <- qnorm(v)
  
  (1 - theta^2)^(-1/2) * exp(-(x^2 + y^2 - 2*theta*x*y)/(2*(1-theta^2))) * exp((x^2 + y^2)/2)
}
#######################################################################################
# Functions related with function COMBINED_HR 
#
#######################################################################################

#**************************************************************************************
##-- Case 1 ***************************************************************************
#**************************************************************************************

#######################################################################################
# Function: Sstar 
#
#######################################################################################
# Description: Returns the value of the survival function of S* at point x given the
#              marginal distributions and the bivariate distributions via copula 
#
# x	        Point in which to be evaluated
# dist1     Distribution function of the marginal T1 (pweibull) 
# dist2     Distribution function of the marginal T2 (pweibull) 
# param1    Parameters of the marginal distribution function T1 (pweibull) 
# param2    Parameters of the marginal distribution function T2 (pweibull) 
# dist_biv  Distribution function of the bivariate distribution via copula
#######################################################################################

Sstar <- function(x,dist1,dist2,param1,param2,dist_biv) { 
  y <- if(length(x) == 1) c(x,x) else cbind(x,x)
  
  
  return(
    1
    - do.call(dist1,c(list(q=x),param1))
    - do.call(dist2,c(list(q=x),param2))
    + (pMvdc(y, dist_biv))
  )
}

#**************************************************************************************
##-- Case 1 ***************************************************************************
#**************************************************************************************
Sstar0_func <- function(z,...) Sstar(x=z,...)
Sstar1_func <- function(z,...) Sstar(x=z,...)
fstar0_func <- function(z,T1pdist,T2pdist,T10param,T20param,distribution0) -grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
fstar1_func <- function(z,T1pdist,T2pdist,T11param,T21param,distribution1) -grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)
Sstar0_func_perc <- function(z,perc,...) Sstar(x=z,...) - perc
Sstar1_func_perc <- function(z,perc,...) Sstar(x=z,...) - perc

#**************************************************************************************
##-- Case 3 ***************************************************************************
#**************************************************************************************
# Sstar0_func_case3 <- function(z,theta,b10,b20,beta1,beta2,ST10,ST20){
#   ST10 <- exp(-(z/b10)^beta1)
#   ST20 <- exp(-(z/b20)^beta2)
#   (-log(1+(exp(-theta*ST10)-1)*(exp(-theta*ST20)-1)/(exp(-theta)-1))/theta)
# }
# Sstar1_func_case3 <- function(z,theta,b11,b21,beta1,beta2,ST11,ST21){
#   ST11 <- exp(-(z/b11)^beta1)
#   ST21 <- exp(-(z/b21)^beta2)
#   (-log(1+(exp(-theta*ST11)-1)*(exp(-theta*ST21)-1)/(exp(-theta)-1))/theta)
# }
# Sstar0_func_perc_case3 <- function(z,perc,theta,...) Sstar0_func_case3(z,theta,...) - perc
# Sstar1_func_perc_case3 <- function(z,perc,theta,...) Sstar1_func_case3(z,theta,...) - perc

#######################################################################################
# Function: HRstar_func, HRstar_fstar0_func, HRstar_fstar1_func, HRstar_fstar0_func_g,
#           HRstar_fstar1_func_g, HRstar_fstar0_func_num, HRstar_fstaro0_func_den,
#           HRstar_fstar1_func_num, HRstar_fstar1_func_den
#
#######################################################################################
# Description: Functions to calculate differente summary measures for the HR* 
#
# HRstar_func: average of the HR*
# HRstar_fstar0_func: average of the HR* weighted by f_*^(0)
# HRstar_fstar1_func: average of the HR* weighted by f_*^(1)
# HRstar_fstar0_func_g: gAHR weighted by f_*^(0)
# HRstar_fstar1_func_g: gAHR weighted by f_*^(1)
# HRstar_fstar0_func_num: numerator of AHR weighted by f_*^(0)
# HRstar_fstaro0_func_den: denominator of AHR weighted by f_*^(0)
# HRstar_fstar1_func_num: numerator of AHR weighted by f_*^(1)
# HRstar_fstar1_func_den: denominator of AHR weighted by f_*^(1)
#
# Parameters
# z	             Points of the HR*
# T1pdist        Distribution function of the marginal T1 (pweibull) 
# T2pdist        Distribution function of the marginal T2 (pweibull) 
# T11param       Parameters of marginal endpoint 1 treated group
# T21param       Parameters of marginal endpoint 2 treated group
# distribution1  Distribution treated group
# T10param       Parameters of marginal endpoint 1 control group
# T20param       Parameters of marginal endpoint 2 control group
# distribution0  Distribution control group
#######################################################################################

# MEAN OF HR* WITHOUT WEIGHTING
HRstar_func <- function(z){
  fstar_1 <- (-grad(Sstar, x=z, dist1=T1pdist, dist2=T2pdist, param1=T11param, param2=T21param, dist_biv= distribution1))
  fstar_0 <- (-grad(Sstar, x=z, dist1=T1pdist, dist2=T2pdist, param1=T10param, param2=T20param, dist_biv= distribution0))
  
  Sstar_1 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)
  Sstar_0 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
  
  Lstar_1 <- fstar_1/Sstar_1
  Lstar_0 <- fstar_0/Sstar_0
  
  return(Lstar_1/Lstar_0)
}

# MEAN OF HR* (sAHR) --> f_0
HRstar_fstar0_func <- function(z){
  
  int_fstar_0 <- integrate(fstar0_func,lower=0,upper=1,subdivisions=1000,
                           T1pdist=T1pdist,T2pdist=T2pdist,
                           #T11param=T11param,T21param=T21param,
                           #distribution1=distribution1,
                           T10param=T10param,T20param=T20param,
                           distribution0=distribution0)$value
  
  fstar_1 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1))
  fstar_0 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0))
  
  Sstar_1 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)
  Sstar_0 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
  
  Lstar_1 <- fstar_1/Sstar_1
  Lstar_0 <- fstar_0/Sstar_0
  
  return((Lstar_1/Lstar_0)*fstar_0/int_fstar_0)
  
}
# MEAN OF HR* (sAHR) --> f_1
HRstar_fstar1_func<-function(z){
  
  # int_fstar_1 <- integrate(fstar1_func,lower=0,upper=1,subdivisions=1000)$value
  int_fstar_1 <- integrate(fstar1_func,lower=0,upper=1,subdivisions=1000,
                           T1pdist=T1pdist,T2pdist=T2pdist,
                           T11param=T11param,T21param=T21param,
                           distribution1=distribution1)$value
  
  fstar_1 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1))
  fstar_0 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0))
  
  Sstar_1 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)
  Sstar_0 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
  
  Lstar_1 <- fstar_1/Sstar_1
  Lstar_0 <- fstar_0/Sstar_0
  
  return((Lstar_1/Lstar_0)*fstar_1/int_fstar_1)
}

# MEAN OF HR* (gAHR) --> f_0
HRstar_fstar0_func_g<-function(z,T1pdist,T2pdist,
                               T11param,T21param,distribution1,
                               T10param,T20param,distribution0){
  # int_fstar_0 <- integrate(fstar0_func,lower=0,upper=1,subdivisions=1000)$value
  int_fstar_0 <- integrate(fstar0_func,lower=0,upper=1,subdivisions=1000,
                           T1pdist=T1pdist,T2pdist=T2pdist,
                           #T11param=T11param,T21param=T21param,
                           #distribution1=distribution1,
                           T10param=T10param,T20param=T20param,
                           distribution0=distribution0)$value
  
  fstar_1 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1))
  fstar_0 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0))
  
  Sstar_1 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)
  Sstar_0 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
  
  Lstar_1 <- fstar_1/Sstar_1
  Lstar_0 <- fstar_0/Sstar_0
  
  return(log(Lstar_1/Lstar_0)*fstar_0/int_fstar_0)
  
}

# MEAN OF HR* (gAHR) --> f_1
HRstar_fstar1_func_g<-function(z,T1pdist,T2pdist,
                               T11param,T21param,distribution1,
                               T10param,T20param,distribution0){
  
  int_fstar_1 <- integrate(fstar1_func,lower=0,upper=1,subdivisions=1000,
                           T1pdist=T1pdist,T2pdist=T2pdist,
                           T11param=T11param,T21param=T21param,
                           distribution1=distribution1)$value
  
  fstar_1 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1))
  fstar_0 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0))
  
  Sstar_1 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)
  Sstar_0 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
  
  Lstar_1 <- fstar_1/Sstar_1
  Lstar_0 <- fstar_0/Sstar_0
  
  return(log(Lstar_1/Lstar_0)*fstar_1/int_fstar_1)
}

# MEAN OF HR* (AHR) --> f_0
HRstar_fstar0_func_num<-function(z,T1pdist,T2pdist,
                                 T11param,T21param,distribution1,
                                 T10param,T20param,distribution0){
  # int_fstar_0 <- integrate(fstar0_func,lower=0,upper=1,subdivisions=1000)$value
  int_fstar_0 <- integrate(fstar0_func,lower=0,upper=1,subdivisions=1000,
                           T1pdist=T1pdist,T2pdist=T2pdist,
                           #T11param=T11param,T21param=T21param,
                           #distribution1=distribution1,
                           T10param=T10param,T20param=T20param,
                           distribution0=distribution0)$value
  
  fstar_1 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1))
  fstar_0 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0))
  
  Sstar_1 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)
  Sstar_0 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
  
  Lstar_1 <- fstar_1/Sstar_1
  Lstar_0 <- fstar_0/Sstar_0
  
  return(Lstar_1/(Lstar_0+Lstar_1)*fstar_0/int_fstar_0)
  
}
HRstar_fstar0_func_den <- function(z,T1pdist,T2pdist,
                                   T11param,T21param,distribution1,
                                   T10param,T20param,distribution0){
  # int_fstar_0 <- integrate(fstar0_func,lower=0,upper=1,subdivisions=1000)$value
  int_fstar_0 <- integrate(fstar0_func,lower=0,upper=1,subdivisions=1000,
                           T1pdist=T1pdist,T2pdist=T2pdist,
                           #T11param=T11param,T21param=T21param,
                           #distribution1=distribution1,
                           T10param=T10param,T20param=T20param,
                           distribution0=distribution0)$value
  
  fstar_1 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1))
  fstar_0 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0))
  
  Sstar_1 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)
  Sstar_0 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
  
  Lstar_1 <- fstar_1/Sstar_1
  Lstar_0 <- fstar_0/Sstar_0
  
  return(Lstar_0/(Lstar_0+Lstar_1)*fstar_0/int_fstar_0)
  
}

# MEAN OF HR* (AHR) --> f_1
HRstar_fstar1_func_num<-function(z,T1pdist,T2pdist,
                                 T11param,T21param,distribution1,
                                 T10param,T20param,distribution0){
  
  # int_fstar_1 <- integrate(fstar1_func,lower=0,upper=1,subdivisions=1000)$value
  int_fstar_1 <- integrate(fstar1_func,lower=0,upper=1,subdivisions=1000,
                           T1pdist=T1pdist,T2pdist=T2pdist,
                           T11param=T11param,T21param=T21param,
                           distribution1=distribution1)$value
  
  fstar_1 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1))
  fstar_0 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0))
  
  Sstar_1 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)
  Sstar_0 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
  
  Lstar_1 <- fstar_1/Sstar_1
  Lstar_0 <- fstar_0/Sstar_0
  
  return(Lstar_1/(Lstar_0+Lstar_1)*fstar_1/int_fstar_1)
  
}
HRstar_fstar1_func_den <- function(z,T1pdist,T2pdist,
                                   T11param,T21param,distribution1,
                                   T10param,T20param,distribution0){
  
  #int_fstar_1 <- integrate(fstar1_func,lower=0,upper=1,subdivisions=1000)$value
  int_fstar_1 <- integrate(fstar1_func,lower=0,upper=1,subdivisions=1000,
                           T1pdist=T1pdist,T2pdist=T2pdist,
                           T11param=T11param,T21param=T21param,
                           distribution1=distribution1)$value
  
  fstar_1 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1))
  fstar_0 <- (-grad(Sstar,x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0))
  
  Sstar_1 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)
  Sstar_0 <- Sstar(x=z,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
  
  Lstar_1 <- fstar_1/Sstar_1
  Lstar_0 <- fstar_0/Sstar_0
  
  return(Lstar_0/(Lstar_0+Lstar_1)*fstar_1/int_fstar_1)
  
}



#######################################################################################
# Function: HRstar_function 
#
#######################################################################################
HRstar_function <- function(t,T1pdist,T2pdist,T10param,T20param,T11param,T21param,
                            distribution0,distribution1,copula,
                            b10,b20,b11,b21,theta,p1,p2,beta1,beta2,HR1,HR2,case,tau=seq(0.25,1,0.25)) {
  
  ##-- Densities for both endpoints
  fT10 <- (beta1/b10) * ((t/b10)^(beta1-1)) * (exp(-(t/b10)^beta1))
  fT11 <- (beta1/b11) * ((t/b11)^(beta1-1)) * (exp(-(t/b11)^beta1))
  fT20 <- (beta2/b20) * ((t/b20)^(beta2-1)) * (exp(-(t/b20)^beta2))
  fT21 <- (beta2/b21) * ((t/b21)^(beta2-1)) * (exp(-(t/b21)^beta2))
  
  ##-- Survival for both endpoints
  ST10 <- exp(-(t/b10)^beta1)
  ST11 <- exp(-(t/b11)^beta1)
  ST20 <- exp(-(t/b20)^beta2)
  ST21 <- exp(-(t/b21)^beta2)
  
  ##-- Survival for the composite endpoint
  # No seria posible hacerlo como en el caso I? ---------------------------!!!!!!!!!!!!!!
  if(copula=='Frank'){
    Sstar0 <- (-log(1+(exp(-theta*ST10)-1)*(exp(-theta*ST20)-1)/(exp(-theta)-1))/theta)
    Sstar1 <- (-log(1+(exp(-theta*ST11)-1)*(exp(-theta*ST21)-1)/(exp(-theta)-1))/theta)  
  }else if(copula=='Clayton'){
    Sstar0 <- (ST10^(-theta) + ST20^(-theta) - 1)^{-1/theta}
    Sstar1 <- (ST11^(-theta) + ST21^(-theta) - 1)^{-1/theta}
  }else if(copula=='Gumbel'){
    Sstar0 <- exp(-((-log(ST10))^theta + (-log(ST20))^theta)^(1/theta))
    Sstar1 <- exp(-((-log(ST11))^theta + (-log(ST21))^theta)^(1/theta))      
  }
  
  
  ##-- Density, hazards and hazard ratio for the composite
  # No seria posible hacerlo como en el caso I? ---------------------------!!!!!!!!!!!!!!
  if(copula=='Frank'){
    fstar0 <- (exp(-theta*ST10)*(exp(-theta*ST20)-1)*fT10 + exp(-theta*ST20)*(exp(-theta*ST10)-1)*fT20)/(exp(-theta*Sstar0)*(exp(-theta)-1))
    fstar1 <- (exp(-theta*ST11)*(exp(-theta*ST21)-1)*fT11 + exp(-theta*ST21)*(exp(-theta*ST11)-1)*fT21)/(exp(-theta*Sstar1)*(exp(-theta)-1))
  }else if(copula=='Clayton'){
    fstar0 <- (ST10^(theta+1) * fT10 + ST20^(theta+1) * fT20)/(Sstar0*(ST10^(-theta) + ST20^(-theta) - 1))
    fstar1 <- (ST11^(theta+1) * fT11 + ST20^(theta+1) * fT21)/(Sstar1*(ST11^(-theta) + ST21^(-theta) - 1))
  }else if(copula=='Gumbel'){
    fstar0 <- Sstar0 * log(Sstar0) * ((-log(ST10))^(theta-1) * fT10 * (-ST10)^(-1) + (-log(ST20))^(theta-1)  * fT20 * (-ST20)^(-1))/((-log(ST10))^theta + (-log(ST20))^theta)
    fstar1 <- Sstar1 * log(Sstar1) * ((-log(ST11))^(theta-1) * fT11 * (-ST11)^(-1) + (-log(ST21))^(theta-1)  * fT21 * (-ST21)^(-1))/((-log(ST11))^theta + (-log(ST21))^theta)
  }
  
  ##-- Hazards and hazard ratio for the composite
  Lstar0 <- (fstar0/Sstar0)
  Lstar1 <- (fstar1/Sstar1)
  HRstar <- (Lstar1/Lstar0)
  
  #-- Summary measures for the HR* (see Schempfer 2009)
  HRstar_int <- rowMeans(cbind(HRstar[-1],rev(rev(HRstar)[-1]))) # Mean of HRs for each interval
  fstar0_int <- rowMeans(cbind(fstar0[-1],rev(rev(fstar0)[-1]))) # Mean of fstar0 for each interval
  fstar1_int <- rowMeans(cbind(fstar1[-1],rev(rev(fstar1)[-1]))) # Mean of fstar1 for each interval --> Needed for gAHR
  Lstar0_int <- rowMeans(cbind(Lstar0[-1],rev(rev(Lstar0)[-1]))) # Mean of Lstar0 for each interval
  Lstar1_int <- rowMeans(cbind(Lstar1[-1],rev(rev(Lstar1)[-1]))) # Mean of Lstar1 for each interval
  
  summaryHR <- matrix(NA,nrow=length(tau),ncol=7)                # matrix to store summary HRs
  nHR <- (HR1+HR2)/2                                             # naive HR
  for(tau_i in tau){
    sel <- t<tau_i
    mHR <- mean(HRstar[sel])                                                      # Mean of HR
    sAHR_0 <- sum(HRstar_int[sel]*fstar0_int[sel])/sum(fstar0_int[sel])           # sHR "_0" indicates that f_0 is used instead of f_1 
    
    # gAHR weighted by f0
    gAHR_0 <- exp(sum(log(HRstar_int[sel])*fstar0_int[sel])/sum(fstar0_int[sel])) # gHR "_0" indicates that f_0 is used instead of f_1
    
    # gAHR weighted by f0+f1
    gAHR   <- exp(sum(log(HRstar_int[sel])*(fstar0_int[sel]+fstar1_int[sel]))/sum(fstar0_int[sel]+fstar1_int[sel])) # Alternative gAHR
    
    # AHR weighted by f0
    AHR_0_num <- sum(Lstar1_int[sel]/(Lstar0_int[sel] + Lstar1_int[sel])*fstar0_int[sel])/sum(fstar0_int[sel])
    AHR_0_den <- sum(Lstar0_int[sel]/(Lstar0_int[sel] + Lstar1_int[sel])*fstar0_int[sel])/sum(fstar0_int[sel])
    AHR_0 <- AHR_0_num/AHR_0_den                                                  # AHR "_0" indicates that f_0 is used instead of f_1
    
    # AHR weighted by f0+f1
    AHR_num <- sum(Lstar1_int[sel]/(Lstar0_int[sel] + Lstar1_int[sel])*(fstar0_int[sel]+fstar1_int[sel]))/sum(fstar0_int[sel]+fstar1_int[sel])
    AHR_den <- sum(Lstar0_int[sel]/(Lstar0_int[sel] + Lstar1_int[sel])*(fstar0_int[sel]+fstar1_int[sel]))/sum(fstar0_int[sel]+fstar1_int[sel])
    AHR <- AHR_num/AHR_den                                                  # AHR "_0" indicates that f_0 is used instead of f_1
    
    # Store information
    summaryHR[which(tau==tau_i),] <- c(nHR,mHR,sAHR_0,gAHR_0,AHR_0,gAHR,AHR)
  }
  rownames(summaryHR) <- paste0('tau=',tau)
  colnames(summaryHR) <- c('nHR','mHR','sAHR','gAHR0','AHR0','gAHR','AHR')                                                         # MEAN OF HR* (AHR) --> f_0  
  
  ## -- Another measures
  # RMST
  RMST_0 <- integrate(Sstar, dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0, lower=0,upper=1,subdivisions=1000)$value
  RMST_1 <- integrate(Sstar, dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1, lower=0,upper=1,subdivisions=1000)$value
  
  # Probability of event during follow_up
  pstar_0 <- 1 - Sstar0_func(1,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
  pstar_1 <- 1 - Sstar1_func(1,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)
  
  # Median
  limits <- c(0,10)                                                                         # The first and the last values must be in opposite signs for the function
  Med_0 <- uniroot(Sstar0_func_perc, interval=limits,extendInt="yes", perc=0.5,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)$root  # Find the root (value which equals the function to zero). extendInt="yes": 'interval limits' is extended automatically if necessary
  Med_1 <- uniroot(Sstar1_func_perc, interval=limits,extendInt="yes", perc=0.5,dist1=T1pdist,dist2=T2pdist,param1=T11param,param2=T21param,dist_biv= distribution1)$root  # Find the root (value which equals the function to zero). extendInt="yes": 'interval limits' is extended automatically if necessary
  
  ##########################
  # Probabilities p11,p21 (Although we do not need to calculate the ARE)
  p11 <- 1-exp(-(1/b11)^beta1)
  p21 <- 1-exp(-(1/b21)^beta2)
  
  ##-- Simplified versions of the functions
  selected_times <- t %in% c(0.0001,seq(0.1,1,0.1))  # Only stored 10 points per function
  simplified_HRstar <- HRstar[selected_times]  
  simplified_Sstar0 <- Sstar0[selected_times]  
  simplified_Sstar1 <- Sstar1[selected_times]
  simplified_Lstar0 <- Lstar0[selected_times]  
  simplified_Lstar1 <- Lstar1[selected_times]
  
  
  return(list(HRstar=simplified_HRstar,
              RMST=c(RMST_0=RMST_0,RMST_1=RMST_1),
              pstar=c(pstar_0=pstar_0,pstar_1=pstar_1),
              Med=c(Med_0=Med_0,Med_1=Med_1),
              probs=c(p1=p1,p2=p2,p11=p11,p21=p21),
              scale=c(b10=b10,b20=b20,b11=b11,b21=b21),
              summaryHR=summaryHR,
              Sstar0=simplified_Sstar0,Sstar1=simplified_Sstar1,
              Lstar0=simplified_Lstar0,Lstar1=simplified_Lstar1))
}

#######################################################################################
# Function: COMBINED_HR 
#
#######################################################################################
# Description: It computes the Combined HR(t) (HR_Star) funtion 
# 
# rho0	  Spearman's coefficient between T1 and T2 in control group
# rho1	  Spearman's coefficient between T1 and T2 in treatment group
# beta1	  Shape parameter for a Weibull law for the relevant event
# beta2   Shape parameter for a Weibull law for the additional event 
# HR1     Hazard Ratio for a Weibull law for the relevant event
# HR2     Hazard Ratio for a Weibull law for the additional event
# p1      Proportion of the relevant event expected in group zero
# p2      Proportion of the additional event expected in group zero
# case    Censoring case -- > 1 to 4
# copula  Copula used:
#            Archimedean: "Frank" (default), "Gumbel" or "Clayton"
#            Elliptical: "Normal" or "T"
#            Extreme Value: "Galambos", "HuslerReiss", "Gumbel", "Tawn" or "Tev"
#            Other: "FGM" or "Plackett"
#######################################################################################


COMBINED_HR <- function(t,rho0,rho1=rho0, beta1, beta2, HR1, HR2, p1, p2, case = 1, copula="Frank", rhoType='Spearman',tau=seq(0.25,1,0.25)){ 
  
  
  ###### 1. ELECTION OF THE COPULA
  copula0 <- CopulaSelection(copula,rho0,rhoType)
  theta <- copula0[[2]]   
  which.copula0 <- copula0[[1]]
  which.copula1 <- CopulaSelection(copula,rho1,rhoType)[[1]]  

  ###### 2. ELECTION OF THE MARGINAL DISTRIBUTIONS
  MarginSelec <- MarginalsSelection(beta1,beta2,HR1,HR2,p1,p2,case,theta,copula=copula)
  T1dist   <- MarginSelec[[1]]
  T2dist   <- MarginSelec[[2]]
  T1pdist  <- MarginSelec[[3]]
  T2pdist  <- MarginSelec[[4]]
  T10param <- MarginSelec[[5]]
  T20param <- MarginSelec[[6]]
  T11param <- MarginSelec[[7]]
  T21param <- MarginSelec[[8]]
  
  ###### 3. Scale parameters
  b10 <- T10param[[2]]
  b20 <- T20param[[2]]
  b11 <- T11param[[2]]
  b21 <- T21param[[2]]
  
  ###### 4. Bivariate distribution in control and treatment groups
  distribution0 <- mvdc(copula = which.copula0, margins = c(T1dist, T2dist),paramMargins = list(T10param, T20param))
  distribution1 <- mvdc(copula = which.copula1, margins = c(T1dist, T2dist),paramMargins = list(T11param, T21param))
  
  ###### 5. Compute summary hazard ratios
  resHR <- HRstar_function(t,T1pdist,T2pdist,T10param,T20param,T11param,T21param,
                           distribution0,distribution1,copula,
                           b10,b20,b11,b21,theta,p1,p2,beta1,beta2,HR1,HR2,case)
  
  return(resHR)
}





#######################################################################################
# Function: CopulaSelection 
#
#######################################################################################
# Description: Constructs a copula class object from the family given and the
#              the corresponding dependence parameter grom the given correlation
# 
# copula  Copula given:
#            Archimedean: "Frank" (default), "Gumbel" or "Clayton"
#            Elliptical: "Normal" or "T"
#            Extreme Value: "Galambos", "HuslerReiss", "Gumbel", "Tawn" or "Tev"
#            Other: "FGM" or "Plackett"
# rho	  Spearman's coefficient between the 2 marginal distributions
#######################################################################################

CopulaSelection <- function(copula,rho,rhoType='Spearman'){
  
  if(rhoType=='Spearman'){
    theta <- switch(copula,
                    Frank =       iRho(frankCopula(1),rho),
                    Gumbel =      iRho(gumbelCopula(2),rho),
                    Clayton =     iRho(claytonCopula(1),rho),
                    FGM =         iRho(fgmCopula(1),rho),
                    Normal =      iRho(normalCopula(0.5),rho),
                    'T' =         iRho(normalCopula(0.5),rho), # iRho(tCopula(0.5),rho), --> see Details in ?iRho
                    Galambos =    iRho(galambosCopula(0.5),rho),
                    HuslerReiss = iRho(huslerReissCopula(0.5),rho),
                    Tawn =        iRho(tawnCopula(0.5),rho),
                    Tev =         iRho(tevCopula(0.5),rho),
                    Plackett =    iRho(plackettCopula(0.5),rho))
  }else{
    theta <- switch(copula,
                    Frank =       iTau(frankCopula(1),rho),
                    Gumbel =      iTau(gumbelCopula(2),rho),
                    Clayton =     iTau(claytonCopula(1),rho),
                    FGM =         iTau(fgmCopula(1),rho),
                    Normal =      iTau(normalCopula(0.5),rho),
                    'T' =         iTau(normalCopula(0.5),rho), # iRho(tCopula(0.5),rho), --> see Details in ?iRho
                    Galambos =    iTau(galambosCopula(0.5),rho),
                    HuslerReiss = iTau(huslerReissCopula(0.5),rho),
                    Tawn =        iTau(tawnCopula(0.5),rho),
                    Tev =         iTau(tevCopula(0.5),rho),
                    Plackett =    iTau(plackettCopula(0.5),rho))
  }
  
  which.copula <- switch(copula,
                         Frank =       archmCopula(family = "frank", dim = 2, param = theta),
                         Gumbel =      archmCopula(family = "gumbel", dim = 2, param = theta),
                         Clayton =     archmCopula(family = "clayton", dim = 2, param = theta),
                         FGM =         fgmCopula(dim = 2, param = theta),
                         Normal =      normalCopula(dim = 2, param = theta),
                         'T' =         tCopula(dim = 2, param = theta),
                         Galambos =    galambosCopula(param = theta),
                         HuslerReiss = huslerReissCopula(param = theta),
                         Tawn =        tawnCopula(param = theta),
                         Tev =         tevCopula(param = theta),
                         Plackett =    plackettCopula(param = theta ))
  
  
  return(c(which.copula,theta))
}
#######################################################################################
# Function: MarginalsSelection 
#
#######################################################################################
# Description: Returns the family distribution and parameters of the marginals
#              (ONLY WEIBULL DISTRIBUTIONS SO FAR) 
#
# beta1	  Shape parameter for a Weibull law for the relevant event
# beta2   Shape parameter for a Weibull law for the additional event 
# HR1     Hazard Ratio for a Weibull law for the relevant event
# HR2     Hazard Ratio for a Weibull law for the additional event
# p1      Proportion of the relevant event expected in group zero
# p2      Proportion of the additional event expected in group zero
# case    Censoring case -- > 1 (default) or 3
# theta   Dependence parameter for the bivariate distribution in control group
#######################################################################################

MarginalsSelection <- function(beta1,beta2,HR1,HR2,p1,p2,case,theta,copula='Frank') 
{
  # Scale parameters for group 0 b10,b20
  
  # dcopula <- dfrank
  dcopula <- get(paste0('d',copula))
  ## -- Case 1 -----------------------
  if(case==1) {
    b10 <- 1/((-log(1-p1))^(1/beta1))
    b20 <- 1/((-log(1-p2))^(1/beta2))
  ## -- Case 2 -----------------------
  } else if (case==2) {  
    
    Fb10 <- function(b10,p1){
      integral<-integrate(function(u) {
        sapply(u, function(u) {
          integrate(function(v){dcopula(u,v,theta)}, lower=1-exp((b10*(-log(1-u))^(1/beta1))^beta2*log(1-p2)), upper=1)$value
        })
      }, lower=0 , upper=1-exp(-1/b10^beta1))$value
      return(integral-p1) 
    }
    limits <- c(0.00001,10000)                                         # The first and the last values must be in opposite signs for the function
    b10 <- try(uniroot(Fb10, interval=limits,p1=p1)$root,silent=TRUE)  # Find the root (value which equals the function zero)
    b20 <- 1/(-log(1-p2))^(1/beta2)
    if(class(b10)=='try-error'){
      dcopula <- dFrank
      b10 <- uniroot(Fb10, interval=limits,p1=p1)$root
      dcopula <- get(paste0('d',copula))
      limits <- c(0.8,1.2)*b10 
      b10 <- uniroot(Fb10, interval=limits,p1=p1)$root
    }
  ## -- Case 3 -----------------------
  } else if (case==3) {
    
    b10 <- 1/((-log(1-p1))^(1/beta1))
    
    Fb20 <- function(b20,p2) {
      integral<-integrate(function(v) {
        sapply(v,function(v) { 
          integrate(function(u){dcopula(u,v,theta)},lower=1-exp((b20*(-log(1-v))^(1/beta2))^beta1*log(1-p1)),upper=1)$value
        })
      }, 
      lower=0 , upper=1-exp(-1/b20^beta2))$value
      return(integral-p2)
    }
    limits <- c(0.00001,10000) 
    b20 <- try(uniroot(Fb20, interval=limits,p2=p2)$root,silent=TRUE)
    if(class(b20)=='try-error'){
      # First attemp: to approximate for limits based on frank copula
      dcopula <- dFrank
      b20 <- uniroot(Fb20, interval=limits,p2=p2)$root
      dcopula <- get(paste0('d',copula))
      #limits <- c(0.8,1.2)*b20 # --> Petaba
      limits <- c(0.5,2)*b20
      #print(Fb20(limits[1]*0.5/0.8,p2))
      #print(Fb20(limits[2]*2/1.2,p2))
      b20 <- try(uniroot(Fb20, interval=limits,p2=p2)$root,silent=TRUE)
      
      # Second attemp: search for factible limits
      if(class(b20)=='try-error'){
        FFB20 <- c()
        VALUES <- unique(c(seq(0.1,10,0.1),seq(10,100,1),seq(100,1000,10),seq(1000,10000,100)))
        for (i in 1:length(VALUES)) {aux <- try(Fb20(VALUES[i],p2),silent=TRUE); FFB20[i] <- ifelse(class(aux)=='try-error',NA,aux)}
        POS_LIMITS <- c(neg=NA,pos=NA)
        POS_MIN <- NA
        VAL_MIN <- 1000
        j <- 1
        while(any(is.na(POS_LIMITS)) & j<=length(FFB20)){
          if(is.na(FFB20[j])){
            POS_LIMITS <- c(neg=NA,pos=NA)
          }else if(FFB20[j]<=0){
            POS_LIMITS[1] <- j
            if(abs(FFB20[j])<VAL_MIN){VAL_MIN <- abs(FFB20[j]);POS_MIN <- j}
          }else{
            POS_LIMITS[2] <- j
            if(abs(FFB20[j])<VAL_MIN){VAL_MIN <- abs(FFB20[j]);POS_MIN <- j}
          }
          
          j <- j + 1
        }
        if(!any(is.na(POS_LIMITS))){
          limits <- VALUES[POS_LIMITS]
          b20 <- uniroot(Fb20, interval=limits,p2=p2)$root
        }else{
          b20 <- VALUES[POS_MIN]
          cat('Imposed value to b20:',b20,'Parameters:',beta1,beta2,HR1,HR2,p1,p2,case,theta,copula,'\n',file = "log_marginalselection.txt",append=TRUE)
        }
        
      }
    }
  ## -- Case 4 -----------------------
  } else if (case==4) {
    
    # We need to create x[1] and x[1] (we assign NA's)
    x <- c(NA,NA)
    
    # We need to change the name of variables as (b10=x[1],b20=[2]) to execute 'multiroot'
    # Compute b10
    Fb10 <- function(b10,b20,p1){
      b10-> x[1]
      b20-> x[2]
      integral<-integrate(function(u) {
        sapply(u, function(u) {
          integrate(function(v){dcopula(u,v,theta)}, lower=1-exp(-(x[1]*(-log(1-u))^(1/beta1)/x[2])^beta2), upper=1)$value
        })
      }, lower= 0 , upper=1-exp(-1/x[1]^beta1))$value
      return(integral-p1)
    }
    
    # Compute b20
    Fb20<-function(b10,b20,p2) {
      b10-> x[1]
      b20-> x[2]
      integral<-integrate(function(v) {
        sapply(v,function(v) {
          integrate(function(u){dcopula(u,v,theta)},lower=1-exp(-(x[2]*(-log(1-v))^(1/beta2)/x[1])^beta1), upper=1)$value
        })
      },
      lower= 0, upper=1-exp(-1/x[2]^beta2))$value
      return(integral-p2)
    }
    
    model <- function(x){
      c(Fb10(x[1],x[2],p1), Fb20(x[1],x[2],p2))
    }
    
    (sol <- multiroot(f = model, start = c(1,1)))
    # (sol <- multiroot(f = model, start = c(0.002,0.002), 
    #                   positive = TRUE))
    
    sol<-as.data.frame(sol[1])
    b10<-sol[1,]
    b20<-sol[2,]
    
  }
  
  # Scale parameters for group 1 b11,b21 (Although we do not need to compute the scale parameters for group 1 (b11, b21) to calculate the ARE)
  b11 <- b10/HR1^(1/beta1)
  b21 <- b20/HR2^(1/beta2)
  
  # Probabilities p11,p21 (Although we do not need to calculate the ARE)
  p11 <- 1-exp(-(1/b11)^beta1)
  p21 <- 1-exp(-(1/b21)^beta2)
  
  T1dist<-"weibull"
  T2dist<-"weibull"
  
  T1pdist<-pweibull
  T2pdist<-pweibull
  
  T10param<-list(shape = beta1, scale = b10)
  T20param<-list(shape = beta2, scale = b20)
  
  T11param<-list(shape = beta1, scale = b11)
  T21param<-list(shape = beta2, scale = b21)
  
  return(list(T1dist,T2dist,T1pdist,T2pdist,T10param,T20param,T11param,T21param,p11,p21))	
}


