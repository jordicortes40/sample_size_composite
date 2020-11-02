#######################################################################################
# Function: ARE 
#
#######################################################################################
# Description: It computes ARE (Assymptotic Relative Efficiency)
# Parameters:
# rho0	                Spearman's coefficient between T1 and T2 in control group
# rho1	                Spearman's coefficient between T1 and T2 in treatment group
# beta1	                Shape parameter for a Weibull law for the relevant event
# beta2                 Shape parameter for a Weibull law for the additional event 
# HR1                   Hazard Ratio for a Weibull law for the relevant event
# HR2                   Hazard Ratio for a Weibull law for the additional event
# p1                    Proportion of the relevant event expected in group zero
# p2                    Proportion of the additional event expected in group zero
# case                  Censoring case (1,2,3 or 4)
# copula                Copula used:
#                          Archimedean: "Frank" (default), "Gumbel" or "Clayton"
#                          Elliptical: "Normal"
#                          Other: "FGM"
# rho_thype             Type of correlation (Spearman or Kendall)
#######################################################################################
ARE <- function(rho0,rho1=rho0, beta1, beta2, HR1, HR2, p1, p2, case = case, copula = copula, rhoType='Spearman'){ 
  
  ############################################################
  ##-- 1. SELECTION OF THE COPULA
  copula0 <- CopulaSelection(copula,rho0,rhoType)
  theta <- copula0[[2]]   
  which.copula0 <- copula0[[1]]
  which.copula1 <- CopulaSelection(copula,rho1,rhoType)[[1]]  
  ############################################################
  
  
  
  ############################################################
  ##-- 2. SELECTION OF THE MARGINAL DISTRIBUTIONS
  MarginSelec <- MarginalsSelection(beta1,beta2,HR1,HR2,p1,p2,case,theta,copula=copula)
  T1dist <- MarginSelec[[1]]
  T2dist <- MarginSelec[[2]]
  T1pdist <- MarginSelec[[3]]
  T2pdist <- MarginSelec[[4]]
  T10param <- MarginSelec[[5]]
  T20param <- MarginSelec[[6]]
  T11param <- MarginSelec[[7]]
  T21param <- MarginSelec[[8]]
  ############################################################
  
  ############################################################
  ###### 3. ARE EXPRESSION FOLLOWING: 
  # Gomez G, Lagakos SW. Statistical considerations when using a composite endpoint for comparing treatment groups. 
  # Stat Med. 2013; 32:719-738 (pages 2 and 3).
  
  # Bivariate distribution in control and treatment groups
  distribution0 <- mvdc(copula = which.copula0, margins = c(T1dist, T2dist), paramMargins = list(T10param, T20param))
  distribution1 <- mvdc(copula = which.copula1, margins = c(T1dist, T2dist), paramMargins = list(T11param, T21param))
  
  if(case==1|case==3) {
    
    # Inside the integral in the numerator --> Does not work in a function apart (I do not why)
    inside_integral <- function(t){
      
      Sstar0 <- Sstar(x = t,dist1 = T1pdist,dist2 = T2pdist, param1 = T10param,param2 = T20param, dist_biv = distribution0)
      Sstar1 <- Sstar(x = t,dist1 = T1pdist,dist2 = T2pdist, param1 = T11param,param2 = T21param, dist_biv = distribution1)
      
      fstar0 <- (-grad(Sstar,x = t, dist1 = T1pdist,dist2 = T2pdist, param1 = T10param, param2 = T20param,dist_biv = distribution0))
      fstar1 <- (-grad(Sstar,x = t, dist1 = T1pdist,dist2 = T2pdist, param1 = T11param, param2 = T21param,dist_biv = distribution1))
      
      Lstar0 <- (fstar0/Sstar0)
      Lstar1 <- (fstar1/Sstar1)
      
      HRstar <- (Lstar1/Lstar0)
      
      logHRstar <- log(HRstar)
      
      ##-- Numerical issues
      fstar0[fstar0<0] <- 0
      logHRstar[is.na(logHRstar) | logHRstar==Inf | logHRstar== -Inf] <- 0
      
      return(logHRstar*fstar0)
    }
    
    # Integral in the numerator
    integral <- integrate(inside_integral,lower=0,upper=1,subdivisions=10000,stop.on.error = FALSE) 
    numerator <-(integral$value)^2
    
    # Denominator
    Sstar0_1 <- Sstar(x=1,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0) 
    ST10_1 <- 1-do.call(T1pdist,c(q=1,T10param)) 
    denominator <- ((log(HR1))^2)*(1-Sstar0_1)*(1-ST10_1)
    
    # ARE value
    AREstarT <- (numerator/denominator)
    
    # If the integral is not computed, we assign a missing value
    if(integral$message!="OK") {AREstarT <- NA}
    
  } else
    if(case==2|case==4) {
      
      # Computation of the scale parameter values: b10, b20 
      
      if(case==2) {
        
        # Compute b20
        b20 <- 1/(-log(1-p2))^(1/beta2)
        
        # Compute b10
        Fb10 <- function(b10,p1){
          integral<-integrate(function(u) {
            sapply(u, function(u) {
              integrate(function(v) (	(theta*(1-exp(-theta))*exp(-theta*(u+v)))/ (exp(-theta)+  exp(-theta*(u+v))  - exp(-theta*u)-exp(-theta*v))^2 )    , lower=0, upper=   exp((b10*(-log(u))^(1/beta1))^beta2*log(1-p2)) 	)$value
            })
          }, lower= exp(-1/b10^beta1), upper=1)$value
          return(integral-p1)
        }
        
        limits <- c(0.00001,10000)                        # The first and the last values must be in opposite signs for the function
        b10 <- uniroot(Fb10, interval=limits,p1=p1)$root  # Find the root (value which equals the function to zero)
        
      }
      
      if(case==4) {
        # We need to create x[1] and x[2] to run 'multiroot' function (library: rootSolve) (NA's initially assigned)
        x<-NA
        y<-NA
        x[1]<-x
        x[2]<-y
        
        # We need to change the name of variables as (b10=x[1],b20=[2]) to execute 'multiroot'
        # Compute b10
        Fb10 <- function(b10,b20,p1){
          b10 -> x[1]
          b20 -> x[2]
          integral<-integrate(function(u) {
            sapply(u, function(u) {
              integrate(function(v) (	(theta*(1-exp(-theta))*exp(-theta*(u+v)))
                                      /(exp(-theta)+  exp(-theta*(u+v))  - exp(-theta*u)-exp(-theta*v))^2 )    , lower=0,
                        upper=   exp((x[1]*(-log(u))^(1/beta1))^beta2*(  -1/(x[2]^beta2)       )) 	)$value
            })
          }, lower= exp(-1/x[1]^beta1), upper=1)$value
          return(integral-p1)
        }
        
        # Compute b20
        Fb20 <- function(b10,b20,p2) {
          b10 -> x[1]
          b20 -> x[2]
          integral<-integrate(function(v) {
            sapply(v,function(v) {
              integrate(function(u)((theta*(1-exp(-theta))*exp(-theta*(u+v)))
                                    /(exp(-theta)+exp(-theta*(u+v))-exp(-theta*u)-exp(-theta*v))^2),lower=0,
                        upper=exp(-((((-log(v))^(1/beta2))*x[2])/x[1])^beta1))$value
            })
          },
          lower= exp(-(1/x[2])^beta2), upper=1)$value
          return(integral-p2)
        }
        
        model <- function(x){
          c(Fb10(x[1],x[2],p1), Fb20(x[1],x[2],p2))
        }
        
        (sol <- multiroot(f = model, start = c(1, 1)))
        
        sol<-as.data.frame(sol[1])
        b10<-sol[1,]
        b20<-sol[2,]
        
      }
      
      
      ## Computation of the numerator
      # Note: Only marginal Weibull distributions for fT10, fT20, ST10, ST20.
      
      fT10 <- function(t) (beta1/b10) * ( (t/b10)^(beta1-1) ) * (exp(-(t/b10)^beta1))
      ST10 <- function(t) exp(-(t/b10)^beta1)
      
      fT20 <- function(t) (beta2/b20) * ( (t/b20)^(beta2-1) ) * (exp(-(t/b20)^beta2))
      ST20 <- function(t) exp(-(t/b20)^beta2)
      
      
      # Sstar0 and fstar0 for any copula
      Sstar0 <- function(t) Sstar0 <- Sstar(x=t,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0)
      fstar0 <- function(t) fstar0 <-(-grad(Sstar,x=t,dist1=T1pdist,dist2=T2pdist,param1=T10param,param2=T20param,dist_biv= distribution0))
      
      #########################################
      
      aux21 <- function(t,y) theta*exp(-theta*(ST10(t)+y))*(1-exp(-theta))/(exp(-theta)-exp(-theta*ST10(t))-exp(-theta*y)+exp(-theta*(ST10(t)+y)))^2
      
      aux22 <- function(u) {integrate(aux21,0, ST20(u),t=u,subdivisions=10000)$value} # t=u indicates that we are derivating respect to the other variable in aux21(t,y). That is, respect to y.
      
      lambdaC10 <- function(t) aux22(t)*fT10(t)/Sstar0(t)
      
      lambdaC11 <- function(t) HR1*lambdaC10(t)
      
      aux23 <- function(x,t) theta*exp(-theta*(x+ST20(t)))*(1-exp(-theta))/(exp(-theta)-exp(-theta*x)-exp(-theta*ST20(t))+exp(-theta*(x+ST20(t))))^2
      aux24 <- Vectorize(function(u){integrate(aux23,0,ST10(u),t=u,subdivisions=10000)$value}) #t=u indicates that we are derivating respect to the other variable in aux23(x,t). That is, respect to x.
      
      lambdaC20 <- function(t) aux24(t)*fT20(t)/Sstar0(t)
      
      lambdaC21 <- function(t) HR2*lambdaC20(t)
      
      # EVALUATION OF LambdaC20 BEFORE COMPUTATION (IT MAY FAIL IN CASES 2/4 FOR BETAS = 0.5 BECAUSE IT IS NOT ALWAYS EVALUABLE AT T=0)
      LambdaC20_check <- tryCatch(LambdaC20 <- function(t) integrate(lambdaC20,lower=0,upper=t,subdivisions=10000)$value,error = function(e) e)
      
      # WHENEVER LambdaC20 FAILS, WE INCREASE THE LOWER LIMIT OF INTEGRATION
      lower_LambdaC20 <- 0
      while(inherits(LambdaC20_check, "error")=="TRUE" ){
        lower_LambdaC20 <- lower_LambdaC20 + 0.001
        LambdaC20_check <- tryCatch(LambdaC20 <- function(t) integrate(lambdaC20,lower=0+lower_LambdaC20,upper=t,subdivisions=10000)$value,error = function(e) e)
      }
      
      # Double integral directly to calculate LambdaC20.
      LambdaC20 <- Vectorize(function(t) integrate(lambdaC20,lower=0+lower_LambdaC20,upper=t,subdivisions=10000)$value)
      
      
      # Computation of the hazards for both groups
      Lstar0 <- function(t) lambdaC10(t) + lambdaC20(t)
      Lstar1 <- function(t) lambdaC11(t) + lambdaC21(t)
      
      # Computation of HRstar
      HRstar <- function(t) Lstar1(t)/Lstar0(t)
      logHRstar <- function(t) log(Lstar1(t)/Lstar0(t))
      
      # temp3 <- function(t) logHRstar(t)*fstar0(t)
      # https://stackoverflow.com/questions/43189512/double-integral-in-r
      temp3 <- Vectorize(function(t) logHRstar(t)*fstar0(t))
      
      # EVALUATION OF temp4 BEFORE COMPUTATION (IT MAY FAIL IN CASES 2/4 FOR BETAS = 0.5 BECAUSE IT IS NOT ALWAYS EVALUABLE AT T=0)
      temp4_check <- tryCatch(temp4 <- integrate(temp3,0,1,subdivisions=10000)$value,error = function(e) e)
      
      # WHENEVER temp4 FAILS, WE INCREASE THE LOWER LIMIT OF INTEGRATION
      lower_temp4 <- 0
      while(inherits(temp4_check, "error")=="TRUE" & lower_temp4 < 1){ # add "& lower_temp4 < 1"
        lower_temp4 <- lower_temp4 + 0.001
        temp4_check<-tryCatch(temp4 <- integrate(temp3, lower_temp4, 1, subdivisions=10000)$value,error = function(e) e)
      }
      
      # Double integral directly
      temp4 <- integrate(temp3,0+lower_temp4, 1, subdivisions=10000)$value
      numerator <- (temp4)^2
      
      
      ## Computation of PROBT1UNC
      PROBT1UNC_temp_num <- function(t) exp(-HR2*LambdaC20(t))*Sstar0(t)*lambdaC10(t)
      PROBT1UNC_temp_den <- function(t) exp(-LambdaC20(t))*1/2 + exp(-HR2*LambdaC20(t))*1/2
      
      # PROBT1UNC_temp <- function(t) PROBT1UNC_temp_num(t)/PROBT1UNC_temp_den(t)
      # https://stackoverflow.com/questions/43189512/double-integral-in-r
      PROBT1UNC_temp <- Vectorize(function(t){PROBT1UNC_temp_num(t)/PROBT1UNC_temp_den(t)})
      
      PROBT1UNC_int_check <- tryCatch(integrate(PROBT1UNC_temp,lower=0, upper=1,subdivisions=10000)$value, error = function(e) e)
      
      ############################################
      # WE EVALUATE THE FUNCTION PROBT1UNC_int BECAUSE IT MAY FAIL IN CASES 2/4 FOR BETAS = 0.5 PROBABLY DUE TO THE LOWER LIMITS OF THE INTERGATES lambdaC20 AND temp4.
      # WHEN IT FAILS, WE SEARCH FOR THE MINIMUM EVALUABLE LIMIT OF INTEGARTION FOR lambdaC20 AND temp4; AND WE SET A LOWER LIMITS OF 0.001 FOR THE REST OF INTEGRATES
      # TO ENSURE CONVERGENCE. 
      
      lower_PROBT1UNC_int <- 0
      inc_lower <- 0
      
      while(inherits(PROBT1UNC_int_check, "error")=="TRUE"){
        lower_PROBT1UNC_int<-0.001
        inc_lower<-inc_lower+0.001
        
        aux22 <- function(u) integrate(aux21,0.001, ST20(u),t=u,subdivisions=10000)$value # t=u indicates that we are derivating respect to the other variable in aux21(t,y). That is, respect to y.
        
        lambdaC10 <- function(t) aux22(t)*fT10(t)/Sstar0(t)
        lambdaC11 <- function(t) HR1*lambdaC10(t)
        
        aux23 <- function(x,t) theta*exp(-theta*(x+ST20(t)))*(1-exp(-theta))/(exp(-theta)-exp(-theta*x)-exp(-theta*ST20(t))+exp(-theta*(x+ST20(t))))^2
        aux24 <- function(u) integrate(aux23,0.001,ST10(u),t=u,subdivisions=10000)$value #t=u indicates that we are derivating respect to the other variable in aux23(x,t). That is, respect to x.
        
        lambdaC20 <- function(t) aux24(t)*fT20(t)/Sstar0(t)
        lambdaC21 <- function(t) HR2*lambdaC20(t)
        
        LambdaC20 <- function(t) integrate(lambdaC20,lower=lower_LambdaC20 + inc_lower,upper=t,subdivisions=10000)$value
        
        Lstar0 <- function(t) lambdaC10(t) + lambdaC20(t)
        Lstar1 <- function(t) lambdaC11(t) + lambdaC21(t)
        
        HRstar <- function(t) Lstar1(t)/Lstar0(t)
        logHRstar <- function(t) log(Lstar1(t)/Lstar0(t))
        
        # temp3 <- function(t) logHRstar(t)*fstar0(t)
        # https://stackoverflow.com/questions/43189512/double-integral-in-r
        temp3 <- Vectorize(function(t) logHRstar(t)*fstar0(t))
        
        # Double integral directly
        temp4 <- integrate(temp3,lower_temp4 + inc_lower,1,subdivisions=10000)$value
        numerator <- (temp4)^2
        
        PROBT1UNC_temp_num <- function(t) exp(-HR2*LambdaC20(t))*Sstar0(t)*lambdaC10(t)
        PROBT1UNC_temp_den <- function(t) exp(-LambdaC20(t))*1/2 + exp(-HR2*LambdaC20(t))*1/2
        # PROBT1UNC_temp <- function(t) PROBT1UNC_temp_num(t)/PROBT1UNC_temp_den(t)
        # https://stackoverflow.com/questions/43189512/double-integral-in-r
        PROBT1UNC_temp <- Vectorize(function(t) PROBT1UNC_temp_num(t)/PROBT1UNC_temp_den(t))
        
        PROBT1UNC_int_check<- tryCatch(integrate(PROBT1UNC_temp,lower=0.001, upper=1,subdivisions=10000)$value, error = function(e) e)
      } 
      
      
      PROBT1UNC_int <- integrate(PROBT1UNC_temp,lower=lower_PROBT1UNC_int, upper=1,subdivisions=10000)$value
      
      ############################################
      # ARE VALUE:
      AREstarT <- numerator/((log(HR1)^2) * PROBT1UNC_int * (1-Sstar0(1)))
      AREstarT 
    }
  
  return(AREstarT)
}


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

# Sstar<-function(x,dist1,dist2,param1,param2,dist_biv) { 
#   y <- if(length(x) == 1) c(x,x) else cbind(x,x)
#   return(
#     1
#     - do.call(dist1,c(list(q=x),param1))
#     - do.call(dist2,c(list(q=x),param2))
#     + (pMvdc(y, dist_biv))
#   )
# }




#ARE_check<-ARE(rho0=0.5,beta1=1,beta2=1,HR1=0.975,HR2=0.747,p1=0.035,p2=0.028,case=3,copula="Frank")
#ARE_check<-ARE(rho0=0.5,beta1=1,beta2=1,HR1=0.7,HR2=0.7,p1=0.1,p2=0.1,case=3,copula="Frank")

ARE.ARRAY <- function(x,case,copula,rhoType){
  ARE_value <- tryCatch(ARE(rho0=x[7],beta1=x[1], beta2=x[2], HR1=x[5], HR2=x[6], p1=x[3], p2=x[4], 
                            case = case, copula=copula, rhoType=rhoType)
                        ,error = function(e) e)
  res <- ifelse(inherits(ARE_value, "error")=="TRUE",NA,ARE_value)
  return(res)
}