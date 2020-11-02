print_time <- function(t0,i,total){
  t1 <- difftime(Sys.time(),t0,units='secs')
  t1_num <- as.numeric(t1)
  expect_duration_sec <- t1/i * total
  expect_duration_min <- formatC(expect_duration_sec/60,digits=1,format='f')
  expect_remaining <- t1/i * (total-i)
  expect_finish <- as.character(Sys.time() + expect_remaining)
  cat('Iteration:',i,'/',total,'Expected Duration (minutes):',expect_duration_min,'Finish:',expect_finish,'\n')
}

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

f.log.rank <- function(t,tau,Treatment){
  STATUS <- t<tau
  t[t>=tau] <- tau 
  
  ##-- Long-Rank Test
  LR <- survdiff(Surv(t,STATUS)~Treatment)
  return(1 - pchisq(LR$chisq,1))
}

f.AHR <- function(t,tau,Treatment){
  STATUS <- t<tau
  t[t>=tau] <- tau 
  
  ##-- AHR Test
  fit <- ahrWKM(tau,formula=Surv(Y, D) ~ Z, data=data.frame(Y=t, D=STATUS, Z=Treatment))
  # fit <- avgHR(tau, data.frame(Y=t, D=STATUS, Z=Treatment), formula=Surv(Y, D) ~ Z)
  return(2*(1-pnorm(abs(fit$Z.hr))))
}

f.Cox <- function(t,tau,Treatment){
  STATUS <- t<tau
  t[t>=tau] <- tau 
  
  ##-- Cox Model
  fit <- coxph(formula=Surv(t,STATUS) ~ Treatment)
  return(fit$coef)
}

sum_events <- function(x,tau) sum(x<tau)


fc <- function(x) formatC(x,digits=2,format='f')
