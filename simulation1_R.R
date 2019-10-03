

rm(list=ls())


library(survival)
library(dplyr)
library(pec)
library(twang)


resultsDIR<- "/netapp2/home/rw751/sentinel/results/"

nstudy <- 10^5
lambda <- 2*10^-5 #1/852            #event rate drives incidence of outcomes, re-scaled in years

time.periods<- 100

laterslice1<- .25
laterslice2<- .50
laterslice3<- .75

a0<-  -3.5
a1<-   1.00
a2<-   3.00
a3<-   3.00

b1<-   1.50
b2<-   4.80
b3<-   1.30


bE<- log(2)



for(sim in 1:1000){

  x1<- runif(nstudy, 0, 1)
  x2<- runif(nstudy, 0, 1)
  x3<- runif(nstudy, 0, 1)


  logodds <-  a0 + a1*x1 + a2*x2 + a3*x3  
  odds    <-  exp(logodds)
  pscore  <-  odds / (odds + 1)
  e       <-  rbinom(nstudy, 1, pscore)

  iptw<- e*(1/pscore) + (1-e)*(1/(1-pscore))
  
  data<- NULL
  for(counter in 0:1){
    
    countfact<- counter
    
    if(countfact==1) e<- 1-e

    LP <- b1*x1 + b2*x2 + b3*x3 + bE*e  

    DRS<- b1*x1 + b2*x2 + b3*x3

    u   <- runif(nstudy, 0, 1)
    time<- -log(u) / (lambda*exp(LP))

    c_driver <- .75 + (6.5*runif(nstudy, 0, 1))
    v = runif(nstudy, 0, 1)
    c_time = -log(v) / (lambda*exp(c_driver))
    c_time = 10000000000  #if no censoring
  
    event <- rep(1, nstudy)
    censored <- c_time < time
    event[censored] <- 0
    time[censored] <- c_time[censored]
    
    id<- rep(1:nstudy)
    
    data<- as.data.frame(rbind(data, cbind(id, x1, x2, x3, e, iptw, time, censored, event, countfact)))

  }

  data1<- data
  #data1<- subset(data, countfact==0)  #if we want full counterfactual pop or not

  data1<- data[order(data1$time),]
  
  # total number of events over time
  n.events<- cumsum(data1$event)
  
  # number of events per time period
  k      <- max(n.events) / time.periods    #N of evnts per interval of event times in full countefactual cohort

  # number of events in each interval after censoring
  p_uncens<- max(n.events) / dim(data1)[1]   #proportion of individual having events
  k<- k * p_uncens

  # population still at risk after 25%, 50%, and 75% of events
  later1<- n.events > (k*laterslice1*time.periods)
  later2<- n.events > (k*laterslice2*time.periods)
  later3<- n.events > (k*laterslice3*time.periods)
  
  
  # basline propensity score (estimated and true)
  bps<- glm(e~x1+x2+x3, data=data1, family=binomial)$fitted.values
  Ee<- pscore
  
  # propensity score in those at risk after 25% have outcome
  ps.later1.fit<- glm(e~x1+x2+x3, data=data1[later1,], family=binomial)
  ps.later1<- predict(ps.later1.fit, newdata=data1, type='response')
  
  # propensity score in those at risk after 50% have outcome
  ps.later2.fit<- glm(e~x1+x2+x3, data=data1[later2,], family=binomial)
  ps.later2<- predict(ps.later2.fit, newdata=data1, type='response')
  
  # propensity score in those at risk after 75% have outcome
  ps.later3.fit<- glm(e~x1+x2+x3, data=data1[later3,], family=binomial)
  ps.later3<- predict(ps.later3.fit, newdata=data1, type='response')
  
  data2<- as.data.frame(cbind(data1, bps, ps.later1, ps.later2, ps.later3, Ee, n.events, later1, later2, later3))
  
  count<- 0
  est1<- NULL
  est2<- NULL
  est3<- NULL
  est4<- NULL
  est5<- NULL
  est6<- NULL
  est7<- NULL
  est8<- NULL
  est9<- NULL
  est10<- NULL
  
  time.intervals<- NULL
  ps.tv.cum<- NULL
  time.cut<- NULL
  time.ps<- NULL
  names.time.ps<- NULL
  data.T<- NULL
  data.T0.T<- NULL
  data.T.cum<- NULL
  data3<- NULL
  
  time.ps<- as.data.frame(matrix(nrow=dim(data2)[1], ncol=time.periods))
  names.time.ps<- NULL
  
  for(slice in 1:time.periods){
  
    count<- count+1
    
    if(count==1){ 
      first <- 0
      last  <- first + k
      mid  <-  first + (k*0.5)	                      #midpoint of the time interval
    }
    if(count>1){  
      first<- (slice-1)*k + 1              #the beginning of the time interval; *recall that k is number of events in each interval
      last <- (first-1) + k                             #end of the time interval
      mid  <- (first-1) + (k*0.5)	                      #midpoint of the time interval
    }
    
    time.intervals[count]<- last
    
    data.T              <- subset(data2, n.events >= first)
    data.T$t            <- data.T$n.events
    test1               <- data.T$n.events > last
    data.T$event[test1] <- 0
    data.T$t[test1]     <- last+1
    
    data.T.mid<- subset(data.T, data.T$n.events >= mid)  #subset of data who are at risk after the midpoint of the interval
    
    data.T0.T              <- data2
    data.T0.T$t            <- data.T0.T$n.events
    test2                  <- data.T0.T$n.events > last
    data.T0.T$event[test2] <- 0
    data.T0.T$t[test2]     <- last+1
    
    
    ##########################################
    ##
    ##            Time Varying PS
    ##
    ##########################################
    
    # when running only actual population 
    ps.model<- glm(e ~ x1 + x2 + x3 +
                     bps + ps.later1 + ps.later2 + ps.later3 + 
                     bps*ps.later1 + bps*ps.later2 + bps*ps.later3 + ps.later1*ps.later2 + ps.later1*ps.later3 + ps.later2*ps.later3
                     , family=binomial, data=data.T.mid)
    
    ps.tv<- predict(ps.model, newdata=data.T,    type="response")
    
    time.ps[,count]<- predict(ps.model, newdata=data2, type='response')
    names.time.ps<- c(names.time.ps, paste("ps", count, sep = "") )
    
    #time.ps.temp<- as.data.frame(time.ps[ ,colSums(is.na(time.ps)) == 0])
    #names(time.ps.temp)<- names.time.ps
    
    #data3<- as.data.frame(cbind(data3, time.ps.temp))
    
    dummy.tv<- data.T0.T$t >= first
    time.cut[dummy.tv]<- last
    
    #if(count==1) start<- rep(0, dim(data.T)[1])
    #if(count>1)  start<- rep(unique(time.cut)[count-1], dim(data.T)[1])
    start<- rep(first, dim(data.T)[1])
    
    data.T<-    as.data.frame(cbind(data.T,    ps.tv, start))
    
    data.T.cum<- as.data.frame(rbind(data.T.cum, data.T))
    
    
    #############################
    ##
    ##       Analysis
    ##
    #############################
  
    # unadjusted
    fit1 <- coxph(Surv(t, event) ~ e, data = data.T)          #at time interval T
    fit2 <- coxph(Surv(t, event) ~ e, data = data.T0.T)       #from baseline to time interval T
    
    # outcome adjusted conditional estimate at time interval T
    fit3 <- coxph(Surv(t, event) ~ e + x1 + x2 + x3, data = data.T)    #at time interval T
    fit4 <- coxph(Surv(t, event) ~ e + x1 + x2 + x3, data = data.T0.T) #from baseline to time interval T  
    
    
    ##################################################
    ##
    ##        Propensity score Analysis
    ##
    ##################################################
    
    # baseline propensity score conditional estimate at time interval T
    fit5 <- coxph(Surv(t, event) ~ e + bps + bps*bps + bps*bps*bps, data = data.T)    #at time interval T 
    fit6 <- coxph(Surv(t, event) ~ e + bps + bps*bps + bps*bps*bps, data = data.T0.T) #from baseline to time interval T 
    
    # baseline propensity score  marginal estimate (weighting)
    fit7 <- coxph(Surv(t, event) ~ e, weights=iptw, data=data.T)
    fit8 <- coxph(Surv(t, event) ~ e, weights=iptw, data=data.T0.T)
    
    # time-dependent propensity score conditional estimate at time interval T
    fit9 <- coxph(Surv(t, event) ~        e + pspline(ps.tv, df=4), data=data.T)

    
    est1[count] <- summary(fit1)$coef[1,1] 
    est2[count] <- summary(fit2)$coef[1,1] 
    est3[count] <- summary(fit3)$coef[1,1] 
    est4[count] <- summary(fit4)$coef[1,1]
    est5[count] <- summary(fit5)$coef[1,1]
    est6[count] <- summary(fit6)$coef[1,1]
    est7[count] <- summary(fit7)$coef[1,1]
    est8[count] <- summary(fit8)$coef[1,1]
    est9[count] <- summary(fit9)$coef[1,1]

  }
  
  est1<- as.data.frame(t(est1))
  est2<- as.data.frame(t(est2))
  est3<- as.data.frame(t(est3))
  est4<- as.data.frame(t(est4))
  est5<- as.data.frame(t(est5))
  est6<- as.data.frame(t(est6))
  est7<- as.data.frame(t(est7))
  est8<- as.data.frame(t(est8))
  est9<- as.data.frame(t(est9))

  
  if(sim==1){write.table(est1,file=paste(resultsDIR, paste("results_test1",".csv",sep=""), sep=""),row.names=FALSE,col.names=T,sep=",",append=TRUE)}
        else{write.table(est1,file=paste(resultsDIR, paste("results_test1",".csv",sep=""), sep=""),row.names=FALSE,col.names=F,sep=",",append=TRUE)}
  
  if(sim==1){write.table(est2,file=paste(resultsDIR, paste("results_test2",".csv",sep=""), sep=""),row.names=FALSE,col.names=T,sep=",",append=TRUE)}
        else{write.table(est2,file=paste(resultsDIR, paste("results_test2",".csv",sep=""), sep=""),row.names=FALSE,col.names=F,sep=",",append=TRUE)}
  
  if(sim==1){write.table(est3,file=paste(resultsDIR, paste("results_test3",".csv",sep=""), sep=""),row.names=FALSE,col.names=T,sep=",",append=TRUE)}
        else{write.table(est3,file=paste(resultsDIR, paste("results_test3",".csv",sep=""), sep=""),row.names=FALSE,col.names=F,sep=",",append=TRUE)}
  
  if(sim==1){write.table(est4,file=paste(resultsDIR, paste("results_test4",".csv",sep=""), sep=""),row.names=FALSE,col.names=T,sep=",",append=TRUE)}
        else{write.table(est4,file=paste(resultsDIR, paste("results_test4",".csv",sep=""), sep=""),row.names=FALSE,col.names=F,sep=",",append=TRUE)}
  
  if(sim==1){write.table(est5,file=paste(resultsDIR, paste("results_test5",".csv",sep=""), sep=""),row.names=FALSE,col.names=T,sep=",",append=TRUE)}
        else{write.table(est5,file=paste(resultsDIR, paste("results_test5",".csv",sep=""), sep=""),row.names=FALSE,col.names=F,sep=",",append=TRUE)}
  
  if(sim==1){write.table(est6,file=paste(resultsDIR, paste("results_test6",".csv",sep=""), sep=""),row.names=FALSE,col.names=T,sep=",",append=TRUE)}
        else{write.table(est6,file=paste(resultsDIR, paste("results_test6",".csv",sep=""), sep=""),row.names=FALSE,col.names=F,sep=",",append=TRUE)}
  
  if(sim==1){write.table(est7,file=paste(resultsDIR, paste("results_test7",".csv",sep=""), sep=""),row.names=FALSE,col.names=T,sep=",",append=TRUE)}
        else{write.table(est7,file=paste(resultsDIR, paste("results_test7",".csv",sep=""), sep=""),row.names=FALSE,col.names=F,sep=",",append=TRUE)}
  
  if(sim==1){write.table(est8,file=paste(resultsDIR, paste("results_test8",".csv",sep=""), sep=""),row.names=FALSE,col.names=T,sep=",",append=TRUE)}
        else{write.table(est8,file=paste(resultsDIR, paste("results_test8",".csv",sep=""), sep=""),row.names=FALSE,col.names=F,sep=",",append=TRUE)}
  
  if(sim==1){write.table(est9,file=paste(resultsDIR, paste("results_test9",".csv",sep=""), sep=""),row.names=FALSE,col.names=T,sep=",",append=TRUE)}
        else{write.table(est9,file=paste(resultsDIR, paste("results_test9",".csv",sep=""), sep=""),row.names=FALSE,col.names=F,sep=",",append=TRUE)}
  

}


data1<-read.csv(paste(resultsDIR1, "results_test1.csv", sep=""))
data2<-read.csv(paste(resultsDIR1, "results_test2.csv", sep=""))
data3<-read.csv(paste(resultsDIR1, "results_test3.csv", sep=""))
data4<-read.csv(paste(resultsDIR1, "results_test4.csv", sep=""))
data5<-read.csv(paste(resultsDIR1, "results_test5.csv", sep=""))
data6<-read.csv(paste(resultsDIR1, "results_test6.csv", sep=""))
data7<-read.csv(paste(resultsDIR1, "results_test7.csv", sep=""))
data8<-read.csv(paste(resultsDIR1, "results_test8.csv", sep=""))
data9<-read.csv(paste(resultsDIR1, "results_test9.csv", sep=""))

est1<- apply(data1, 2, mean)
est2<- apply(data2, 2, mean)
est3<- apply(data3, 2, mean)
est4<- apply(data4, 2, mean)
est5<- apply(data5, 2, mean)
est6<- apply(data6, 2, mean)
est7<- apply(data7, 2, mean)
est8<- apply(data8, 2, mean)
est9<- apply(data9, 2, mean)
est10<- NULL
for(i in 1:length(est9)){
  est10[i]<- mean(est9[1:i])
}









