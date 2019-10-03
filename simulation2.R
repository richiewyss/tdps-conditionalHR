
#### start timer ####
t_start <- Sys.time()

rm(list=ls())

setwd("/netapp2/home/rw751/sentinel/newest/")

resultsDIR<- "/netapp2/home/rw751/sentinel/results/plasmode/"


library(survival)
library(plyr)
library(dplyr)
library(pec)
library(twang)


######################################################################################################
##
##    Specify survival and censoring models (Optional: Model Survival and Censoring in real data)
##
######################################################################################################

###############################   Cox model
keep <- c("TIME", "EVENT", "EXPOSURE", covars_bross)
smod1 <- coxph(Surv(TIME,  EVENT) ~ ., data = subset(concurrent, select = keep), x=TRUE)

## extract survival and censoring objects from fitted models
Xmat<- as.matrix(smod1$x)
n<- nrow(Xmat)
beta1<- coef(smod1)
fit <- survfit(smod1)
s0  <- fit$surv      # survival curve for average patient
ts  <- fit$time
nts <- length(ts)


####################################################################################################
##
##      This is where global parameters and function begin
##
####################################################################################################

#parameters
effectOR<- 2.0
MMOut=1
MMExp=1
exposedPrev<- NULL
eventRate<- 0.05
size<- 40000  #dim(concurrent)[1]
cut.points<- 20
strength.conf<- 1
censoring<- 0
nsim<- 1000
iter<- 1
scenario<- 8



######################################################################################################
##
##              Simulating treatment assignment (this step is optional)
##
######################################################################################################
  
keep <- c("EXPOSURE", covars_bross)
data.temp<- subset(concurrent, select = keep)
coef.ps<- glm(EXPOSURE~., family=binomial, data=data.temp)$coef[-1]
  
Xa<- as.vector(Xmat[,-1] %*% coef.ps)
  
exposedPrev<- mean(concurrent$EXPOSURE)
fn <- function(d) mean(plogis(d + Xa)) - exposedPrev
exposuredelta <- uniroot(fn, lower = -10, upper = 10)$root  ##Richie comment: set alpha0 to get exposure prevalence that you specify.
Ee <- plogis(exposuredelta + Xa)
EXPOSURE<- rbinom(dim(Xmat)[1], 1, Ee)
Xmat[,1]<- EXPOSURE
  
######################################################################################################
##
##                Setting Event Rate & correlation between PS and DRS
##
######################################################################################################

## find event rate in base cohort (if everyone was followed to end of study)
Xb <- as.vector(Xmat %*% beta1)   #linear predictor
s0end <- min(s0)
if(is.null(eventRate)) eventRate <- 1-mean(s0end^exp(Xb - mean(Xb)) ) 

## set correlation between PS and linear predictor of hazard
bnew1 <- replace(MMOut*beta1, names(beta1) == 'EXPOSURE', log(effectOR))
Xbnew1 <- as.vector(Xmat %*% bnew1)

## estimate PS to evaluate corrlation between PS and Linear Predictor for Outcome (only needed if step 3a is commented out)
#fit.ps<- glm(EXPOSURE ~ ., family=binomial, data=as.data.frame(Xmat))
#Ee<- fit.ps$fitted.values
#coef.ps<- fit.ps$coef[-1]

###################################################################################################
##
##            Setting correlation between PS and DRS for bross covs
##
###################################################################################################

# correlation as seen in data (cor=-0.38)
if(strength.conf==1){
  bnew1.update<- bnew1
  Xbnew1<- Xbnew1
  cor(Ee, Xbnew1, method="spearman")
}

# moderate correlation between propensity and risk (cor=-0.19)
if(strength.conf==2){
  bnew1.update<- bnew1
  beta.new<- bnew1[-1]
  coef1.ps<- coef.ps
  test1<- (beta.new < 0 & coef1.ps < 0) | (beta.new > 0 & coef1.ps > 0)
  beta.new[test1]<- beta.new[test1]*-1
  bnew1.update[-1]<- beta.new
  bnew1.update[15:40]<- -1*bnew1.update[15:40]
  Xbnew1 <- as.vector(Xmat %*% bnew1.update)
  cor(Ee, Xbnew1, method="spearman")
}

# weak correlation between propensity and risk (cor=-0.08)
if(strength.conf==3){
  bnew1.update<- bnew1
  beta.new<- bnew1[-1]
  coef1.ps<- coef.ps
  test1<- (beta.new < 0 & coef1.ps < 0) | (beta.new > 0 & coef1.ps > 0)
  beta.new[test1]<- beta.new[test1]*-1
  bnew1.update[-1]<- beta.new
  bnew1.update[15:47]<- -1*bnew1[15:47]
  Xbnew1 <- as.vector(Xmat %*% bnew1.update)
  cor(Ee, Xbnew1, method="spearman")
}

# very low correlation between propensity and risk (cor=-0.02)
if(strength.conf==4){
  bnew1.update<- bnew1
  beta.new<- bnew1[-1]
  coef1.ps<- coef.ps
  test1<- (beta.new < 0 & coef1.ps < 0) | (beta.new > 0 & coef1.ps > 0)
  beta.new[test1]<- beta.new[test1]*-1
  bnew1.update[-1]<- beta.new
  bnew1.update[14:39]<- -1*bnew1[14:39]
  Xbnew1 <- as.vector(Xmat %*% bnew1.update)
  cor(Ee, Xbnew1, method="spearman")
}

  
####################################################################################################
##
##        Simulate Survival curves for each individual at each time point
##
####################################################################################################
  
## find delta value needed to get approximate desired event rate under new parameters
sXend <- s0end^(exp(Xbnew1 - mean(Xb)))               
fn <- function(d) mean(sXend^d) - (1 - eventRate)
delta <- uniroot(fn, lower = 0, upper = 20)$root

# setup matrix of individual survival curve under new parameters
Sx <- matrix(unlist(lapply(s0, function(s) s^(delta*exp(Xbnew1 - mean(Xb))))), nrow = n)
  
######################################################################################################
##
##        Setting Exposure Prevalence
##
######################################################################################################

## Setting size of exposed and unexposed to get desired exposure prevalence
if(is.null(exposedPrev)) exposedPrev <- mean(Xmat[,1])

##created treated and control groups to sample separately from treated and control to match exposed prevalence specified above.
data.treated<- subset(Xmat, Xmat[,1]==1) 
data.control<- subset(Xmat, Xmat[,1]==0) 
n.treated<- dim(data.treated)[1]
n.control<- dim(data.control)[1]
size.treated<- round(size*exposedPrev)
size.control<- round(size*(1-exposedPrev))


###############################################################################################################
##
##    Start Simulation and Simulate Event and Censoring Times based on survival and censoring curves 
##
###############################################################################################################

#### sample and simulate
for(iter in 1:nsim) {
  
  time0<- Sys.time()
  set.seed(iter)
  
  # Sampling index
  #index <- sample(n, size, replace = TRUE)                                    # when exposure prevalence is same as empirical data
  index.treated <- sample(n, size.treated, replace = TRUE, prob=Xmat[,1])      # when setting exposure prevalence
  index.control <- sample(n, size.control, replace = TRUE, prob=(1-Xmat[,1]))
  index<- c(index.treated, index.control)

  ####################################################
  ## Time to event and survival probability 
 
  # event time
  u <- runif(size, 0, 1)
  w1.event <- apply(Sx[index,] < u, 1, function(x) which(x)[1])  # the first time survival drops below u
  stime <- ts[w1.event]
  w2.event <- Sx[index, nts] > u                                 # for any individuals with survival that never drops below u,
  stime[w2.event] <- max(ts) + 1                                 # replace with arbitrary time beyond last observed event/censoring time

  ##########################################################
  ## censoring time
  
  # no censoring directly
  if(censoring==0){
    ctime<- rep(max(ts), size)
  }

  # with censoring (sampling from actual time distribution)
  if(censoring==1){
    # random uniform censoring
    #ctime<-  sample(1:length(ts), size, replace=TRUE)
    
    # empirical distribution censoring times
    ctime<- sample(concurrent$TIME, size, replace=FALSE)
    index.ctime<- ctime == 0
    ctime[index.ctime]<- 1
  }
  
  ###########################################################
  ## put it together
  tnew <- pmin(stime, ctime)
  ynew <- ifelse(stime == tnew, 1, 0)
  xcovs<- Xmat[index,]
  #xcovs<- concurrent[,c('EXPOSURE',all_covars)]
  
  finaldata<- as.data.frame(cbind(tnew, ynew, xcovs))
  names(finaldata)<- c("TIME", "EVENT", names(finaldata)[-c(1:2)])

  ## percent event and censored before end of follow-up
  p.event.before <- 1-mean(w2.event)                                    #percent having event before censoring 
  p.event.after <- mean(finaldata$EVENT)                                #percent having event after censoring
  p.cens.after  <- mean(ifelse(ctime < stime & ctime<max(ts), 1, 0))    #percent censored after events
  p.event.cens<- 1-p.event.after/p.event.before                         #percentage of events that were censored
  
  data.cens<- as.data.frame(cbind(p.event.before, p.event.after, p.cens.after, p.event.cens))
    
  ###################################################################################################
  ##
  ##    Create discrete time indicators AND final data set AND counterfactual data set
  ##
  ###################################################################################################
  
  # cutting time based on percentage of outcomes
  time.quant<- unique(quantile(finaldata[finaldata$EVENT==1,]$TIME, (1:cut.points)/cut.points), label=(1:cut.points), include.lowest=TRUE)
  time.quant[which(time.quant==max(time.quant))]<- max(finaldata$TIME)
  time.quant<- round(time.quant)
  time.cat<- as.numeric(cut(finaldata$TIME, time.quant))
  time.cat[is.na(time.cat)]<- 0
  time.cat<- time.cat+1

  # final dataset
  finaldata2<- as.data.frame(cbind(tnew, ynew, time.cat, xcovs))
  names(finaldata2)<- c("TIME", "EVENT", "TIME.Cat", names(finaldata2)[-c(1:3)])
  row.names(finaldata2)<- seq(1,dim(finaldata)[1],1)
  finaldata2$id<- row.names(finaldata2)
  finaldata2$TIME.0<- rep(0, dim(finaldata2)[1])
  
  ####################################################################################################
  ##
  ##      Estimate Propensity score
  ##
  ####################################################################################################
    
  # basline propensity score (estimated and true)
  keep <- c("EXPOSURE", covars_bross)
  data.temp<- subset(finaldata2, select = keep)
  ps<- glm(EXPOSURE~., data=data.temp, family=binomial)$fitted.values

  finaldata2<- as.data.frame(cbind(finaldata2, ps))

  ####################################################################################################
  ##
  ##   Analysis
  ##
  ####################################################################################################
  
  # unadjusted
  fit1 <- coxph(Surv(TIME.0, TIME, EVENT) ~ EXPOSURE, data = finaldata2)   
  est1<-    summary(fit1)$coef[1,1]

  # outcome regression
  fit2 <- coxph(Surv(TIME.0, TIME, EVENT) ~ ., data = finaldata2[,-match(c("TIME.Cat", "id"), names(finaldata2))])  
  est2<-    summary(fit2)$coef[1,1]

  # 1-to-1 matched analysis
  ratio=1
  data.temp<- finaldata2
  source("/netapp2/home/rw751/sentinel/analysis_matching2.R", local=TRUE)
  m.data1<- matched.data
  m.data1.cond.total<- matched.data.cond
  
  fit3 <- coxph(Surv(TIME, EVENT) ~ EXPOSURE + strata(m.set), data = m.data1.cond.total)  
  est3<-    summary(fit3)$coef[1,1]
  
  fit4<- coxph(Surv(TIME, EVENT) ~ EXPOSURE, data = m.data1) 
  est4<-    summary(fit4)$coef[1,1]

  ##################################################
  ##
  ##          PS regression
  ##
  ##################################################
    
  fit5<- coxph(Surv(TIME, EVENT) ~ EXPOSURE + pspline(ps, df=4), data = finaldata2)

  est5<-    summary(fit5)$coef[1,1]

  estimates<- as.data.frame(cbind(est1, est2, est3, est4, est5))
  names(estimates)<- c('crude', 'reg', 'match.cond', 'match', 'psreg1')
    
    
  ####################################################
  ## 
  ##            All discrete times
  ##
  ####################################################
  
  d.time<- length(unique(time.cat))
  discrete.times<- as.data.frame(cbind(d.time))
    
    
  #################################################
  ##
  ##      Analysis 2: Estimating Time Dependent PS
  ##
  #################################################
    
  source("/netapp2/home/rw751/sentinel/time_PS.R", local=TRUE)

  ################################################################################
  ##
  ##              Adjustment for Time-Dependent PS in outcome model
  ##
  ################################################################################
    
  names(time.ps)<- names.time.ps
  
  ## creating data for time-varying variables
  finaldata3<- as.data.frame(cbind(finaldata2, time.ps))
    
  ## Splitting data for time-dependent PS
  time.quant.new<- time.quant[seq(0, max(finaldata2$TIME.Cat), by=2)]
  finaldata3.split<- survSplit(data=finaldata3, cut=time.quant.new, end="TIME", start="TIME0", event="EVENT", id="idnew", episode="num")
  
  psnew<- rep(NA, dim(finaldata3.split)[1])
  
  for(ii in 1:length(unique(finaldata3.split$TIME0))){
    
    ps.temp<- paste('ps', ii, sep="")
    index.temp<- finaldata3.split$TIME0 == unique(finaldata3.split$TIME0)[ii]
    psnew[which(index.temp)]<- finaldata3.split[which(index.temp), ps.temp]
    
  }
    
  finaldata3.split$psnew<- psnew
  
  fit6<- coxph(Surv(TIME0, TIME, EVENT) ~ EXPOSURE + pspline(psnew, df=4), data = finaldata3.split)
  est6<- summary(fit10)$coef[1,1]

} 
  
  
