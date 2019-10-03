*options nomacrogen nosymbolgen;

*ods html close; 
*ods html;

%global hr n tx_intercept event_rate k laterslice1 laterslice2;
   %let hr = 2;                     *true conditional HR;
   %let n  = (10**5);               *size of cohort;
   %let tx_intercept= -3.5;         *drives the percent in treatment arm;
   %let event_rate   = 2*(10**-5);  *drives incidence of outcomes, re-scaled in years;
   %let k = (&n ) / 100;            *N of evnts per interval of time;
   %let laterslice1 = 25;
   %let laterslice2 = 50;
   %let laterslice3 = 75;

%macro step1;

data cohort (keep=tx time covar1-covar3 pscore iptw link event drs);
   length tx event link 3;
      link=1; total_events = 0;
   call streaminit(0);      *seed for random number generator;
   n = &N;                  *n of people in each simulated cohort;
   HR = &HR;                *true conditional hazard ratio;
   logHR = log(HR);
   lamda = &event_rate;     *baseline hazard (if covars are all zero);

   array covars {3} covar1-covar3;
   array b_coef {3} b1-b3;
   array c_coef {3} c1-c3;

   b1= 1.00; b2= 3.00; b3= 3.00;
   c1= 1.50; c2= 4.80; c3= 1.30;


do id = 1 to n;
   logodds = &tx_intercept;
   do j = 1 to 3;
      covars{j} = rand("uniform");
      logodds = logodds + (b_coef{j} * covars{j}); *logit of prob of treatmnt;
   end;
   odds    = exp(logodds);
   pscore  = odds / (odds + 1);
   tx      = rand("Bernoulli",pscore); *the treatment is assigned here;

   if tx = 1 then iptw = 1/pscore;  else iptw = 1 / (1-pscore);

   drs = 0;
   do j = 1 to 3;
      drs = drs + (c_coef{j} * covars{j}); *disease risk score;
   end;
   LP   = logHR*tx + drs;      *linear predictor for mechanism generating outcome event times;
   u    = rand("uniform");
   time = -log(u) / (lamda*exp(LP));    *exponential (from Weibull,eta=1, Austin Stat Med 2013,p2840);
   logtime=log(time);

   c_driver = .75 + (6.5*rand("uniform"));
   v = rand("uniform");
   c_time = -log(v) / (lamda*exp(c_driver));
   if c_time < time then event=0;
   else event =1;
   time = min(time,c_time);
   total_events = total_events + event;
 if id = n then call symput ('total_events', total_events);
   label time = 'time to event'
         pscore= 'true propensity score'
         drs   = 'disease risk score ignoring tx';
output cohort;
end; run;

proc sort data=cohort; by time; run;

data cohort;
 set cohort; by link;
 length n_event 6;
 retain n_event 0;
 p_uncens = &total_events / &n;
 k = &k*p_uncens;
 n_event = n_event + event;
 if n_event > (k * &laterslice1) then tx_later1 = tx;
    else tx_later1=.;
 if n_event > (k * &laterslice2) then tx_later2 = tx;
    else tx_later2=.;
 if n_event > (k * &laterslice3) then tx_later3 = tx;
    else tx_later3=.;
run;

proc logistic data = cohort descending;
     model tx = covar1-covar3;
     output out=cohort(drop=_level_) pred = bps;
     title 'bps'; run;
proc logistic data = cohort descending;
     model tx_later1 = covar1-covar3;
     output out=cohort(drop= _level_) pred = pslater1;
     title 'pslater1'; run;
proc logistic data = cohort descending;
     model tx_later2 = covar1-covar3;
     output out=cohort(drop= _level_) pred = pslater2;
     title 'pslater2'; run;
proc logistic data = cohort descending;
     model tx_later3 = covar1-covar3;
     output out=cohort(drop= _level_) pred = pslater3;
     title 'pslater3'; run;
data cohort; set cohort;
     bps_sq = bps**2; bps_cubed = bps**3;
     pslater1_sq = pslater1**2; pslater1_cubed = pslater1**3;
     pslater2_sq = pslater2**2; pslater2_cubed = pslater2**3;
     pslater3_sq = pslater3**2; pslater3_cubed = pslater3**3; run;
%mend step1;

%macro step2(slice,K);

data time_t;
   set cohort;
   length event tx_begin tx_mid tx_final 3;
   t = n_event;
   first = (&slice-1)*k + 1;
   last  = (first-1) + k;
   mid   = (first-1) + (k*0.5);
   if n_event < first then delete;
   if n_event > last then do;
      event=0;
      t = last+1;
   end;
   tx_mid   = tx;
   tx_begin = tx;
   tx_final = tx;
   if n_event <= mid then tx_mid = .;
   if n_event <= last then tx_final = .;
run;

ods output parameterestimates=parms1;
proc phreg data=time_t;
   model t*event(0) = tx covar1-covar3;
   title "Covariate adjusted HRc at time t";
   data covadj_piece; set parms1;
        if parameter='tx'; link = 1;
        covadj_lnhr = estimate;
        keep covadj_lnhr link; run;

ods output parameterestimates=parms2;
proc phreg data=time_t;
   weight iptw;
   model t*event(0) = tx;
   title "IPTW estimate of HRm at time t (no robust variance: COVS dropped to reduce comp burden)"; run;
   data iptwadj_piece; set parms2;
         if parameter='tx'; link = 1;
         iptwadj_lnhr = estimate;
         keep iptwadj_lnhr link; run;

ods output parameterestimates=parms3;
proc phreg data=time_t;
   model t*event(0) = tx bps bps_sq bps_cubed;
   title "HRc estimated with adjustment for a polynomial function of the usual baseline PS";  run;
   data bpsadj_piece; set parms3;
          if parameter='tx'; link = 1;
          bpsadj_lnhr = estimate;
          keep bpsadj_lnhr link; run;

/*Make the time-varying PS*/
ods trace off;
ods output parameterestimates=ps_parms;
proc logistic data = time_t descending;
     model tx_mid =   covar1-covar3
                      bps pslater1 pslater2 pslater3
                      bps*pslater1 bps*pslater2 bps*pslater3 pslater1*pslater2 pslater1*pslater3 pslater2*pslater3;
     output out=time_t(drop=_level_)   pred = ps;
     title 'make the midpoint PS';    run;
data ps_parms; set ps_parms; keep variable waldchisq; run;
proc transpose data=ps_parms out=ps_parms;
  id variable;
  run;
data ps_parms; set ps_parms;
  drop _name_ _label_ intercept;
  link=1; run;
run;

data time_t; set time_t;
     ps_sq = ps**2; ps_cubed = ps**3; 
run;

/*Adjust the HRc estimate for the time-varying PS */
ods output parameterestimates=parms4;
proc phreg data=time_t;
   model t*event(0) = tx ps ps_sq ps_cubed;
   title "HRc estimated with adjustment for a polynomial function of the time-varying PS";
   title2 "Slice = &slice"; run;
   data psadj_piece; set parms4;
        if parameter='tx'; link = 1;
        psadj_lnhr = estimate;
        keep psadj_lnhr link; 
run;

data r;
 merge covadj_piece psadj_piece bpsadj_piece iptwadj_piece ps_parms; by link;
 true_HR = &HR;
 replication = &replication;
 pctl = &slice;
 n = &n;
 k=&k;
 laterslice1 = &laterslice1;
 laterslice2 = &laterslice2;
 laterslice3 = &laterslice3;
  dif =  covadj_lnhr - psadj_lnhr;
    absdif=  abs(dif);
  bias = log(true_hr) - psadj_lnhr;
    absbias = abs(bias);
  label bias = "truth minus psadj_lnhr"
        dif  = "covadj minus psadj";
run;

libname x 'i:stat\sim_res\bigsim\';
proc datasets;
 append base = simres40 new=r force nowarn; run;
 run;
/*
proc datasets;
 append base = x.tvps_3_12_18_smalln new=r force nowarn; run;
 run;
*/
%mend step2;

%macro repeat;
  %step1;
  %do slice =  1 %to 95 %by 1;
     %step2(&slice,&k);
  %end;
%mend repeat;

/*proc printto log='temp.log' new; run;*/

%macro again;
 %do replication = 1  %to 2 %by 1;
  %repeat;
 %end;
%mend again;

%again;


proc sort data=simres40; by pctl replication; run;

*ods html close; 
*ods html;

proc corr;
 var dif absdif bias absbias;
 with pctl;
title 'correlation of percentile with bias and with difference between covadj and psadj';
run;

proc means n nmiss mean stderr p50 p5 p95 min max;
 var dif absdif bias absbias;
title 'dif and bias';
run;


proc tabulate data= simres40;
class pctl;
var replication psadj_lnhr bpsadj_lnhr covadj_lnhr iptwadj_lnhr;
table (pctl), (replication*n) (covadj_lnhr psadj_lnhr bpsadj_lnhr iptwadj_lnhr)*(mean*f= 6.4 stderr*f=5.4)
      /box='Time period (t)';
label psadj_lnhr = "(2) log HRc at t, adjusted by time-varying PS"
      bpsadj_lnhr ="(3) log HRc at t, adjusted by baseline PS"
      covadj_lnhr ="(1) log HRc at t, adjusted by individl covars"
      iptwadj_lnhr = "(4) log HRm at t, adjusted by IPTW"
      pctl = '% dead by end of t'
      replication = 'Reps'
;
    keylabel p50="Median" ;
title "Table 1. Four HR estimates: 3 of conditnl logHR (true HRc=&HR) and 1 of marginal logHR";
title2 "CENSORED analyses of treatment effect on mortality, N=1M, 50% censored randomly";
title3 'Cox regr adj by (1) individl covars (2) time-varying PS (3) baseline PS (4) IPTW';
title4 'Each row has logHR estimates for a time period including the deaths of 1% of cohort';
footnote1 '12 baseline covars are indepndnt and uniformly distributed, and the corr of baseline PS with DRS: .84';
footnote2 "All 12 baseline  covars increase risk and also increase prob. of tx, initially 25% are treated";
footnote3 "12 covar coefficients for predictor of tx: .25, .75, .75, .25, .75, .75, .20, .70, .60, .30, .80, .90";
footnote4 "coeffs for predtr of expontl time-to-event:.35, 1.2, .30, .35, 1.2, .30, .35, 1.2, .30, .35, 1.2, .30";
run;

proc summary nway data=simres40;
     class pctl;
     var covadj_lnhr psadj_lnhr bpsadj_lnhr iptwadj_lnhr;
     output out=means mean=covadj_lnhr psadj_lnhr bpsadj_lnhr iptwadj_lnhr;
run;

proc print data=means; run;

data todate;
  set means;
  retain covsum pssum bpssum iptwsum wtsum 0;
     wt = sqrt(_freq_);
     wtsum = wtsum + wt;
     covsum = covsum   + (covadj_lnhr  *wt);
     pssum = pssum     + (psadj_lnhr   *wt);
     bpssum = bpssum   + (bpsadj_lnhr  *wt);
     iptwsum = iptwsum + (iptwadj_lnhr *wt);
  cov_lnhr_td  = covsum  /wtsum;     cov_hr_td = exp(cov_lnhr_td);
  ps_lnhr_td   = pssum   /wtsum;      ps_hr_td = exp(ps_lnhr_td);
  bps_lnhr_td  = bpssum  /wtsum;     bps_hr_td = exp(bps_lnhr_td);
  iptw_lnhr_td = iptwsum /wtsum;    iptw_hr_td = exp(iptw_lnhr_td);
     row = _n_; 
run;

proc print data=todate; run;

proc print data=todate split='*';
 id row;
 var pctl cov_hr_td ps_hr_td bps_hr_td iptw_hr_td;
 format cov_hr_td ps_hr_td bps_hr_td iptw_hr_td 6.2;
 label
      row = 'Time (t)'
      pctl = '% of cohort*depleted*by events*by end of t'
      cov_hr_td = '(1) HRc from*start thru t*adjustd by*covariates'
      ps_hr_td  = '(2) HRc from*start thru t*adjustd by*time-varying PS'
      bps_hr_td  = '(3) HRc from*start thru t*adjustd by*baseline PS'
      iptw_hr_td  = '(4) HRm from*start thru t*estimatd by*IPTW' ;
title "Table 2. Four HR estimates: 3 of conditional HR (HRc=&HR) and 1 of marg HR";
title2 "CENSORED analyses of treatment effect on mortality, N=1M, 50% censored randomly";
title3 'Cox regr adj by (1) individ covars (2) time-varying PS (3) baseline PS (4) IPTW,';
title4 "for follow-up from the start thru 20 possible endpoints";
footnote1 '12 baseline covars are independnt and uniformly distr, and the corr of baseline PS with DRS: .845';
footnote2 "All 12 baseline  covars increase risk and also increase prob. of tx, initially 25% are treated";
footnote3 "12 covar coefficients for predictor of tx: .25, .75, .75, .25, .75, .75, .20, .70, .60, .30, .80, .90";
footnote4 "coefs for predictr of expontl time-to-evnt:.35, 1.2, .30, .35, 1.2, .30, .35, 1.2, .30, .35, 1.2, .30";
run;


