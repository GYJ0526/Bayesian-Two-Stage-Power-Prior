libname NOAA 'C:\Users\yeongjin.gwon\OneDrive - University of Nebraska Medical Center\01 Academics\01 Research\04 Spatial Modeling\02 NOAA\1. Manuscript\3. Northern Rockies and Plains\Output';

/*******************************************/
/* Overall analysis for general population */
/*******************************************/
proc import out=NOAA.test
  datafile= "C:\Users\yeongjin.gwon\OneDrive - University of Nebraska Medical Center\01 Academics\01 Research\04 Spatial Modeling\02 NOAA\1. Manuscript\3. Northern Rockies and Plains\Output\NorthernRockies.csv"
  dbms=csv replace;
  getnames=yes;
  datarow=2;
run;

/* Descriptive statistics */
proc means data=NOAA.test min p25 p50 p75 max mean maxdec=3;
  class State Year;
  var TD CD RD;
run;

proc means data=NOAA.test min p25 p50 p75 max mean maxdec=3;
  class Year;
  var Drought SPEIB06 SPEIB12 EDDIB06 EDDIB12;
run;

proc means data=NOAA.test min p25 p50 p75 max mean maxdec=3;
  class Year;
  var Dmod Dsev SPEIM06 SPEIM12 SPEIS06 SPEIS12 EDDIM06 EDDIM12 EDDIS06 EDDIS12;
run;


/* Bayesian Two-stage modeling */
%macro BayesTwostage(OUTCOME,D1,D2,Drought);
data test;
  set NOAA.test;
  idc+1;
  by State County Year Month;
  if first.County then idc=1;
  Lpop=log(Pop);
run;

/* Stage1. Quasi-Poisson regression model */
ods exclude all;
proc glimmix data=test;
  by county;
  class Year;
  effect spl=spline(idc/ basis=bspline details);
  model &OUTCOME=&D1 &D2 spl Year TempA TempA*TempA / link=log offset=Lpop solution;
  _variance_ = _mu_;
  random _residual_;
  ods output ParameterEstimates=_para_est_USDM ConvergenceStatus=ConvergeStatus;
run;
ods exclude none;

data OUT;
  merge _para_est_USDM ConvergeStatus;
  by County;
  if Status=0 then output;
run;

proc sort data=OUT out=OUT_USDM_M;
  by Effect;
run;

data Moderate;
  set OUT_USDM_M(where=(Effect=&Drought));
  SE2=0;
  if Estimate<=-2.0 or Estimate>=1.0 then delete;
  if StdErr>=1.0 then delete;
  Est=round(Estimate,0.001);
  SE=round(StdErr,0.001);
  SE2=SE**2;
  if SE=.  then delete;
  keep County Est SE SE2 Probt;
run;

/* Stage2. Bayesian linear regression */
proc mcmc data=Moderate outpost=postoutB nmc=50000 nbi=10000 nthin=5 seed=1234 propcov=quanew DIC maxtune=30;
 
  parms alpha phi;
  prior alpha~normal(0,var=100);
  prior phi~normal(0,var=5);
  random muK~normal(0,var=1) subject=County;
  xalpha=alpha+muK*exp(phi);
  model Est~normal(xalpha,var=SE2);

run;

data NOAA.out;
  set postoutB;
  ORR=exp(alpha);
run;

proc means data=NOAA.out mean stddev p5 p95 min max maxdec=3;
  var ORR;
run;

%mend;

%BayesTwostage(OUTCOME=CD,D1=Dmod,D2=Dsev,Drought="Dmod");
%BayesTwostage(OUTCOME=CD,D1=Dmod,D2=Dsev,Drought="Dsev");
%BayesTwostage(OUTCOME=CD,D1=SPEIM06,D2=SPEIS06,Drought="SPEIM06");
%BayesTwostage(OUTCOME=CD,D1=SPEIM06,D2=SPEIS06,Drought="SPEIS06");
%BayesTwostage(OUTCOME=CD,D1=SPEIM12,D2=SPEIS12,Drought="SPEIM12");
%BayesTwostage(OUTCOME=CD,D1=SPEIM12,D2=SPEIS12,Drought="SPEIS12");
%BayesTwostage(OUTCOME=CD,D1=EDDIM06,D2=EDDIS06,Drought="EDDIM06");
%BayesTwostage(OUTCOME=CD,D1=EDDIM06,D2=EDDIS06,Drought="EDDIS06");
%BayesTwostage(OUTCOME=CD,D1=EDDIM12,D2=EDDIS12,Drought="EDDIM12");
%BayesTwostage(OUTCOME=CD,D1=EDDIM12,D2=EDDIS12,Drought="EDDIS12");

/************************/
/* Power Prior Approach */
/************************/
%macro BayesPowerPrior(OUTCOME,D1,D2,Drought);
data test;
  set NOAA.test;
  idc+1;
  by State County Year Month;
  if first.County then idc=1;
  Lpop=log(Pop);
run;

proc means data=test min p25 p50 p75 max mean maxdec=3;
  class Year;
  var TD CD RD Drought SPEIB06 SPEIB12 Dmod Dsev SPEIM06 SPEIM12 SPEIS06 SPEIS12;
run;

/* Construction of historical data */
data NOAA.hist;
  set test;
  if Year<=2011 then output;
run;

ods exclude all;
proc glimmix data=NOAA.hist;
  by county;
  class Year;
  effect spl=spline(idc/ basis=bspline details);
  model &OUTCOME=&D1 &D2 spl Year TempA TempA*TempA / link=log solution;
  _variance_ = _mu_;
  random _residual_;
  ods output ParameterEstimates=_para_est_USDM ConvergenceStatus=ConvergeStatus;
run;
ods exclude none;

data OUT;
  merge _para_est_USDM ConvergeStatus;
  by County;
  if Status=0 then output;
run;

proc sort data=OUT out=OUT_USDM_M;
  by Effect;
run;

data NOAA.hist;
  set OUT_USDM_M(where=(Effect=&Drought));
  Group=1;
  SE2=0;
  if Estimate<=-2.0 or Estimate>=1.0 then delete;
  if StdErr>=1.0 then delete;
  Est=round(Estimate,0.001);
  SE=round(StdErr,0.001);
  if SE=.  then delete;
  SE2=SE**2;
  keep County Est SE SE2 Probt Group;
run;

/* Construction of current data */
data NOAA.curr;
  set test;
  if Year>=2012 then output;
run;

ods exclude all;
proc glimmix data=NOAA.curr;
  by county;
  class Year;
  effect spl=spline(idc/ basis=bspline details);
  model &OUTCOME=&D1 &D2 spl Year TempA TempA*TempA / link=log solution;
  _variance_ = _mu_;
  random _residual_;
  ods output ParameterEstimates=_para_est_USDM ConvergenceStatus=ConvergeStatus;
run;
ods exclude none;

data OUT;
  merge _para_est_USDM ConvergeStatus;
  by County;
  if Status=0 then output;
run;

proc sort data=OUT out=OUT_USDM_M;
  by Effect;
run;

data NOAA.curr;
  set OUT_USDM_M(where=(Effect=&Drought));
  Group=2;
  SE2=0;
  if Estimate<=-2.0 or Estimate>=1.0 then delete;
  if StdErr>=1.0 then delete;
  Est=round(Estimate,0.001);
  SE=round(StdErr,0.001);
  if SE=.  then delete;
  SE2=SE**2;
  keep County Est SE SE2 Probt Group;
run;

data NOAA.final;
  set NOAA.hist NOAA.curr;
run;

/* Stage2. Bayesian power prior method */
proc mcmc data=NOAA.final outpost=postoutB nmc=50000 nbi=10000 nthin=5 seed=1234 propcov=quanew DIC maxtune=30;
 
  parms alpha phi a0;
  prior alpha~normal(0,var=100);
  prior phi~normal(0,var=5);
  prior a0~beta(1,1);

  if Group=1 then do;
    xalpha=a0*alpha;
	SEA=SE2+exp(phi);
  end;
  else if Group=2 then do;
    xalpha=alpha+muK*exp(phi);
    SEA=SE2;
  end;

  random muK~normal(0,var=1) subject=County;
  model Est~normal(xalpha,var=SEA);

run;

data NOAA.out;
  set postoutB;
  ORR=exp(alpha);
  A0=a0;
run;

proc means data=NOAA.out mean stddev p5 p95 min max maxdec=3;
  var ORR;
run;

%mend;

%BayesPowerPrior(OUTCOME=CD,D1=Dmod,D2=Dsev,Drought="Dmod");
%BayesPowerPrior(OUTCOME=CD,D1=Dmod,D2=Dsev,Drought="Dsev");
%BayesPowerPrior(OUTCOME=CD,D1=SPEIM06,D2=SPEIS06,Drought="SPEIM06");
%BayesPowerPrior(OUTCOME=CD,D1=SPEIM06,D2=SPEIS06,Drought="SPEIS06");
%BayesPowerPrior(OUTCOME=CD,D1=SPEIM12,D2=SPEIS12,Drought="SPEIM12");
%BayesPowerPrior(OUTCOME=CD,D1=SPEIM12,D2=SPEIS12,Drought="SPEIS12");
%BayesPowerPrior(OUTCOME=CD,D1=EDDIM06,D2=EDDIS06,Drought="EDDIM06");
%BayesPowerPrior(OUTCOME=CD,D1=EDDIM06,D2=EDDIS06,Drought="EDDIS06");
%BayesPowerPrior(OUTCOME=CD,D1=EDDIM12,D2=EDDIS12,Drought="EDDIM12");
%BayesPowerPrior(OUTCOME=CD,D1=EDDIM12,D2=EDDIS12,Drought="EDDIS12");

/* Run the model for different outcome */
%BayesTwostage(OUTCOME=RD,D1=Dmod,D2=Dsev,Drought="Dmod");
%BayesTwostage(OUTCOME=RD,D1=Dmod,D2=Dsev,Drought="Dsev");
%BayesTwostage(OUTCOME=RD,D1=SPEIM06,D2=SPEIS06,Drought="SPEIM06");
%BayesTwostage(OUTCOME=RD,D1=SPEIM06,D2=SPEIS06,Drought="SPEIS06");
%BayesTwostage(OUTCOME=RD,D1=SPEIM12,D2=SPEIS12,Drought="SPEIM12");
%BayesTwostage(OUTCOME=RD,D1=SPEIM12,D2=SPEIS12,Drought="SPEIS12");
%BayesTwostage(OUTCOME=RD,D1=EDDIM06,D2=EDDIS06,Drought="EDDIM06");
%BayesTwostage(OUTCOME=RD,D1=EDDIM06,D2=EDDIS06,Drought="EDDIS06");
%BayesTwostage(OUTCOME=RD,D1=EDDIM12,D2=EDDIS12,Drought="EDDIM12");
%BayesTwostage(OUTCOME=RD,D1=EDDIM12,D2=EDDIS12,Drought="EDDIS12");

%BayesPowerPrior(OUTCOME=RD,D1=Dmod,D2=Dsev,Drought="Dmod");
%BayesPowerPrior(OUTCOME=RD,D1=Dmod,D2=Dsev,Drought="Dsev");
%BayesPowerPrior(OUTCOME=RD,D1=SPEIM06,D2=SPEIS06,Drought="SPEIM06");
%BayesPowerPrior(OUTCOME=RD,D1=SPEIM06,D2=SPEIS06,Drought="SPEIS06");
%BayesPowerPrior(OUTCOME=RD,D1=SPEIM12,D2=SPEIS12,Drought="SPEIM12");
%BayesPowerPrior(OUTCOME=RD,D1=SPEIM12,D2=SPEIS12,Drought="SPEIS12");
%BayesPowerPrior(OUTCOME=RD,D1=EDDIM06,D2=EDDIS06,Drought="EDDIM06");
%BayesPowerPrior(OUTCOME=RD,D1=EDDIM06,D2=EDDIS06,Drought="EDDIS06");
%BayesPowerPrior(OUTCOME=RD,D1=EDDIM12,D2=EDDIS12,Drought="EDDIM12");
%BayesPowerPrior(OUTCOME=RD,D1=EDDIM12,D2=EDDIS12,Drought="EDDIS12");

/***********************/
/* Rural areas fitting */
/***********************/
%macro NOAA2stageRural(OUTCOME,D1,D2,V1,V2,OUTPUT1,OUTPUT2);

data tests0;
  set NOAA.test;
  if Urban=0 then output;
run;

ods exclude all;
proc glimmix data=tests0;
  by State County;
  class Year;
  effect spl=spline(Month / naturalcubic knotmethod=equal(3));
  model &OUTCOME = &D1 &D2 spl Year TempA TempA*TempA / link = log solution;
  _nu =  1 / exp(_phi_);
  _variance_ = (1 / _nu) / ((_mu_) ** (1 / _nu));
  _z = 0;
  do i = 0 to 100;
    _z = _z + (_mu_ ** i) / fact(i) ** _nu;
  end;
  _prob = (_mu_ ** &OUTCOME) / (fact(&OUTCOME) ** _nu) * (_z ** (-1));
  _logl_ = log(_prob);
  ods output ParameterEstimates=_para_est_S0 ConvergenceStatus=ConvergeStatus;
run;
ods exclude none;

data OUT;
  merge _para_est_S0 ConvergeStatus;
  by County;
  if Status=0;
run;

proc sort data=OUT out=OUTS0;
  by Effect;
run;

data Moderate;
  set OUTS0(where=(Effect=&V1));
  SE2=0;
  if Estimate<=-2.0 or Estimate>=0.7 then delete;
  if StdErr>=1.0 then delete;
  Est=round(Estimate,0.001);
  SE=round(StdErr,0.001);
  SE2=SE**2;
  if SE=0  then delete;
  keep State County Est SE SE2 Probt;
run;

data Severe;
  set OUTS0(where=(Effect=&V2));
  SE2=0;
  if Estimate<=-2.0 or Estimate>=0.7 then delete;
  if StdErr>=1.0 then delete;
  Est=round(Estimate,0.001);
  SE=round(StdErr,0.001);
  SE2=SE**2;
  if SE=0  then delete;
  keep State County Est SE SE2 Probt;
run;

/*************************/
/* Stage2. Meta Analysis */
/*************************/
ods graphics on; 
proc mcmc data=Moderate outpost=postoutMod plots=none diagnostics=none stats=all
		  nmc=100000 nbi=50000 nthin=5 seed=1234 propcov=quanew DIC maxtune=30;

  by State;
  ods output PostSummaries=PostSummaries; 
  ods output PostIntervals=PostIntervals; 
  parms alpha phi;
  prior alpha~normal(0,var=100);
  prior phi~normal(0,var=5);
  random muK~normal(0,var=1) subject=County;

  xalpha=alpha+muK*exp(phi);
  model Est~normal(xalpha,var=SE2);

run;
ods graphics off;

data test(keep=State Parameter Mean StdDev HPDLower HPDUpper);
  merge Postsummaries Postintervals;
  by State Parameter;
run;

proc sort data=test out=test2;
  by Parameter;
run;

proc export data=test2
  outfile="Urban0-CMP-&OUTCOME.-BayesMetaState-&OUTPUT1..csv"
  dbms=csv
  replace;
run;

ods graphics on; 
proc mcmc data=Severe outpost=postoutSev plots=none diagnostics=none stats=all
		  nmc=100000 nbi=50000 nthin=5 seed=1234 propcov=quanew DIC maxtune=30;

  by State;
  ods output PostSummaries=PostSummaries; 
  ods output PostIntervals=PostIntervals; 
  parms alpha phi;
  prior alpha~normal(0,var=100);
  prior phi~normal(0,var=5);
  random muK~normal(0,var=1) subject=County;

  xalpha=alpha+muK*exp(phi);
  model Est~normal(xalpha,var=SE2);

run;
ods graphics off;

data test(keep=State Parameter Mean StdDev HPDLower HPDUpper);
  merge Postsummaries Postintervals;
  by State Parameter;
run;

proc sort data=test out=test2;
  by Parameter;
run;

proc export data=test2
  outfile="Urban0-CMP-&OUTCOME.-BayesMetaState-&OUTPUT2..csv"
  dbms=csv
  replace;
run;

%mend;

%NOAA2stageRural(OUTCOME=RD,D1=Dmod,D2=Dsev,V1="Dmod",V2="Dsev",OUTPUT1=USDMmoderate,OUTPUT2=USDMsevere);
%NOAA2stageRural(OUTCOME=CD,D1=Dmod,D2=Dsev,V1="Dmod",V2="Dsev",OUTPUT1=USDMmoderate,OUTPUT2=USDMsevere);
%NOAA2stageRural(OUTCOME=TD,D1=Dmod,D2=Dsev,V1="Dmod",V2="Dsev",OUTPUT1=USDMmoderate,OUTPUT2=USDMsevere);
%NOAA2stageRural(OUTCOME=RD,D1=SPEIM06,D2=SPEIS06,V1="SPEIM06",V2="SPEIS06",OUTPUT1=SPEI06moderate,OUTPUT2=SPEI06severe);
%NOAA2stageRural(OUTCOME=CD,D1=SPEIM06,D2=SPEIS06,V1="SPEIM06",V2="SPEIS06",OUTPUT1=SPEI06moderate,OUTPUT2=SPEI06severe);
%NOAA2stageRural(OUTCOME=TD,D1=SPEIM06,D2=SPEIS06,V1="SPEIM06",V2="SPEIS06",OUTPUT1=SPEI06moderate,OUTPUT2=SPEI06severe);
%NOAA2stageRural(OUTCOME=RD,D1=SPEIM12,D2=SPEIS12,V1="SPEIM12",V2="SPEIS12",OUTPUT1=SPEI12moderate,OUTPUT2=SPEI12severe);
%NOAA2stageRural(OUTCOME=CD,D1=SPEIM12,D2=SPEIS12,V1="SPEIM12",V2="SPEIS12",OUTPUT1=SPEI12moderate,OUTPUT2=SPEI12severe);
%NOAA2stageRural(OUTCOME=TD,D1=SPEIM12,D2=SPEIS12,V1="SPEIM12",V2="SPEIS12",OUTPUT1=SPEI12moderate,OUTPUT2=SPEI12severe);

/***********************/
/* Urban areas fitting */
/***********************/
%macro NOAA2stageUrban(OUTCOME,D1,D2,V1,V2,OUTPUT1,OUTPUT2);

data tests1;
  set NOAA.test;
  if Urban=1 then output;
run;

ods exclude all;
proc glimmix data=tests1;
  by State County;
  class Year;
  effect spl=spline(Month / naturalcubic knotmethod=equal(3));
  model &OUTCOME = &D1 &D2 spl Year TempA TempA*TempA / link = log solution;
  _nu =  1 / exp(_phi_);
  _variance_ = (1 / _nu) / ((_mu_) ** (1 / _nu));
  _z = 0;
  do i = 0 to 100;
    _z = _z + (_mu_ ** i) / fact(i) ** _nu;
  end;
  _prob = (_mu_ ** &OUTCOME) / (fact(&OUTCOME) ** _nu) * (_z ** (-1));
  _logl_ = log(_prob);
  ods output ParameterEstimates=_para_est_S1 ConvergenceStatus=ConvergeStatus;
run;
ods exclude none;

data OUT;
  merge _para_est_S1 ConvergeStatus;
  by County;
  if Status=0;
run;

proc sort data=OUT out=OUTS1;
  by Effect;
run;

data Moderate;
  set OUTS1(where=(Effect=&V1));
  SE2=0;
  if Estimate<=-2.0 or Estimate>=0.7 then delete;
  if StdErr>=1.0 then delete;
  Est=round(Estimate,0.001);
  SE=round(StdErr,0.001);
  SE2=SE**2;
  if SE=0  then delete;
  keep State County Est SE SE2 Probt;
run;

data Severe;
  set OUTS1(where=(Effect=&V2));
  SE2=0;
  if Estimate<=-2.0 or Estimate>=0.7 then delete;
  if StdErr>=1.0 then delete;
  Est=round(Estimate,0.001);
  SE=round(StdErr,0.001);
  SE2=SE**2;
  if SE=0  then delete;
  keep State County Est SE SE2 Probt;
run;

/*************************/
/* Stage2. Meta Analysis */
/*************************/
ods graphics on; 
proc mcmc data=Moderate outpost=postoutMod plots=none diagnostics=none stats=all
		  nmc=100000 nbi=50000 nthin=5 seed=1234 propcov=quanew DIC maxtune=30;

  by State;
  ods output PostSummaries=PostSummaries; 
  ods output PostIntervals=PostIntervals; 
  parms alpha phi;
  prior alpha~normal(0,var=100);
  prior phi~normal(0,var=5);
  random muK~normal(0,var=1) subject=County;

  xalpha=alpha+muK*exp(phi);
  model Est~normal(xalpha,var=SE2);

run;
ods graphics off;

data test(keep=State Parameter Mean StdDev HPDLower HPDUpper);
  merge Postsummaries Postintervals;
  by State Parameter;
run;

proc sort data=test out=test2;
  by Parameter;
run;

proc export data=test2
  outfile="Urban1-CMP-&OUTCOME.-BayesMetaState-&OUTPUT1..csv"
  dbms=csv
  replace;
run;

ods graphics on; 
proc mcmc data=Severe outpost=postoutSev plots=none diagnostics=none stats=all
		  nmc=100000 nbi=50000 nthin=5 seed=1234 propcov=quanew DIC maxtune=30;

  by State;
  ods output PostSummaries=PostSummaries; 
  ods output PostIntervals=PostIntervals; 
  parms alpha phi;
  prior alpha~normal(0,var=100);
  prior phi~normal(0,var=5);
  random muK~normal(0,var=1) subject=County;

  xalpha=alpha+muK*exp(phi);
  model Est~normal(xalpha,var=SE2);

run;
ods graphics off;

data test(keep=State Parameter Mean StdDev HPDLower HPDUpper);
  merge Postsummaries Postintervals;
  by State Parameter;
run;

proc sort data=test out=test2;
  by Parameter;
run;

proc export data=test2
  outfile="Urban1-CMP-&OUTCOME.-BayesMetaState-&OUTPUT2..csv"
  dbms=csv
  replace;
run;

%mend;

%NOAA2stageUrban(OUTCOME=RD,D1=Dmod,D2=Dsev,V1="Dmod",V2="Dsev",OUTPUT1=USDMmoderate,OUTPUT2=USDMsevere);
%NOAA2stageUrban(OUTCOME=CD,D1=Dmod,D2=Dsev,V1="Dmod",V2="Dsev",OUTPUT1=USDMmoderate,OUTPUT2=USDMsevere);
%NOAA2stageUrban(OUTCOME=TD,D1=Dmod,D2=Dsev,V1="Dmod",V2="Dsev",OUTPUT1=USDMmoderate,OUTPUT2=USDMsevere);
%NOAA2stageUrban(OUTCOME=RD,D1=SPEIM06,D2=SPEIS06,V1="SPEIM06",V2="SPEIS06",OUTPUT1=SPEI06moderate,OUTPUT2=SPEI06severe);
%NOAA2stageUrban(OUTCOME=CD,D1=SPEIM06,D2=SPEIS06,V1="SPEIM06",V2="SPEIS06",OUTPUT1=SPEI06moderate,OUTPUT2=SPEI06severe);
%NOAA2stageUrban(OUTCOME=TD,D1=SPEIM06,D2=SPEIS06,V1="SPEIM06",V2="SPEIS06",OUTPUT1=SPEI06moderate,OUTPUT2=SPEI06severe);
%NOAA2stageUrban(OUTCOME=RD,D1=SPEIM12,D2=SPEIS12,V1="SPEIM12",V2="SPEIS12",OUTPUT1=SPEI12moderate,OUTPUT2=SPEI12severe);
%NOAA2stageUrban(OUTCOME=CD,D1=SPEIM12,D2=SPEIS12,V1="SPEIM12",V2="SPEIS12",OUTPUT1=SPEI12moderate,OUTPUT2=SPEI12severe);
%NOAA2stageUrban(OUTCOME=TD,D1=SPEIM12,D2=SPEIS12,V1="SPEIM12",V2="SPEIS12",OUTPUT1=SPEI12moderate,OUTPUT2=SPEI12severe);

