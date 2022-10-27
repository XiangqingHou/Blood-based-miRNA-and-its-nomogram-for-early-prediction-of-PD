proc import datafile="D:\Desk\work\work_PPMI\meta_data.xlsx" dbms=xlsx out=work;quit;
proc import datafile="D:\Desk\work\work_PPMI\smallRNA.supp.3_ESM.xlsx" dbms=xlsx out=work1;quit;
data work_1;set work;keep HudAlphaSampleName PATNO PATNO_visit sex month ethnicity race education_level;run;
proc sort data=work_1;by HudAlphaSampleName;run;
proc sort data=work1(rename=HudsonAlpha_sample_name=HudAlphaSampleName);by HudAlphaSampleName;run;
data work2;merge work1(in=aa) work_1;by HudAlphaSampleName;if aa;
drop HudsonAlpha_library_ID Sample_ID Sequencing_batch_ID Sequencing_plate_number HudAlphaSampleName 
Most_likely_primary_diagnosis Age_at_diagnosis Age_binning UPDRS_score_IV
Visit_prior_medication_use  Visit_post_medication_use PATNO;
run;
/**********************************************1444 patients contribute to 4171 sample*/
data work3;set work2;
age=age_at_consent+0;
score_t=UPDRS_score_total+0;
score1=UPDRS_score_I+0;
score2=UPDRS_score_II+0;
score3=UPDRS_score_III+0;
stage=Hoehn_and_Yahr_staging+0;
if Visit='BL' then time=0;
if Visit='V02' then time=6;
if Visit='V04' then time=12;
if Visit='V06' then time=24;
if Visit='V08' then time=36;
if Biogroup='Prodromal' or Biogroup='Control' or Biogroup='PD';
drop age_at_consent UPDRS_score_total UPDRS_score_I UPDRS_score_II UPDRS_score_III;
if missing(score_t) then delete;/*1416 to 3806*/
run;
%idrpt(data=work3,id=Patient_number)
/**************************keep patients with baseline information PD HY stage in 0,1,2***************************/
data work_bl;set work3;
if Visit='BL';
if Biogroup='PD' and stage=3 then delete;
if Biogroup='PD' and stage=4 then delete;
if Biogroup='PD' and stage=5 then delete;
keep Patient_number;
run;
proc sort data=work_bl;by Patient_number;
proc sort data=work3;by Patient_number;
data work4;merge work_bl(in=aa) work3;by Patient_number;if aa;run;
/*1284 patients contribute to 3550 visits*/
proc sort data=work4;by time;run;
proc freq data=work4;tables stage;by time;run;
data work5;set work4;
if Ethnicity='Not Hispanic or Latino' then Ethnicity1=1;
if Ethnicity='Hispanic or Latino' then Ethnicity1=2;
if missing(Ethnicity1) then Ethnicity1=3;
if Education_Level='Less than 12 years' or Education_Level='0 years' then Education_Level1=1;
if Education_Level='12-16 years' then Education_Level1=2;
if Education_Level='Greater than 16 years' then Education_Level1=3;
if missing(Education_Level1) then Education_Level1=4;
if gender='Male' then sex1=1;else sex1=2;
if Genetic_status='GBA+' then gene=1;
if Genetic_status='LRRK2+' then gene=2;
if Genetic_status='LRRK2-/SNCA-/GBA-' then gene=3;
if Genetic_status='SNCA+' then gene=4;
if missing(Genetic_status) then gene=5;
patno=tranwrd(Sequencing_sample_ID,"-",".");
if Biogroup='Control' then outcome=0;else outcome=1;
keep Patient_number patno Ethnicity1 Education_Level1 sex1 gene age outcome Clinical_site score_t stage time Biogroup;/**/
run;
proc sort data=work5;by Clinical_site;run;
data split_id;set work5;by Clinical_site;if first.Clinical_site;keep Clinical_site;run;
data model1_1;set split_id;run;
%macro num(train_pct,validate_pct);
proc sql noprint;
select count(*) into :t_nobs from model1_1;
quit;
%let train_obs=%sysevalf(&train_pct.*&t_nobs.);
%let validate_obs=%sysevalf(&validate_pct.*&t_nobs.);
%put &train_obs.;
%put &validate_obs.;
%mend;
%num(.7,.3);
%macro partition(train_obs,validate_obs);
proc surveyselect data=model1_1 out=split seed=123
group=(&train_obs.,&validate_obs.);
run;
data training validate;set split;
if groupid=1 then do;drop groupid;output training;end;
if groupid=2 then do;drop groupid;output validate;end;
run;
%mend;
%partition(23,10);
proc sort data=training;by Clinical_site;
proc sort data=Validate;by Clinical_site;
proc sort data=work5;by Clinical_site;
data tmp_train;merge training(in=aa) work5;by Clinical_site;if aa;train=1;run;/*2211*/
data clinical;set tmp_train;keep patient_number score_t stage time patno outcome;run;
/*Perform feature selection using R software*/
/*Merge clinical datasets and difference miRNA datasets*/
proc import datafile="D:\Desk\Figure_ROC\diff_miRNA.csv" dbms=csv out=diff_miRNA;quit;
data diff_miRNA1;set diff_miRNA;id=scan(patno,3,'.')||'_'||scan(patno,4,'.');run;
data tmp_train1;set tmp_train;id=scan(patno,3,'.')||'_'||scan(patno,4,'.');run;
data tmp_Validate1;set tmp_Validate;id=scan(patno,3,'.')||'_'||scan(patno,4,'.');run;
proc sort data=tmp_train1;by id;
proc sort data=tmp_Validate1;by id;
proc sort data=diff_miRNA1;by id;
data train;merge diff_miRNA1 tmp_train1(in=aa);by id;if aa;run;
data external;merge diff_miRNA1 tmp_Validate1(in=aa);by id;if aa;run;
data train1;set train;if time=0;age1=age;run;/*794*/
data External1;set External;if time=0;age1=age;run;/*490*/
/****************************************/
data model1_1;set train1;run;
%macro num(train_pct,validate_pct);
proc sql noprint;
select count(*) into :t_nobs from model1_1;
quit;
%let train_obs=%sysevalf(&train_pct.*&t_nobs.);
%let validate_obs=%sysevalf(&validate_pct.*&t_nobs.);
%put &train_obs.;
%put &validate_obs.;
%mend;
%num(.7,.3);
%macro partition(train_obs,validate_obs);
proc surveyselect data=model1_1 out=split seed=123
group=(&train_obs.,&validate_obs.);
run;
data training validate;set split;
if groupid=1 then do;drop groupid;output training;end;
if groupid=2 then do;drop groupid;output validate;end;
run;
%mend;
%partition(556,238);
/********************************/
/*training sets n=556*/
/*validation sets n=238*/
/*testing sets n=490*/
/********************************/
