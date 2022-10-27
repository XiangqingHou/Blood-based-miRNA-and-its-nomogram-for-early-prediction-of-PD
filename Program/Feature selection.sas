%let result=hsa_miR_4301 hsa_miR_106b_5p hsa_miR_18a_5p hsa_miR_190a_5p hsa_miR_21_5p hsa_miR_22_3p
hsa_miR_3200_5p hsa_miR_3613_5p hsa_miR_3928_3p hsa_miR_423_5p hsa_miR_4433b_5p hsa_miR_4677_5p
hsa_miR_485_3p hsa_miR_548b_5p hsa_miR_628_3p hsa_miR_654_5p hsa_miR_6815_5p;
ods html  gpath="D:\Desk\Figure_ROC" dpi=300;
ods graphics /width=14cm height=10cm  
outputfmt=jpg;
proc glmselect data=training plots=(criterion ase) seed=123;
class Education_Level1(ref='1');
model outcome=&result gender Education_Level1 age/selection=elasticnet(stop=none choose=aic) cvmethod=random(8);
output out=out2 p=predelasticnet;
run;
