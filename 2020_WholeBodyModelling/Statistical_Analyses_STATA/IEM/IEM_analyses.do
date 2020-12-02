*STATA script performing the analyses of prediction accuracies of models regarding biomarker in IEMs

*Johannes Hertel, Oct 2018, email: johannes.hertel@nuigalway.ie

version 14
clear
clear mata
clear matrix
set more off
capture log close
set maxvar 32000

cd A:\Hertel\Finalized\Thiele_2018_HH\Statistical_Analyses_STATA\IEM //adapt directory as suitable

log using "prediction_accuracy_analyses_080520.log", replace

*******************************************
*Data preparation: Import from Excel-files*
*******************************************

foreach j in REC HTA HV{
	import excel "A:\Hertel\Finalized\Thiele_2018_HH\Statistical_Analyses_STATA\IEM\IEM_predictions_Overview_rc_08_05_20.xlsx", firstrow //adapt directory as suitable
	egen id=seq(), by(G)
	drop G N
	keep id *`j'
	gen iem_bio=IEM+Biomarker
	egen test=seq(), by(iem_bio)
	drop if test>1
	drop test
	save `j'_pred.dta, replace
	gen agree_`j'=insilico_`j'-invivo_`j'
	replace agree_`j'=1 if agree_`j'==-1 | agree_`j'==-2 | agree_`j'==2
	recode agree_`j' (0=1) (1=0)
	tab agree_`j'
	egen mean_agree_`j'=mean(agree_`j'), by(IEM_`j')
	save `j'_pred.dta, replace
	collapse agree_`j',by(IEM_`j')
	rename IEM_`j' IEM
	save `j'_IEM.dta, replace
	clear
	}
	
**************************************************************
*Analyses of agreements on the level of reactions/metabolites*
**************************************************************	

use REC_pred.dta
merge 1:1 iem_bio using HV_pred.dta
rename _merge merge1
merge 1:1 iem_bio using HTA_pred.dta
rename _merge merge2

*Cohen's Kappa: Metric of agreement between in-vivo measurements and prediction from modelling

foreach j in REC HTA HV{
	tab invivo_`j' insilico_`j'
	kap  invivo_`j' insilico_`j'
	local CI_l=r(kappa)-invnormal(0.975)*r(se)
	local CI_u=r(kappa)+invnormal(0.975)*r(se)
	display "95%-CI lower bound:"
	display `CI_l'
	display "95%-CI upper bound:"
	display `CI_u'
	}

*Model-fit and p-values from multinomial logistic regressions

foreach j in REC HTA HV{
	replace insilico_`j'=insilico_`j'+1
	mlogit invivo_`j' i.insilico_`j',
	display "p-value:"
	display chi2tail(2,e(chi2))
	}
clear
	
*********************************************
*Analyses of agreements on the level of IEMs*
*********************************************	

use REC_IEM.dta
merge 1:1 IEM using HV_IEM.dta
rename _merge merge1
merge 1:1 IEM using HTA_IEM.dta
rename _merge merge2

*Non-parametric sign test comparing agreement percentages in IEMs for every pair of model
signtest agree_REC=agree_HV
signtest agree_REC=agree_HTA
signtest agree_HV=agree_HTA

*IEMs with REC better HV
tab IEM if agree_REC>agree_HV & agree_REC!=.
*IEMs with REC better HTA
tab IEM if agree_REC>agree_HTA & agree_REC!=.
*IEMs with REC worse HV
tab IEM if agree_REC<agree_HV & agree_REC!=.
*IEMs with REC worse HTA
tab IEM if agree_REC<agree_HTA & agree_REC!=.
*IEMs with REC equal HV
tab IEM if agree_REC==agree_HV & agree_REC!=.
*IEMs with REC equal HTA
tab IEM if agree_REC==agree_HTA & agree_REC!=.

log close

clear
