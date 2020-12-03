*STATA script performing the analyses of prediction accuracies of BMRs

*Johannes Hertel, Nov 2018, email: johannes.hertel@uni.lu


clear
clear mata
clear matrix
set more off
capture log close

cd A:\Hertel\Finalized\Thiele_2018_HH\BMR

log using "BMR_prediction_analyses.log", replace

import excel "A:\Hertel\Finalized\Thiele_2018_HH\Statistical_Analyses_STATA\BMR\BMR_predictions_Oct_2018_1.xlsx", sheet("BMR calculation - parameter ide") firstrow cellrange(A3:P16) 

*calculate model-fit and p-values on the training data set n=13
foreach j of varlist Muscleatphydrolysisadjustment-P{
	reg MeasuredBMRkcal `j'
	test `j'
	display r(p)
	ci2 MeasuredBMRkcal `j', corr
	}

reg MeasuredBMRkcal P M
test P
display r(p)
test M
display r(p)
*Graphs:
cd A:\Hertel\Finalized\Thiele_2018_HH\BMR\figures\raw
*Training data set:
scatter MeasuredBMRkcal Muscleatphydrolysisadjustment || lfit MeasuredBMRkcal Muscleatphydrolysisadjustment, saving(fig_5A1.gph, replace)
scatter MeasuredBMRkcal M || lfit MeasuredBMRkcal M, saving(fig_5A2.gph, replace)
scatter MeasuredBMRkcal P || lfit MeasuredBMRkcal P, saving(fig_5A3.gph, replace)	
	
*calculate model-fit and p-values on the test data set n=	
clear

import excel "A:\Hertel\Finalized\Thiele_2018_HH\BMR\Kopie von BMR_predictions_Oct_2018_1.xlsx", sheet("BMR calculation - parameter val") firstrow cellrange(A5:K33)
destring(WBMbasedpredictionofBMR BasedonoriginalHarrisBened), replace force 
egen sex=group(Gender)

reg  BMRkcal  WBMbasedpredictionofBMR sex 	
predict WBM_score
ci2 WBM_score BMRkcal, corr
reg  BMRkcal  BasedonoriginalHarrisBened sex
predict Base_score
ci2 Base_score BMRkcal, corr

reg  BMRkcal BasedonoriginalHarrisBened sex WBMbasedpredictionofBMR
test WBMbasedpredictionofBMR
display r(p)
test BasedonoriginalHarrisBened 
display r(p)

*Graphs:
cd A:\Hertel\Finalized\Thiele_2018_HH\BMR\figures\raw
*Test data set
scatter BMRkcal WBM_score || lfit BMRkcal WBM_score, saving(fig_5B1.gph, replace)
scatter BMRkcal Base_score || lfit BMRkcal Base_score, saving(fig_5B2.gph, replace)
 


log close
clear
