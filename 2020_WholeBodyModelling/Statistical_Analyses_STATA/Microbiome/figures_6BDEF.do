*STATA script performing the analyses related to Figure 6B/D/E/F

*Johannes Hertel, Nov/Dec 2018, email: johannes.hertel@uni.lu

clear
clear mata
clear matrix
set more off
set maxvar 32767
capture log close

cd A:\Hertel\Finalized\Thiele_2018_HH\Microbiome_modeling\
log using figures_6BDEF.log, replace

import delimited A:\Hertel\Finalized\Thiele_2018_HH\Microbiome_modeling\bacteroidia_clostridia_ratio.csv
save bacteroidia_clostridia_ratio.dta, replace
clear
import delimited A:\Hertel\Finalized\Thiele_2018_HH\Microbiome_modeling\fold_changes_modeling.csv
rename v1 sample
save fold_changes_modeling.dta, replace
merge 1:1 sample using bacteroidia_clostridia_ratio.dta
drop _merge

*Liver Sulphotransferase
mfp reg liver_pcsf bacteroidia_clostridia_ratio
*p-value:
display Ftail(e(df_m), e(df_r),e(F))
predict sulpho_score
*graph
cd A:\Hertel\Finalized\Thiele_2018_HH\Microbiome_modeling\figures\raw
scatter liver_pcsf bacteroidia_clostridia_ratio || fpfit sulpho_score bacteroidia_clostridia_ratio, saving(figure_6D, replace)

*Liver alcohol dehydrogenase
mfp reg liver_alcd2if bacteroidia_clostridia_ratio
*p-value:
display Ftail(e(df_m), e(df_r),e(F))
predict alco_score
*graph
cd A:\Hertel\Finalized\Thiele_2018_HH\Microbiome_modeling\figures\raw
scatter liver_alcd2if bacteroidia_clostridia_ratio || fpfit alco_score bacteroidia_clostridia_ratio, saving(figure_6E, replace)

*Colon carboxylic acid:CoA ligase
mfp reg colon_hmr_0156 bacteroidia_clostridia_ratio
*p-value:
display Ftail(e(df_m), e(df_r),e(F))
predict carbo_CoA_ligase
*graph
cd A:\Hertel\Finalized\Thiele_2018_HH\Microbiome_modeling\figures\raw
scatter colon_hmr_0156 bacteroidia_clostridia_ratio || fpfit carbo_CoA_ligase bacteroidia_clostridia_ratio, saving(figure_6F, replace)

clear
import excel "A:\Hertel\Finalized\Thiele_2018_HH\Microbiome_modeling\AbsoluteFluxes_Microbiota_Germfree.xlsx", sheet("Sheet2") cellrange(A2:W149) firstrow
rename A sample
gen d_Braindopamineproduction=-(BraindopamineproductionG-Braindopamineproduction)
gen d_Brainserotoninproduction=-BrainserotoninproductionG+Brainserotoninproduction
gen d_Brainadrenalineproduction=-BrainadrenalineproductionG+Brainadrenalineproduction
gen d_BrainGABAproduction=-BrainGABAproductionG+BrainGABAproduction
gen d_Brainhistamineproduction=-BrainhistamineproductionG+Brainhistamineproduction
gen d_Brainkynurenicacid=-BrainkynurenicacidG+Brainkynurenicacidproduction
gen d_Brainnorepinephrine=-BrainnorepinephrineG+Brainnorepinephrineproduction
gen d_BrainLkynurenine=-BrainnorepinephrineG+BrainLkynurenineproduction

ttest BraindopamineproductionG=Braindopamineproduction
*Bonferoni correction:
display r(p)*8
ttest BrainserotoninproductionG=Brainserotoninproduction
*Bonferoni correction:
display r(p)*8
ttest BrainadrenalineproductionG=Brainadrenalineproduction
*Bonferoni correction:
display r(p)*8
ttest BrainkynurenicacidG=Brainkynurenicacidproduction
*Bonferoni correction:
display r(p)*8
ttest BrainGABAproductionG=BrainGABAproduction
*Bonferoni correction:
display r(p)*8
ttest BrainhistamineproductionG=Brainhistamineproduction
*Bonferoni correction:
display r(p)*8
ttest BrainnorepinephrineG=Brainnorepinephrineproduction
*Bonferoni correction:
display r(p)*8
ttest BrainnorepinephrineG=BrainLkynurenineproduction
*Bonferoni correction:
display r(p)*8
local i=1
label define met 1 "Dopamine" 2 "Serotonin" 3 "Adrenaline" 4 "GABA" 5 "Histamine" 6 "Kynurenic acid" 7 "Norepinephrine" 8 "L-Kynurenine"
gen var_name=.
label values var_name met

gen mean_=.
gen CI_l=.
gen CI_u=.
gen p_val=.
foreach j of varlist d_*{
	mean `j'
	matrix V=e(V)
	matrix B=e(b)
	replace mean=B[1,1] in `i'
	replace CI_l=B[1,1]+invttail(e(df_r),0.975)*sqrt(V[1,1]) in `i'
	replace CI_u=B[1,1]-invttail(e(df_r),0.975)*sqrt(V[1,1]) in `i'
	replace var_name=`i' in `i'
	local i=`i'+1
	}
cd A:\Hertel\Finalized\Thiele_2018_HH\Microbiome_modeling\figures\raw	
graph twoway bar mean_ var_name, horizontal || rcap CI_u CI_l var_name, horizontal saving(figure4B_raw, replace)

clear

import excel "A:\Hertel\Finalized\Thiele_2018_HH\Microbiome_modeling\Flux_Sulfotransferase_vs_Clostridioides.xlsx", sheet("Sheet1") firstrow
rename A sample
mfp reg LiverSulfotransferase Clostridioides
display Ftail(e(df_m), e(df_r),e(F))

*Graph S3
cd A:\Hertel\Finalized\Thiele_2018_HH\Microbiome_modeling\figures\raw
scatter LiverSulfotransferase Clostridioides || function y=.0009936*(((x+2.00000002337e-07)*100)^(-1)-12.49682479)-.2439572*(((x+2.00000002337e-07)*100)^(-.5)-3.535084834)+4.024158, range(0 0.006) saving(figure_S3, replace)

clear
log close
