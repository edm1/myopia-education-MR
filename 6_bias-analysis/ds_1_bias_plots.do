//Neil Davies 15/03/17, updated 01/08/17 for the full release
//This runs the bias component plots for the myopia project:

cap prog drop hetero_test
prog def hetero_test

args path outcome exposure iv 

macro shift 3
local covar="`*'"

cap drop _const
cap gen _const  = 1

di "outcome=`outcome'"
di "exposure=`exposure'"
di "instrumen=`iv'"

gmm (`outcome' - {xb1:`exposure'  _const})  ///
	(`outcome' - {xb2:`exposure'  _const}) , ///
	instruments(1:`exposure' ) ///
	instruments(2:`iv') ///
	winit(unadjusted,independent) onestep  ///
	vce(robust) ///
	deriv(1/xb1 = -1) ///
	deriv(2/xb2 = -1)
drop _const

local outcome2=substr("`outcome'",1,16)
est sto results_`outcome2'

lincom _b[xb1:`exposure']-_b[xb2:`exposure']
local het_p=2*(1-normal(abs(r(estimate)/r(se))))

regsave `exposure' using "`path'", detail(all) pval ci replace addvar(het_p,`het_p')

end

import delimited "/Volumes/Height_BMI_and_schooling/MR Myopia education/rawdata/myopia_education_toNeil_170728/phenotypes_alleleScores_170728.forNeil.tsv", encoding(ISO-8859-1)


destring tdi_log northing easting birthweight breastfed,replace force
gen male=(sex_genetic=="Male")
//Normalize continious variables

foreach i in northing easting tdi_log age birthweight ave_mse_clean eduyears_clean pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10{
	egen z_`i'=std(`i')
	}

foreach i in z_northing z_easting z_tdi_log z_age male breastfed z_birthweight z_pc1 z_pc2 z_pc3 z_pc4 z_pc5 z_pc6 z_pc7 z_pc8 z_pc9 z_pc10{
	hetero_test "results/myopia_bias_plots_`i'" `i' ave_mse_clean myopia_as_dosage 
	}

foreach i in z_northing z_easting z_tdi_log z_age male breastfed z_birthweight z_pc1 z_pc2 z_pc3 z_pc4 z_pc5 z_pc6 z_pc7 z_pc8 z_pc9 z_pc10{
	hetero_test "results/educ_bias_plots_`i'" `i' eduyears_clean ea_as_dosage 
	}

use "results/myopia_bias_plots_z_northing",clear
foreach i in  z_easting z_tdi_log z_age male breastfed z_birthweight z_pc1 z_pc2 z_pc3 z_pc4 z_pc5 z_pc6 z_pc7 z_pc8 z_pc9 z_pc10{
	append using "results/myopia_bias_plots_`i'"
	}
foreach i in z_northing z_easting z_tdi_log z_age male breastfed z_birthweight z_pc1 z_pc2 z_pc3 z_pc4 z_pc5 z_pc6 z_pc7 z_pc8 z_pc9 z_pc10{
	append using "results/educ_bias_plots_`i'"
	}

//Plot the results:
mkmat coef ci_lower ci_upper if strrpos(var,"mse")!=0 & strrpos(var,"xb2")!=0, mat(results_iv) rownames(_estimates_name )
mkmat coef ci_lower ci_upper if strrpos(var,"mse")!=0 & strrpos(var,"xb1")!=0, mat(results_ols) rownames(_estimates_name )

#delimit ;
coefplot (matrix(results_iv[,1]), ms(T) mc(blue) lc(blue) ciopts(lc(blue)) lc(blue) ci((results_iv[,2] results_iv[,3])) offset(.05)) ///
	     (matrix(results_ols[,1]), ms(S) mc(red) ciopts(lc(red)) lstyle(p1) lc(red) ci((results_ols[,2] results_ols[,3]))  offset(-.05)) , ///
	xtitle(Bias value)  ytitle(Covariate)  xline(0) leg(off) ylabel(,format(%9.1f))  graphregion(color(white)) plotregion(color(white)) grid(none) ///
	coeflabels(results_z_northing="Northing"
results_z_easting="Easting"
results_z_tdi_log="Log Townsend Deprivation Index"
results_z_age="Age"
results_male="Male"
results_breastfed="Breastfed"
results_z_birthweight="Birthweight"
results_z_pc1="1st genetic PC"
results_z_pc2="2nd genetic PC"
results_z_pc3="3rd genetic PC"
results_z_pc4="4th genetic PC"
results_z_pc5="5th genetic PC"
results_z_pc6="6th genetic PC"
results_z_pc7="7th genetic PC"
results_z_pc8="8th genetic PC"
results_z_pc9="9th genetic PC"
results_z_pc10="10th genetic PC");

graph export "results/bias_mse.eps",replace fontface("Times New Roman");


mkmat coef ci_lower ci_upper if strrpos(var,"edu")!=0 & strrpos(var,"xb2")!=0, mat(results_iv) rownames(_estimates_name );
mkmat coef ci_lower ci_upper if strrpos(var,"edu")!=0 & strrpos(var,"xb1")!=0, mat(results_ols) rownames(_estimates_name );
	
coefplot (matrix(results_iv[,1]), ms(T) mc(blue) lc(blue) ciopts(lc(blue)) lc(blue) ci((results_iv[,2] results_iv[,3])) offset(.05)) ///
	     (matrix(results_ols[,1]), ms(S) mc(red) ciopts(lc(red)) lstyle(p1) lc(red) ci((results_ols[,2] results_ols[,3]))  offset(-.05)) , ///
	xtitle(Bias value)  ytitle(Covariate)  xline(0) leg(off) ylabel(,format(%9.1f))    graphregion(color(white)) plotregion(color(white)) grid(none) ///
	coeflabels(results_z_northing="Northing"
	results_z_easting="Easting"
	results_z_tdi_log="Log Townsend Deprivation Index"
	results_z_age="Age"
	results_male="Male"
	results_breastfed="Breastfed"
	results_z_birthweight="Birthweight"
	results_z_pc1="1st genetic PC"
	results_z_pc2="2nd genetic PC"
	results_z_pc3="3rd genetic PC"
	results_z_pc4="4th genetic PC"
	results_z_pc5="5th genetic PC"
	results_z_pc6="6th genetic PC"
	results_z_pc7="7th genetic PC"
	results_z_pc8="8th genetic PC"
	results_z_pc9="9th genetic PC"
	results_z_pc10="10th genetic PC");
graph export "results/bias_educ.eps",replace fontface("Times New Roman");


	
