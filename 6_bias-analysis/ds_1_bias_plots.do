//Neil Davies 15/03/17
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

use "/Volumes/Height_BMI_and_schooling/MR Myopia education/workingdata/data.dta",clear

rename eid n_eid
joinby n_eid using "/Volumes/Height_BMI_and_schooling/MR Myopia education/rawdata/exclusions.dta",unmatched(both)

drop if _m==3

//Normalize continious variables

foreach i in northing easting tdi_log age birthweight ave_mse_clean eduyears_clean pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10{
	egen z_`i'=std(`i')
	}

foreach i in z_northing z_easting z_tdi_log z_age sex_genetic breastfed z_birthweight z_pc1 z_pc2 z_pc3 z_pc4 z_pc5 z_pc6 z_pc7 z_pc8 z_pc9 z_pc10{
	hetero_test "results/myopia_bias_plots_`i'" `i' ave_mse_clean myopia_as_dosage 
	}

foreach i in z_northing z_easting z_tdi_log z_age sex_genetic breastfed z_birthweight z_pc1 z_pc2 z_pc3 z_pc4 z_pc5 z_pc6 z_pc7 z_pc8 z_pc9 z_pc10{
	hetero_test "results/educ_bias_plots_`i'" `i' eduyears_clean ea_as_dosage 
	}

use "results/myopia_bias_plots_z_tdi_log",clear
foreach i in z_northing z_easting z_tdi_log z_age sex_genetic breastfed z_birthweight z_pc1 z_pc2 z_pc3 z_pc4 z_pc5 z_pc6 z_pc7 z_pc8 z_pc9 z_pc10{
	append using "results/myopia_bias_plots_`i'"
	}
foreach i in z_northing z_easting z_tdi_log z_age sex_genetic breastfed z_birthweight z_pc1 z_pc2 z_pc3 z_pc4 z_pc5 z_pc6 z_pc7 z_pc8 z_pc9 z_pc10{
	append using "results/educ_bias_plots_`i'"
	}
	
//Plot the results:
mkmat coef ci_lower ci_upper if strrpos(var,"mse")!=0 & strrpos(var,"xb2")!=0, mat(results_iv) rownames(_estimates_name )
mkmat coef ci_lower ci_upper if strrpos(var,"mse")!=0 & strrpos(var,"xb1")!=0, mat(results_ols) rownames(_estimates_name )
	
coefplot (matrix(results_iv[,1]), ms(T) mc(blue) lc(blue) ciopts(lc(blue)) lc(blue) ci((results_iv[,2] results_iv[,3])) offset(.05)) ///
	     (matrix(results_ols[,1]), ms(S) mc(red) ciopts(lc(red)) lstyle(p1) lc(red) ci((results_ols[,2] results_ols[,3]))  offset(-.05)) , ///
	xtitle(Bias value)  ytitle(Covariate)  xline(0) leg(off) ylabel(,format(%9.1f))  plotregion(lc(white))
graph export "results/bias_mse.eps",replace fontface("Times New Roman")

mkmat coef ci_lower ci_upper if strrpos(var,"edu")!=0 & strrpos(var,"xb2")!=0, mat(results_iv) rownames(_estimates_name )
mkmat coef ci_lower ci_upper if strrpos(var,"edu")!=0 & strrpos(var,"xb1")!=0, mat(results_ols) rownames(_estimates_name )
	
coefplot (matrix(results_iv[,1]), ms(T) mc(blue) lc(blue) ciopts(lc(blue)) lc(blue) ci((results_iv[,2] results_iv[,3])) offset(.05)) ///
	     (matrix(results_ols[,1]), ms(S) mc(red) ciopts(lc(red)) lstyle(p1) lc(red) ci((results_ols[,2] results_ols[,3]))  offset(-.05)) , ///
	xtitle(Bias value)  ytitle(Covariate)  xline(0) leg(off) ylabel(,format(%9.1f))  plotregion(lc(white))
graph export "results/bias_educ.eps",replace fontface("Times New Roman")


	
