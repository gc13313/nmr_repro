/*******************************************************************************
SCRIPT OF RESTRICTED CUBIC SPLINES TO ASSESS NON-LINEARITY


*******************************************************************************/

/*******************************************************************************
AGE AT MENARCHE
*******************************************************************************/


sum age_menarche, d


_pctile age_menarche, percentile(5 27.5 50 72.5 95)
gen p1 = r(r1)
gen p2 = r(r2)
gen p3 = r(r3)
gen p4 = r(r4)
gen p5 = r(r5)

list p1 p2 p3 p4 p5 if _n==1

sum age_menarche, d

gen age_menarche_cat=0 if age_menarche<9
replace age_menarche_cat=1 if age_menarche>8 & age_menarche<11
replace age_menarche_cat=2 if age_menarche>10 & age_menarche<13
replace age_menarche_cat=3 if age_menarche>12 & age_menarche<15
replace age_menarche_cat=4 if age_menarche>14 & age_menarche<17
replace age_menarche_cat=5 if age_menarche>16 & age_menarche<19
replace age_menarche_cat=6 if age_menarche>18 & age_menarche!=.

 label define age_menarche_cat ///
 0 "<9" ///
 1 "9-10" ///
 2 "11-12" ///
 3 "13-14" ///
 4 "15-16" ///
 5 "17-18" ///
 6 ">18", replace
label values age_menarche_cat age_menarche_cat

fre age_menarche_cat if sex=="Female"


mkspline agemenarche_spline = age_menarche, cubic nknots(3) displayknots 

recode education (4=0)

** FOR THE GRAPHS
gen age_menarche_gr=age_menarche if age_menarche>7 & age_menarche<19

*********************************************************************************
* LOOP FOR OUTCOMES
*********************************************************************************
 
 tempname memhold
	postfile `memhold'  str30 outcome rmse_sp df_sp AIC_sp BIC_sp  rmse_l df_l AIC_l BIC_l   using spline_menarche_3, replace
	
 *in a loop
 foreach x of varlist Albumin	ApoA1	Cholines DHA Gln Gly GlycA HDL_CE  HDL_FC  HDL_L HDL_PL LA ///
	LDL_FC 	PUFA Phosphatidylc Phosphoglyc  Sphingomyelins Total_CE Total_FC Total_PL  Unsaturation {
		
 **SPLINE MODEL
 regress `x' agemenarche_spline* i.education i.bodysize_age10 c.age
estimates store `x'_sp

local rmse_sp=e(rmse)

predictnl `x'_mean_sp = _b[_cons]+_b[agemenarche_spline1]*agemenarche_spline1+_b[agemenarche_spline2]*agemenarche_spline2  + _b[age]*60 ///
, se(`x'_se_sp)


generate `x'_lci_sp = `x'_mean_sp - 1.96*`x'_se_sp
generate `x'_uci_sp = `x'_mean_sp + 1.96*`x'_se_sp
 estat ic
 
 
		matrix M = r(S)
		matrix list M
		
		local df_sp = M[1,4]
		local AIC_sp = M[1,5]
		local BIC_sp = M[1,6]
		

**LINEAR MODEL
 regress `x' age_menarche i.education i.bodysize_age10 c.age
estimates store `x'_l

local rmse_l=e(rmse)

predictnl `x'_mean_l = _b[_cons]+_b[age_menarche]*age_menarche  + _b[age]*60 ///
, se(`x'_se_l)

generate `x'_lci_l = `x'_mean_l - 1.96*`x'_se_l
generate `x'_uci_l = `x'_mean_l + 1.96*`x'_se_l
 estat ic

		matrix M = r(S)
		matrix list M
		
		local df_l = M[1,4]
		local AIC_l = M[1,5]
		local BIC_l = M[1,6]
		

**GRAPH
twoway ///
 (line `x'_mean_sp age_menarche_gr, sort lcolor(ebblue))  ///
(rarea `x'_lci_sp `x'_uci_sp age_menarche_gr, color("ebblue%20") sort ) ///
 (line `x'_mean_l age_menarche_gr, sort lcolor(cranberry))  ///
(rarea `x'_lci_l `x'_uci_l age_menarche_gr, color("cranberry%10") sort ) ///
,  ytitle("`x'") xtitle("Age at menarche (years)") xlab(8 (2) 18, grid) ylab(, grid) ///
 legend(off) ///
 name(`x'_menarche, replace) saving(`x'_menarche, replace)
	
 	
		post `memhold' ("`x'") (`rmse_sp') (`df_sp') (`AIC_sp')  (`BIC_sp') (`rmse_l') (`df_l') (`AIC_l')  (`BIC_l') 
	
}
	
postclose `memhold'

/*******************************************************************************
AGE AT MENOPAUSE
*******************************************************************************/

use reprotraits_NMR_marker_data_ds.dta, clear


sum age_menopause, d


_pctile age_menopause, percentile(5 27.5 50 72.5 95)
gen p1 = r(r1)
gen p2 = r(r2)
gen p3 = r(r3)
gen p4 = r(r4)
gen p5 = r(r5)

list p1 p2 p3 p4 p5 if _n==1


sum age_menopause, d

gen age_menopause_cat=1 if age_menopause!=. & age_menopause<40
replace age_menopause_cat=2 if age_menopause>=40 & age_menopause<45
replace age_menopause_cat=3 if age_menopause>=45 & age_menopause<50
replace age_menopause_cat=4 if age_menopause>=50 & age_menopause<52
replace age_menopause_cat=5 if age_menopause>=52 & age_menopause<55
replace age_menopause_cat=6 if age_menopause>=55 & age_menopause!=.

label define age_menopause_cat 1 "<40" 2 "40-44" 3 "45-49" 4 "50-51" 5 "52-54" 6 "55+", replace
label values age_menopause_cat age_menopause_cat

fre age_menopause_cat if age_menopause!=.

mkspline agemeno_spline = age_menopause, cubic nknots(4) displayknots 

list agemeno_spline* age_menopause in 1/20

recode education (4=0)

** REMOVE EXTREMES ON GRAPH
gen age_menopause_gr=age_menopause if age_menopause>29 & age_menopause<66

 
*********************************************************************************
* LOOP FOR OUTCOMES
*********************************************************************************

 
 tempname memhold
	postfile `memhold'  str30 outcome rmse_sp df_sp AIC_sp BIC_sp  rmse_l df_l AIC_l BIC_l   using spline_menopause_4, replace
	
 *in a loop
 foreach x of varlist Albumin	ApoA1	Cholines DHA Gln Gly GlycA HDL_CE  HDL_FC  HDL_L HDL_PL LA ///
	LDL_FC 	PUFA Phosphatidylc Phosphoglyc  Sphingomyelins Total_CE Total_FC Total_PL  Unsaturation ///
	Citrate LDL_L LDL_PL Total_TG HDL_TG VLDL_TG {
		
 **SPLINE MODEL
 regress `x' agemeno_spline* i.education i.bodysize_age10 c.age
estimates store `x'_sp

local rmse_sp=e(rmse)

predictnl `x'_mean_sp = _b[_cons]+_b[agemeno_spline1]*agemeno_spline1+_b[agemeno_spline2]*agemeno_spline2 +_b[agemeno_spline3]*agemeno_spline3 + _b[age]*60 ///
, se(`x'_se_sp)


generate `x'_lci_sp = `x'_mean_sp - 1.96*`x'_se_sp
generate `x'_uci_sp = `x'_mean_sp + 1.96*`x'_se_sp
 estat ic
 
 
		matrix M = r(S)
		matrix list M
		
		local df_sp = M[1,4]
		local AIC_sp = M[1,5]
		local BIC_sp = M[1,6]
		

**LINEAR MODEL
 regress `x' age_menopause i.education i.bodysize_age10 c.age
estimates store `x'_l

local rmse_l=e(rmse)

predictnl `x'_mean_l = _b[_cons]+_b[age_menopause]*age_menopause + _b[age]*60 ///
, se(`x'_se_l)

generate `x'_lci_l = `x'_mean_l - 1.96*`x'_se_l
generate `x'_uci_l = `x'_mean_l + 1.96*`x'_se_l
 estat ic

		matrix M = r(S)
		matrix list M
		
		local df_l = M[1,4]
		local AIC_l = M[1,5]
		local BIC_l = M[1,6]
		

**GRAPH
twoway ///
 (line `x'_mean_sp age_menopause_gr, sort lcolor(ebblue))  ///
(rarea `x'_lci_sp `x'_uci_sp age_menopause_gr, color("ebblue%20") sort ) ///
 (line `x'_mean_l age_menopause_gr, sort lcolor(cranberry))  ///
(rarea `x'_lci_l `x'_uci_l age_menopause_gr, color("cranberry%10") sort ) ///
,  ytitle("`x'") xtitle("Age at menopause (years)") xlab(30 (5) 65, grid) ylab(, grid) ///
 legend(off) ///
 name(`x'_meno_4, replace) saving(`x'_meno_4, replace)
	
 
		post `memhold' ("`x'") (`rmse_sp') (`df_sp') (`AIC_sp')  (`BIC_sp') (`rmse_l') (`df_l') (`AIC_l')  (`BIC_l') 
	
}
	
postclose `memhold'




/*******************************************************************************
PARITY
*******************************************************************************/

use reprotraits_NMR_marker_data_ds.dta, clear

*5th, 27.5th, 50th, 72.5th and 95th for k = 5	

**age at parity and albumin
drop parity
gen parity=n_2734_0_0
recode parity -3=.
sum parity
sum parity, d


_pctile parity, percentile(5 27.5 50 72.5 95)
gen p1 = r(r1)
gen p2 = r(r2)
gen p3 = r(r3)
gen p4 = r(r4)
gen p5 = r(r5)

list p1 p2 p3 p4 p5 if _n==1


sum parity, d

fre parity


mkspline parity_spline = parity, cubic knots(1 2 3) displayknots 


fre edu
fre education
recode education (4=0)

*********************************************************************************
* LOOP FOR OUTCOMES
*********************************************************************************

 
 tempname memhold
	postfile `memhold'  str30 outcome mean_sp lci_sp uci_sp mean_l lci_l uci_l rmse_sp df_sp AIC_sp BIC_sp  rmse_l df_l AIC_l BIC_l   using spline_parity_3, replace
	
 *in a loop
 foreach x of varlist Albumin	ApoA1	Cholines DHA Gln Gly GlycA HDL_CE  HDL_FC  HDL_L HDL_PL LA ///
	LDL_FC 	PUFA Phosphatidylc Phosphoglyc  Sphingomyelins Total_CE Total_FC Total_PL  Unsaturation {
		
 **SPLINE MODEL
 regress `x' parity_spline* i.education i.bodysize_age10 c.age
estimates store `x'_sp

local rmse_sp=e(rmse)

predictnl `x'_mean_sp = _b[_cons]+_b[parity_spline1]*parity_spline1+_b[parity_spline2]*parity_spline2 + _b[age]*60 ///
, se(`x'_se_sp)

local mean_sp=`x'_mean_sp

generate `x'_lci_sp = `x'_mean_sp - 1.96*`x'_se_sp
generate `x'_uci_sp = `x'_mean_sp + 1.96*`x'_se_sp
 estat ic
 
 
		matrix M = r(S)
		matrix list M
		
		local df_sp = M[1,4]
		local AIC_sp = M[1,5]
		local BIC_sp = M[1,6]
		
local lci_sp=`x'_lci_sp
local uci_sp=`x'_uci_sp

**LINEAR MODEL
 regress `x' parity i.education i.bodysize_age10 c.age
estimates store `x'_l

local rmse_l=e(rmse)

predictnl `x'_mean_l = _b[_cons]+_b[parity]*parity + _b[age]*60 ///
, se(`x'_se_l)

local mean_l=`x'_mean_l

generate `x'_lci_l = `x'_mean_l - 1.96*`x'_se_l
generate `x'_uci_l = `x'_mean_l + 1.96*`x'_se_l
 estat ic

		matrix M = r(S)
		matrix list M
		
		local df_l = M[1,4]
		local AIC_l = M[1,5]
		local BIC_l = M[1,6]
		
local lci_l=`x'_lci_l
local uci_l=`x'_uci_l


 	
		post `memhold' ("`x'") (`mean_sp') (`lci_sp') (`uci_sp') (`mean_l') (`lci_l') (`uci_l') (`rmse_sp') (`df_sp') (`AIC_sp')  (`BIC_sp') (`rmse_l') (`df_l') (`AIC_l')  (`BIC_l') 
	
}
	
postclose `memhold'
