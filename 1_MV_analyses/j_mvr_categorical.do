/*******************************************************************************
SENSITIVITY ANALYSIS WITH EACH EXPOSURE CATGEORISED


*******************************************************************************/

/*******************************************************************************
AGE AT MENARCHE
*******************************************************************************/

use reprotraits_NMR_marker_data_ds.dta, clear

**AGE AT MENARCHE: Tertiles seems appropriate here
xtile tert_age_menarche=age_menarche, n(3)
sum age_menarche if tert_age_menarche==1 
sum age_menarche if tert_age_menarche==2
sum age_menarche if tert_age_menarche==3

**make exclusions
*drop if age_menarche<8
*drop if age_menarche>19 & age_menarche!=.

fre tert_age_menarche
*recode tert_age_menarche (2=0) 
*recode tert_age_menarche (3=2)
recode tert_age_menarche (1=0) (2=1) (3=2)

sum age_menarche if tert_age_menarche==0 
sum age_menarche if tert_age_menarche==1
sum age_menarche if tert_age_menarche==2


cd "$results"

fre tert_age_menarche

*use model_unadj_menarche, clear
tempname memhold
	postfile `memhold'  str30 type_model str30 outcome str30 exposure N_ref N_1 N_2  beta1 se1 lci1 uci1 double pvalue1 beta2 se2 lci2 uci2 double pvalue2  using model2_main_menarche, replace
	
	foreach x of varlist tert_age_menarche {
	foreach y of varlist Total_C - S_HDL_TG_pct {
    regress `y' i.`x' i.edu i.bodysize_age10 c.age

		matrix M = r(table)
		local beta1 = M[1,2]
		local se1 = M[2,2]
		local lci1 = M[5,2]
		local uci1 = M[6,2]
		local pvalue1 = M[4,2]
		
		local beta2 = M[1,3]
		local se2 = M[2,3]
		local lci2 = M[5,3]
		local uci2 = M[6,3]
		local pvalue2 = M[4,3]
	
	
			tab tert_age_menarche if tert_age_menarche==0 &  `y'!=. & edu!=. & bodysize_age10!=. & age!=.
		local N_ref=r(N)
	tab tert_age_menarche if tert_age_menarche==1 &  `y'!=. & edu!=. & bodysize_age10!=. & age!=.
		local N_1=r(N)
	tab tert_age_menarche if tert_age_menarche==2 &  `y'!=. & edu!=. & bodysize_age10!=. & age!=.
		local N_2=r(N)

	
		post `memhold' ("2") ("`y'") ("`x'") (`N_ref') (`N_1') (`N_2') (`beta1') (`se1') (`lci1') (`uci1') (`pvalue1')  (`beta2') (`se2') (`lci2') (`uci2') (`pvalue2')
	
}
	}
postclose `memhold'


/*******************************************************************************
PARITY
*******************************************************************************/
use reprotraits_NMR_marker_data_ds.dta, clear

**AGE AT MENOPAUSE: Quartiles seems appropriate here
sum parity
sum n_2734_0_0

fre n_2734_0_0

gen parity_cont=n_2734_0_0
recode parity_cont -3=.

sum parity_cont
gen cat_parity=0 if parity_cont==0
replace cat_parity=1 if parity_cont==1
replace cat_parity=2 if parity_cont==2
replace cat_parity=3 if parity_cont>=3 & parity_cont!=.


sum parity_cont if cat_parity==0 
sum parity_cont if cat_parity==1
sum parity_cont if cat_parity==2
sum parity_cont if cat_parity==3

cd "$results"

tempname memhold
	postfile `memhold'  str30 type_model str30 outcome str30 exposure N_ref N_1 N_2 N_3  beta1 se1 lci1 uci1 double pvalue1 beta2 se2 lci2 uci2 double pvalue2  beta3 se3 lci3 uci3 double pvalue3 using model2_main_parity, replace
	
	foreach x of varlist cat_parity {
	foreach y of varlist Total_C - S_HDL_TG_pct {
    regress `y' i.`x' i.edu i.bodysize_age10 c.age

		matrix M = r(table)
		local beta1 = M[1,2]
		local se1 = M[2,2]
		local lci1 = M[5,2]
		local uci1 = M[6,2]
		local pvalue1 = M[4,2]
		
		local beta2 = M[1,3]
		local se2 = M[2,3]
		local lci2 = M[5,3]
		local uci2 = M[6,3]
		local pvalue2 = M[4,3]
	
		local beta3 = M[1,4]
		local se3 = M[2,4]
		local lci3 = M[5,4]
		local uci3 = M[6,4]
		local pvalue3 = M[4,4]
	
			tab cat_parity if cat_parity==0 &  `y'!=. & edu!=. & bodysize_age10!=. & age!=.
		local N_ref=r(N)
	tab cat_parity if cat_parity==1 &  `y'!=. & edu!=. & bodysize_age10!=. & age!=.
		local N_1=r(N)
	tab cat_parity if cat_parity==2 &  `y'!=. & edu!=. & bodysize_age10!=. & age!=.
		local N_2=r(N)
	tab cat_parity if cat_parity==3 &  `y'!=. & edu!=. & bodysize_age10!=. & age!=.
		local N_3=r(N)

	
		post `memhold' ("2") ("`y'") ("`x'") (`N_ref') (`N_1') (`N_2') (`N_3') (`beta1') (`se1') (`lci1') (`uci1') (`pvalue1')  (`beta2') (`se2') (`lci2') (`uci2') (`pvalue2')   (`beta3') (`se3') (`lci3') (`uci3') (`pvalue3') 
	
}
	}
postclose `memhold'
/*******************************************************************************
AGE AT MENOPAUSE
*******************************************************************************/

use reprotraits_NMR_marker_data_ds.dta, clear

**AGE AT MENOPAUSE: Quartiles seems appropriate here
xtile quart_age_menopause=age_menopause, n(4)
sum age_menopause if quart_age_menopause==1 
sum age_menopause if quart_age_menopause==2
**# Bookmark #1
sum age_menopause if quart_age_menopause==3
sum age_menopause if quart_age_menopause==4

**make exclusions
*drop if age_menopause<30
*drop if age_menopause>60 & age_menopause!=.

fre quart_age_menopause
recode quart_age_menopause (1=0) (2=1) (3=2) (4=3)

sum age_menopause if quart_age_menopause==0 
sum age_menopause if quart_age_menopause==1
sum age_menopause if quart_age_menopause==2
sum age_menopause if quart_age_menopause==3

cd "$results"

fre quart_age_menopause

tempname memhold
	postfile `memhold'  str30 type_model str30 outcome str30 exposure N_ref N_1 N_2 N_3  beta1 se1 lci1 uci1 double pvalue1 beta2 se2 lci2 uci2 double pvalue2  beta3 se3 lci3 uci3 double pvalue3 using model2_main_menopause, replace
	
	foreach x of varlist quart_age_menopause {
	foreach y of varlist Total_C - S_HDL_TG_pct {
    regress `y' i.`x' i.edu i.bodysize_age10 c.age

		matrix M = r(table)
		local beta1 = M[1,2]
		local se1 = M[2,2]
		local lci1 = M[5,2]
		local uci1 = M[6,2]
		local pvalue1 = M[4,2]
		
		local beta2 = M[1,3]
		local se2 = M[2,3]
		local lci2 = M[5,3]
		local uci2 = M[6,3]
		local pvalue2 = M[4,3]
	
		local beta3 = M[1,4]
		local se3 = M[2,4]
		local lci3 = M[5,4]
		local uci3 = M[6,4]
		local pvalue3 = M[4,4]
	
			tab quart_age_menopause if quart_age_menopause==0 &  `y'!=. & edu!=. & bodysize_age10!=. & age!=.
		local N_ref=r(N)
	tab quart_age_menopause if quart_age_menopause==1 &  `y'!=. & edu!=. & bodysize_age10!=. & age!=.
		local N_1=r(N)
	tab quart_age_menopause if quart_age_menopause==2 &  `y'!=. & edu!=. & bodysize_age10!=. & age!=.
		local N_2=r(N)
	tab quart_age_menopause if quart_age_menopause==3 &  `y'!=. & edu!=. & bodysize_age10!=. & age!=.
		local N_3=r(N)

	
		post `memhold' ("2") ("`y'") ("`x'") (`N_ref') (`N_1') (`N_2') (`N_3') (`beta1') (`se1') (`lci1') (`uci1') (`pvalue1')  (`beta2') (`se2') (`lci2') (`uci2') (`pvalue2')   (`beta3') (`se3') (`lci3') (`uci3') (`pvalue3') 
	
}
	}
postclose `memhold'
