
*********************************************************************************
cd "data\reprotraits"

*nbirths, age_menarche, age_menopause

use reprotraits_NMR_marker_data_ds.dta, clear

fre sex
drop if sex=="Female"
drop nbirths
rename noofchil_cont nbirths

cd "results\reprotraits"

**MODEL 1
tempname memhold
	postfile `memhold'  str30 type_model str30 outcome str30 exposure N  beta se lci uci double pvalue using model1_unadj, replace
	
	foreach x of varlist nbirths {
	foreach y of varlist Total_C - S_HDL_TG_pct {
    regress `y' `x' 

		matrix M = r(table)
		local beta = M[1,1]
		local se = M[2,1]
		local lci = M[5,1]
		local uci = M[6,1]
		local pvalue = M[4,1]
		
		
	
		local N=e(N)

	
		post `memhold' ("1") ("`y'") ("`x'") (`N') (`beta') (`se') (`lci') (`uci') (`pvalue') 
	
}
	}
postclose `memhold'

**MODEL 2
tempname memhold
	postfile `memhold' str30 type_model str30 outcome str30 exposure N  beta se lci uci double pvalue using model2_main, replace
	
	foreach x of varlist nbirths {
	foreach y of varlist Total_C - S_HDL_TG_pct {
    regress `y' `x' i.edu i.bodysize_age10 c.age
		matrix M = r(table)
		local beta = M[1,1]
		local se = M[2,1]
		local lci = M[5,1]
		local uci = M[6,1]
		local pvalue = M[4,1]
		
	
		local N=e(N)

	
		post `memhold'  ("2") ("`y'") ("`x'") (`N') (`beta') (`se') (`lci') (`uci') (`pvalue') 
	
}
	}
postclose `memhold'


**MODEL 3
tempname memhold
	postfile `memhold' str30 type_model str30 outcome str30 exposure N  beta se lci uci double pvalue using model3_overadj, replace
	
	foreach x of varlist nbirths {
	foreach y of varlist Total_C - S_HDL_TG_pct {
		
    regress `y' `x'  i.edu i.bodysize_age10 i.smoking i.alcohol c.bmi c.age
		matrix M = r(table)
		local beta = M[1,1]
		local se = M[2,1]
		local lci = M[5,1]
		local uci = M[6,1]
		local pvalue = M[4,1]
		
	
		local N=e(N)

	
		post `memhold'  ("3") ("`y'") ("`x'") (`N') (`beta') (`se') (`lci') (`uci') (`pvalue') 
	
}
	}
postclose `memhold'



**append these


use  model1_unadj, clear
append using model2_main
append using model3_overadj


gen sex="male"

save "PNC", replace




