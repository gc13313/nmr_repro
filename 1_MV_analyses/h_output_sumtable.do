**MR RES
import excel "output\mr_results.xlsx", sheet("mr_res") firstrow clear

keep if priority==1 
keep abbrev outcome metab_name Category b se pval type_analyses Derived method


save "mr_long", replace


**MV/PNC RES
use "MV_PNC_results.dta", clear

replace type_analyses="NC" if sex=="male"

gen abbrev="NEB" if exposure=="nbirths"
replace abbrev="ANM" if exposure=="age_menopause"
replace abbrev="AAM" if exposure=="age_menarche"

rename beta b
rename pvalue pval
rename biomarkername metab_name
rename derivedmeasure Derived

rename type_model method

keep b se pval type_analyses method outcome exposure abbrev metab_name N Derived

save "mv_nc_long", replace

**MERGE WITH MR WITH MVR/PNC
append using "mr_long"

drop exposure
sort outcome abbrev type_analyses 


order abbrev Category metab_name   outcome type_analyses b se pval
gen lci=b-1.96*se
gen uci=b+1.96*se

replace type_analyses="Multivariable regression" if type_analyses=="MV"
replace type_analyses="Mendelian randomisation" if type_analyses=="MR"
replace type_analyses="Paternal negative control" if type_analyses=="NC"

sort Derived outcome abbrev type_analyses 
replace Category=Category[_n-1] if metab_name==metab_name[_n-1]

gen main_analyses=1 if method=="Inverse variance weighted" | method=="2"
replace main_analyses=2 if main_analyses==.

sort Derived outcome abbrev  main_analyses type_analyses method 

label define main_analyses 1 "Main" 2 "Sensitivity"

label val main_analyses main_analyses

order Category metab_name outcome type_analyses main_analyses method b se lci uci pval Derived
**Age at mennrache
preserve
keep if abbrev=="AAM"
 export excel Category metab_name outcome type_analyses main_analyses method b se lci uci pval Derived using "Suppl summary data results", sheet("Age at menarche") firstrow(var)replace

restore 



**NEB
**Parity/number of children
preserve
keep if abbrev=="NEB"
 export excel Category metab_name outcome type_analyses main_analyses method b se  lci uci pval Derived using "Suppl summary data results", sheet("Parity number of children") firstrow(var) sheetmodify

restore 


**ANM
**Age at menopause
preserve
keep if abbrev=="ANM"
 export excel Category metab_name outcome type_analyses main_analyses method b se  lci uci pval Derived using "Suppl summary data results", sheet("Age at menopause") firstrow(var) sheetmodify
 
 restore 

 