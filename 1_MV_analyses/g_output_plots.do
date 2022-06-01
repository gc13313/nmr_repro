
**Read in annotations file for metabolites
import excel "data\biomarker_annotations_updated (1).xlsx", sheet("Sheet1") firstrow case(lower) clear

rename nameintsv outcome
save biomarker_name_info, replace

**append MVR and PNC res
use "MVR", clear
append using "PNC"

gen type_analyses="MV"

tab exposure, m
tab sex, m
tab type_model, m

tab type_model exposure, m

order exposure outcome sex type_analyses type_model beta se pvalue lci uci N
keep exposure outcome sex type_analyses type_model beta se pvalue lci uci N

**merge with metabolite info
merge m:1 outcome using biomarker_name_info
drop _merge

**save in different formats
save "MV_PNC_results.dta", replace

export excel exposure outcome sex type_analyses type_model beta se pvalue lci uci N ///
using "MV_PNC_results", replace

export delimited using "MV_PNC_results.csv", replace

