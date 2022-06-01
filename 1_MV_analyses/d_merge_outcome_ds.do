cd ""


********************************************************************************
********************        2. Merge with IEU linker IDS       *****************
********************************************************************************
import delimited "linker.csv", clear
rename app IID
tempfile ieulinker
save `ieulinker'
clear

use reprotraits_ds.dta, clear

*** Remove withdrawns
merge 1:1 IID using `ieulinker'
tab _merge
drop if _merge==2 | _merge==1
drop _merge
*drop eid

save reprotraits_data_ieu_ds.dta, replace
use reprotraits_data_ieu_ds.dta, clear

********************************************************************************
********************           4. Merge with NMR data          *****************
********************************************************************************
import delimited "nmr_dat_all.txt", delimiter(space) varnames(1) case(preserve) clear 
rename IID ieu
merge 1:1 ieu using reprotraits_data_ieu_ds.dta, force
tab _merge
drop if _merge==2 
drop _merge
save reprotraits_nmr_data_ds.dta, replace





********************************************************************************
********************          4. Merge with MARKER data        *****************
********************************************************************************
import delimited "markers_all.txt", delimiter(space) varnames(1) case(preserve) clear 
rename IID ieu
merge 1:1 ieu using reprotraits_data_ieu_ds.dta, force
tab _merge
drop if _merge==2 
drop _merge
save reprotraits_marker_data_ds.dta, replace


********************************************************************************
********************    5. Merge with NMR with MARKER data     *****************
********************************************************************************
use reprotraits_nmr_data_ds.dta, clear
merge 1:1 ieu using reprotraits_marker_data_ds.dta, force
drop _merge


rename (*_int) (*)

** making each metabolite numeric
foreach x of varlist Total_C - S_HDL_TG_pct ALB APOA APOB CHOL GLU HDL LDLD TRIG {
replace `x'="." if `x'=="NA"
destring `x', replace
}


save reprotraits_NMR_marker_data_ds.dta, replace















