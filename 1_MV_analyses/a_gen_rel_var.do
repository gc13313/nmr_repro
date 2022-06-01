use "data.37048_gc2", clear

gen FID = n_eid
lab var FID "FID"
gen IID = FID
lab var IID "IID"

/*******************************************************************************
1. Generate confounders/covariables
*******************************************************************************/

** Sex
fre n_31_0_0
clonevar sex = n_31_0_0
fre sex

** Age
clonevar age = n_21003_0_0
egen agegr = cut(age), at(35,40,45,50,55,60,65,75) label
recode agegr (0=1) // only one aged <40, so included in the 40- category
lab var agegr "Age in 5-y groups"

tab n_2724_0_0 agegr, m


** Ethnicity
clonevar ethnic = n_21000_0_0
replace ethnic = 1 if ethnic==1001 | ethnic==1002 | ethnic==1003 // White
replace ethnic = 2 if ethnic==2001 | ethnic==2002 | ethnic==2003 | ethnic==2004 // Mixed
replace ethnic = 3 if ethnic==3001 | ethnic==3002 | ethnic==3003 | ethnic==3004 // Asian or Asian British
replace ethnic = 4 if ethnic==4001 | ethnic==4002 | ethnic==4003  // Black or Black British
lab def ethnic 1 "White" 2 "Mixed" 3 "Asian" 4 "Black" 5 "Chinese" 6 "Others", replace
lab val ethnic ethnic

** Household income
clonevar income = n_738_0_0
replace income = n_738_1_0 if income==. | income<0
replace income = n_738_2_0 if income==. | income<0
lab def income 1"<18,000" 2"18,000-30,999" 3"31,000-51,999" 4"52,000-100,000" 5">100,000"
lab val income income

** Education
clonevar education = n_6138_0_0
replace education = n_6138_0_1 if education==. | education==-3
replace education = n_6138_0_2 if education==. | education==-3
replace education = n_6138_0_3 if education==. | education==-3
replace education = n_6138_0_4 if education==. | education==-3
replace education = n_6138_0_5 if education==. | education==-3
recode education (1=4) (2=3) (3=2) (4=2) (5=2) (6=2) (-7=1)
lab def education 4"College/ university" 3"A level" 2"O level or CSEs or other" 1"None of the above"
lab val education education

** Deprivation index
clonevar dep_index = n_189_0_0
xtile depindex_5 = dep_index, nq(5)
lab var depindex_5 "Quintiles deprivation index"

** Smoking status
clonevar smoking = n_20116_0_0

replace smoking = smoking + 1
lab def smoking 1"Never" 2"Previous" 3"Current", modify
lab val smoking smoking


** Alcohol consumption frequency
clonevar alcohol = n_1558_0_0

** BMI (anthropometry)
clonevar bmi = n_21001_0_0

** BMI (BIA)
clonevar bmi_bia = n_23104_0_0

** Weight (kg)
clonevar weight = n_21002_0_0

** Height (standing height)
clonevar height = n_50_0_0

** Waist circumference
clonevar wc = n_48_0_0

** Month of birth
clonevar month_birth = n_52_0_0

** Month attending assessment centre
clonevar month_att = n_55_0_0

** Skin colour
clonevar skin_colour = n_1717_0_0

** Easy of skin tanning
clonevar skin_tanning = n_1727_0_0

foreach x of varlist ethnic income smoking alcohol skin_colour skin_tanning {
	replace `x'=. if `x'<0
}

/*******************************************************************************
2. Generate maternal variables
*******************************************************************************/


** Ever pregnant 
gen ever_pregnant = 0
replace ever_pregnant = 1 if (n_2734_0_0>0 & n_2734_0_0!=.) | n_2774_0_0==1 // defined as having at least one live birth or ever having pregnancy complication
ta ever_pregnant // 224,381 (84.8%) ever had a pregnancy
*keep if ever_pregnant==1


** Birth weight
gen bw = n_20022_0_0
lab var bw "Birth weight (kg)"

** Currently pregnant
clonevar pregnant = n_3140_0_0
recode pregnant (1=2) (0=1) (2=.) // unsure recoded to missing
lab val pregnant case_control

** Number of live births
clonevar nbirths = n_2734_0_0
clonevar parity = nbirths
fre parity
recode parity -3=.
gen livebirth_cont=parity
recode parity (4/max=4)
fre parity

** Age at first and last live birth
clonevar age_firstbirth = n_2754_0_0
replace age_firstbirth=n_2754_1_0 if  (n_2754_0_0==. | n_2754_0_0<0) & (n_2754_1_0!=. | n_2754_1_0>0)
replace age_firstbirth=n_2754_2_0 if  (n_2754_0_0==. | n_2754_0_0<0) & (n_2754_2_0!=. | n_2754_2_0>0)
replace age_firstbirth=n_2754_3_0 if  (n_2754_0_0==. | n_2754_0_0<0) & (n_2754_3_0!=. | n_2754_3_0>0)
replace age_firstbirth =. if age_firstbirth<0 // did not have live birth recoded as missing

count if age_firstbirth==. & ever_pregnant==1
lab val age_firstbirth age_birth
clonevar age_lastbirth = n_2764_0_0
replace age_lastbirt=. if age_lastbirth<0 // did not have live birth recoded as missing
lab val age_lastbirth age_birth


**Age at menarche
fre n_2714_0_0
clonevar age_menarche = n_2714_0_0
replace age_menarche=n_2714_1_0 if  (n_2714_0_0==. | n_2714_0_0<0) & (n_2714_1_0!=. | n_2714_1_0>0)
replace age_menarche=n_2714_2_0 if  (n_2714_0_0==. | n_2714_0_0<0) & (n_2714_2_0!=. | n_2714_2_0>0)
replace age_menarche=n_2714_3_0 if  (n_2714_0_0==. | n_2714_0_0<0) & (n_2714_3_0!=. | n_2714_3_0>0)
replace age_menarche=. if age_menarche<0 
fre age_menarche
hist age_menarche

**Age at menopause
fre n_3581_0_0
clonevar age_menopause = n_3581_0_0
replace age_menopause=n_3581_1_0 if  (n_3581_0_0==. | n_3581_0_0<0) & (n_3581_1_0!=. | n_3581_1_0>0)
replace age_menopause=n_3581_2_0 if  (n_3581_0_0==. | n_3581_0_0<0) & (n_3581_2_0!=. | n_3581_2_0>0)
replace age_menopause=n_3581_3_0 if  (n_3581_0_0==. | n_3581_0_0<0) & (n_3581_3_0!=. | n_3581_3_0>0)
replace age_menopause=. if age_menopause<0 
fre age_menopause
hist age_menopause

*remove those who had a surgical menopause (not needed)
fre n_2824_0_0
fre n_2724_0_0
*replace age_menopause=. if (n_2824_0_0!=. | n_2824_1_0!=. | n_2824_2_0!=. | n_2824_3_0!=.) & n_2724_0_0==2


sum age_menopause
sum age_menopause if n_2724_0_0==2
sum age_menopause if n_2724_0_0==3
sum age_menopause if n_2724_0_0==1
sum age_menopause if n_2724_0_0==4

sum age_menopause if n_2824_0_0!=.

list age_menopause n_2824_0_0 if n_2824_0_0!=. & age_menopause!=.
count if n_2824_0_0!=. & age_menopause!=.


**WOMEN parity
sum livebirth_cont

****MEN - noofchil
fre n_2405_0_0 
fre n_2405_1_0 
fre n_2405_2_0 
fre n_2405_3_0

clonevar noofchil=n_2405_0_0 
fre noofchil
recode noofchil -3=.
recode noofchil -1=.

replace noofchil=n_2405_1_0 if (n_2405_1_0!=. & n_2405_1_0!=-1 & n_2405_1_!=-3)
replace noofchil=n_2405_2_0 if (n_2405_2_0!=. & n_2405_2_0!=-1 & n_2405_2_!=-3)
replace noofchil=n_2405_3_0 if (n_2405_3_0!=. & n_2405_3_0!=-1 & n_2405_3_!=-3)


fre noofchil
gen noofchil_cont=noofchil
gen noofchil_cat=noofchil
recode noofchil_cat (4/max=4)
fre noofchil

fre noofchil if n_2405_1_0>=0 &  n_2405_1_0!=.

*comparative body shape
fre n_1687_0_0
clonevar bodysize_age10 = n_1687_0_0
replace bodysize_age10=. if n_1687_0_0==-3 | n_1687_0_0==-1
replace bodysize_age10=n_1687_1_0 if bodysize_age10==. & (n_1687_1_0!=. | n_1687_1_0>0)
replace bodysize_age10=n_1687_2_0 if bodysize_age10==. & (n_1687_2_0!=. | n_1687_2_0>0)

sum nbirths

recode nbirths -3=.
recode nbirths -1=.

fre education
gen edu=education
recode edu (1=0) (2=1) (3=2) (4=3)

recode bodysize_age10 (3=0) (-1=.)

save reprotraits_ds.dta, replace

