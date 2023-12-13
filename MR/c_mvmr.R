############################################################################
#                                Set - UP                                  #
############################################################################

##
## Clear the work environment

rm(list = ls())


##
## Setting repository (UoB)

options(repos = c(CRAN ="http://www.stats.bris.ac.uk/R/"))


##
## Setting digits

options(digits = 10)


##
## Library directory

.libPaths()


##
## Updating packages 

#update.packages(ask = "FALSE")


##
## Install packages

#install.packages(c("data.table", "purrr", "devtools", "xlsx", "readxl", "dplyr", "ggplot2"))
#devtools::install_github("MRCIEU/TwoSampleMR")
#devtools::install_github("MRCIEU/MRInstruments")
#devtools::install_github("NightingaleHealth/ggforestplot")


##
## Load libraries

library(data.table)
library(purrr)
library(TwoSampleMR)
library(MRInstruments)
#library(xlsx)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggforestplot)


##
## List available GWASs in MR-Base catalogue

ao <- available_outcomes()


##
## Exposure info file

exp.info <- read_excel("data/exp_info.xlsx", sheet = "iv_exp")


##
## NMR metabolites info file

nmr.info <- read_excel("data/biomarker_annotations_updated (1).xlsx")
  

############################################################################
#                             Extract data                                 #
############################################################################

##
## Exposure data (AAM: N=329345)

aamdat <- read.table("data/expdat.txt") %>% filter(id.exposure == "ST3_41588_2017_BFng3841_MOESM301_ESM")


##
## Covariable data (BMI - GIANT females)

filter(ao, id == "ieu-a-974")

bmidat <- extract_outcome_data(snps = aamdat$SNP, outcomes = 'ieu-a-974')


##
## Outcome data

outdat <- read.table("data/nmr/nmr_female_sumdat_ivs.txt") %>%
    			rename_all(., .funs = list(~ tolower(.))) %>%
    			rename(effect_allele = allele1, other_allele = allele0, eaf = a1freq, pval = p_bolt_lmm_inf) %>%
				  filter(grepl("_female", id)) %>%
    			format_data(.,
    							type = "outcome",
    							phenotype_col = "id",
    							snp_col = "snp",
    							pos_col = "bp"
    						)
		
# Observed number of outcomes (expected: 249)
length(unique(outdat$outcome))


############################################################################
#                                  Run UVMR                                #
############################################################################

uvmr_res <- harmonise_data(aamdat, outdat) %>% 
            mr(., method_list = "mr_ivw") %>%
            mutate(analysis = "MR (unadjusted)")


############################################################################
#                                  Run MVMR                                #
############################################################################

##
## Harmonise data

# First, orient AAM and BMI data to the same effect allele
dat1 <- harmonise_data(aamdat, bmidat) %>% 
    		rename_at(vars(ends_with("outcome")), ~sub("outcome", "bmi", .)) %>%
    		filter(mr_keep == T)

# Second, orient AAM and outcomes data to the same effect allele
dat2 <- select(aamdat, id.exposure:units.exposure) %>%
    		harmonise_data(., outdat) %>% 
    		select(!ends_with("exposure")) %>%
    		filter(mr_keep == T)

# Third, merge data
dat <- merge(dat1, dat2, by = "SNP")


##
## Function for performing MVMR

run_mvmr <- function(out) {

					# Filter data for one outcome
					df <- filter(dat, outcome == out)
					
					# Format datasets to MVMR package
					df <- MVMR::format_mvmr(BXGs = df[, c("beta.exposure", "beta.bmi")],
											BYG = df[, "beta.outcome"],
											seBXGs = df[, c("se.exposure", "se.bmi")],
											seBYG = df[, "se.outcome"],
											RSID = df[, "SNP"]
											) 
					
					# Test for weak instruments (assume no sample overlap between AAM and BMI datasets)
					weakiv <- MVMR::strength_mvmr(r_input = df, gencov = 0) %>%
      								as.data.table(.) %>%
      								rename(cF_AAM = exposure1, cF_BMI = exposure2) %>%
      								mutate(outcome = out)
								
					# Test for horizontal pleiotropy (assume no sample overlap between AAM and BMI datasets)
					pleio <- MVMR::pleiotropy_mvmr(r_input = df, gencov = 0) %>%
    								as.data.table(.) %>%
    								mutate(outcome = out)
								
					# Run MVMR to estimate direct effects
					eff <- MVMR::ivw_mvmr(r_input = df) %>%
  								as.data.table(.) %>%
  								mutate(., exposure = c("Age at menarche", "BMI")) %>%
  								mutate(outcome = out)
					
					# Combine MVMR estimates
					res <- merge(weakiv, pleio, by = "outcome") %>%
						    	merge(., eff, by = "outcome")

}				  
				  

##
## Run MVMR to estimate direct effects 

mvmr_res <- map(unique(dat$outcome), ~run_mvmr(.)) %>%
      			bind_rows %>%
      			mutate(analysis = "MR (adjusted for BMI)") %>%
      			rename(b = Estimate, se = "Std. Error", pval = "Pr(>|t|)") %>%
      			filter(exposure == "Age at menarche")


############################################################################
#                       Combine UVMR and MVMR results                      #
############################################################################
	
##
## Combine UVMR and MVMR results

res <- bind_rows(uvmr_res, mvmr_res) %>%
    		mutate(outcome = gsub("^nmr_dat_female_|_int_imputed[.]txt[.]gz_TMPDAT[.]txt$", "", outcome)) %>%
    		merge(., nmr.info, by.x = "outcome", by.y = "Name in TSV", all.x = T) %>%
    		rename(Derived = "Derived measure", metab_name = "Biomarker name") %>%
    		mutate(
        				Category = ifelse(Group %in% c("Lipoprotein subclasses", "Relative lipoprotein lipid concentrations"), as.character(Subgroup), as.character(Group)),
    				) %>%
        rename(cF = cF_AAM) %>%
        mutate_at(vars(cF, Qstat), ~round(., 0)) %>%
        mutate_at(vars(b, se, pval, Qpval), ~round(., 3)) %>%
        arrange(exposure, metab_name, analysis) %>%
        select(exposure, outcome, analysis, cF, b, se, pval, Qstat, Qpval, Category, Derived) %>%
        mutate(Category = as.character(map(strsplit(.$Category, split = " [(]"), 1))) %>%
        mutate(Category = case_when(
                Category == "Chylomicrons and extremely large VLDL" ~ "QM & extremely large VLDL",
                Category == "Chylomicrons and extremely large VLDL ratios" ~ "QM & extremely large VLDL ratios",
                Category == "Lipoprotein particle concentrations" ~ "Particle concentration",
                Category == "Lipoprotein particle sizes" ~ "Particle size",
                T ~ as.character(Category)
        )) %>%
        mutate(Category = factor(Category, levels=c(
                                    "Amino acids", "Fatty acids", "Glycolysis related metabolites", "Ketone bodies", "Fluid balance", "Inflammation", "Apolipoproteins",
                                    "Cholesterol", "Cholesteryl esters", "Free cholesterol", "Phospholipids", "Triglycerides", "Other lipids", "Total lipids",
                                    "Particle concentration", "Particle size", "Very large HDL", "Large HDL", "Medium HDL", "Small HDL",
                                    "IDL", "Large LDL", "Medium LDL", "Small LDL",
                                    "QM & extremely large VLDL", "Very large VLDL", "Large VLDL", "Medium VLDL", "Small VLDL", "Very small VLDL",
                                    "QM & extremely large VLDL ratios", "Very large VLDL ratios", "Large VLDL ratios", "Medium VLDL ratios", "Small VLDL ratios", "Very small VLDL ratios", "IDL ratios",
                                    "Large LDL ratios", "Medium LDL ratios", "Small LDL ratios", "Very large HDL ratios", "Large HDL ratios", "Medium HDL ratios", "Small HDL ratios"
        )))



##
## Save MR results

writexl::write_xlsx(list(mvmr_res = mvmr_res), path = "./results/mvmr_results.xlsx")


q("no")		
		
		
		
		
		
		
		
 
