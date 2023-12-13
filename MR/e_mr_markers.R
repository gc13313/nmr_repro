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

#install.packages(c("data.table", "purrr", "devtools", "xlsx", "readxl", "dplyr", "ggplot2", "remotes"))
#devtools::install_github("MRCIEU/TwoSampleMR")
#devtools::install_github("MRCIEU/MRInstruments")
#devtools::install_github("NightingaleHealth/ggforestplot")


##
## Load libraries

library(data.table)
library(purrr)
library(TwoSampleMR)
#library(MRInstruments)
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
## Markers info file

markers.info <- read_excel("data/markers_info.xlsx", sheet = 1)


############################################################################
#                       Prepare and harmonise data                         #
############################################################################

##
## Exposure data

expdat <- read.table("data/expdat.txt") %>% filter(analysis == 1)

aamdat <- filter(expdat, abbrev=="AAM" & ukb_included=="yes")


##
## Covariable data (BMI - GIANT females)

filter(ao, id == "ieu-a-974")

bmidat <- extract_outcome_data(snps = aamdat$SNP, outcomes = 'ieu-a-974')


##
## Outcome data

outdat <- read.table("data/markers/markers500_female_sumdat_ivs.txt") %>%
    			rename_all(., .funs = list(~ tolower(.))) %>%
    			rename(effect_allele = allele1, other_allele = allele0, eaf = a1freq, pval = p_bolt_lmm_inf) %>%
				filter(info >= 0.8) %>%
    			format_data(.,
    							type = "outcome",
    							phenotype_col = "id",
    							snp_col = "snp",
    							pos_col = "bp"
    						)
 
length(unique(outdat$outcome)) # Observed number of outcomes (expected: 8 markers)


##
## Harmonise exposure and outcome data

dat <- harmonise_data(expdat, outdat)


############################################################################
#                                  Run MR                                  #
############################################################################

##
## Run MR

res.raw <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))			

 
##
## Format MR results

res <-  res.raw %>%
			merge(., exp.info, by.x = "id.exposure", by.y = "id") %>%
			mutate(outcome = gsub("markers500_female_|_500_female_imputed.txt.gz_TMPDAT.txt", "", outcome)) %>%
			mutate(
					method_abbrev = case_when(method == "Inverse variance weighted" ~ "ivw",
											  method == "MR Egger" ~ "mregger",
											  method == "Weighted median" ~ "wmedian"
				   ),
				   method = factor(method, levels=c("MR Egger", "Weighted median", "Inverse variance weighted")),
				   type_analyses = "MR"
				  ) %>%
			merge(., markers.info, by.x = "outcome", by.y = "var_name")
			

##
## Save MR results

mr.res <- res %>%
 			mutate(outcome = nmr_name) %>%
			arrange(exposure, outcome, method) %>%
        	select(id.exposure, id.outcome, trait, abbrev, outcome, biomarker_name, nsnp, method, b, se, pval) %>%
			mutate(assay = "Blood biochemistry")
        	
writexl::write_xlsx(list(mr_res_markers500 = mr.res), path = "./results/mr_results_markers500.xlsx")


############################################################################
#                                  Run MR                                  #
############################################################################

##
## Import MR results on NMR measures

all.res <- read_excel("results/mr_results.xlsx", sheet = 1) %>%
			filter(outcome %in% markers.info$nmr_name) %>%
			filter(priority == 1) %>%
			mutate(assay = "NMR metabolomics") %>%
			bind_rows(mr.res)


############################################################################
#                          Run Multivariable MR                            #
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
dat3 <- merge(dat1, dat2, by = "SNP")


##
## Function for performing MVMR

run_mvmr <- function(out) {

					# Filter data for one outcome
					df <- filter(dat3, outcome == out)
					
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
				mutate(
						outcome = gsub("markers500_female_|_500_female_imputed.txt.gz_TMPDAT.txt", "", outcome),
						analysis = "MR (adjusted for BMI)"
						) %>%
      			rename(b = Estimate, se = "Std. Error", pval = "Pr(>|t|)") %>%
      			filter(exposure == "Age at menarche") 

			
mvmr_res <- mvmr_res %>%
				mutate(nmr_name = case_when(
											outcome == "ALB" ~ "Albumin",
											outcome == "APOA" ~ "ApoA1", 
											outcome == "APOB" ~ "ApoB", 
											outcome == "CHOL" ~ "Total_C",
											outcome == "GLU" ~ "Glucose", 
											outcome == "HDL" ~ "HDL_C", 
											outcome == "LDLD" ~ "Clinical_LDL_C", 
											outcome == "TRIG" ~ "Total_TG"
											))


##
## Select univariable MR results (IVW) for AAM

uvmr_res <- mutate(res, analysis = "MR (unadjusted)") %>%
			filter(method == "Inverse variance weighted" & abbrev=="AAM" & ukb_included=="yes")

 
##
## Format MR results

p2 <- bind_rows(uvmr_res, mvmr_res) %>% 
				arrange(nmr_name) %>%
				ggforestplot::forestplot(
							  df = .,
							  name = nmr_name,
							  estimate = b,
							  se = se,
							  pvalue = pval,
							  psignif = 0.05/(18*3),
							  xlab = "1-SD increment in metabolite \nper 1 year increase in age at menarche",
							  colour = analysis
							) +
					theme(						
							legend.text = element_text(size = 20),
							legend.title = element_blank(),
							axis.text.x = element_text(size=20),
							axis.title.x = element_text(size=20),
							axis.text.y = element_text(size=20),
							strip.text = element_text(size = 20)
							) +
					scale_color_manual(values=c("blue2", "black")) 

ggexport(p2, filename = "./results/plots/forest_AAM_nmr_MVMR_markers500.png", width = 1200, height = 800)


q("no")					