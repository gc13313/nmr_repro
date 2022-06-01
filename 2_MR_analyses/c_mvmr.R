############################################################################
#                                Set - UP                                  #
############################################################################

##
## Set working directory 

setwd(paste0(Sys.getenv('SCRATCH'), "/repro_health/"))


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

#install.packages(c("data.table", "purrr", "devtools", "dplyr"))
#devtools::install_github("MRCIEU/TwoSampleMR")
#devtools::install_github("MRCIEU/MRInstruments")


##
## Load libraries

library(data.table)
library(purrr)
library(TwoSampleMR)
library(MRInstruments)
library(dplyr)


##
## List available GWASs in MR-Base catalogue

ao <- available_outcomes()


##
## Exposure info file

exp.info <- readxl::read_excel("data/exp_info.xlsx", sheet = "iv_exp")


##
## NMR metabolites info file

nmr.info <- read_excel("data/biomarker_annotations_updated.xlsx")


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
  mutate(analysis = "IVW (unadjusted)")


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
  mutate(analysis = "IVW (adjusted for BMI)") %>%
  rename(b = Estimate, se = "Std. Error", pval = "Pr(>|t|)") %>%
  filter(exposure == "Age at menarche")


##
## Save MR results

writexl::write_xlsx(list(mvmr_res = mvmr_res), path = "./output/mvmr_results.xlsx")

q("no")		
