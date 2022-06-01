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
## Biomarker info file (MR-Base IDs for biomarker panel in 500K UKB participants - males+females)

out.info <- readxl::read_excel("data/ukb_biomarker_panel.xlsx")


##
## Exposure info file

exp.info <- read_excel("data/exp_info.xlsx", sheet = "iv_exp")


############################################################################
#                  NEGATIVE CONTROL MR - PARITY X SKIN COLOUR              #
############################################################################

##
## Exposure data (AAM, NEB, ANM)

expdat <- read.table("data/expdat.txt") %>% 
  mutate_at(., vars(SNP, exposure), ~as.character(.))


##
## Outcome data

outname <- as.list(c("skin_colour", "skin_tanning")) 

nc_outdat <- outname %>%
  map(., ~list.files(path = paste0(Sys.getenv('RDSF_IEU2'), "/p6/094/working/results/GWAS/UKB/mothers"), pattern = ., full.names = T)) %>%
  map(., ~fread(.)) %>%
  map(., ~format_data(.,
                      type = "outcome",
                      snps = expdat$SNP,
                      snp_col = "snp",
                      beta_col = "beta",
                      se_col = "se",
                      eaf_col = "eaf",
                      effect_allele_col = "effect_allele",
                      other_allele_col = "other_allele",
                      pval_col = "pval"
  )) %>%
  map2(., outname, ~mutate(.x, outcome = .y)) %>%
  bind_rows


##
## Harmonise exposure and outcome data

nc_dat <- harmonise_data(expdat, nc_outdat)


##
## Run MR

nc_res <- mr(nc_dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median")) 

writexl::write_xlsx(list(nc_res = nc_res), path = "./output/negcontrol_results.xlsx")

