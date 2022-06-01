############################################################################
#                                Set - UP                                  #
############################################################################


##
## Set working directory 

setwd(paste0(Sys.getenv('SCRATCH'), "repro_health/"))


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

update.packages(ask = "FALSE")


##
## Install packages

#install.packages(c("data.table", "purrr", "devtools", "dplyr", "ggplot2"))
#devtools::install_github("MRCIEU/TwoSampleMR")
#devtools::install_github("MRCIEU/MRInstruments")
#devtools::install_github("NightingaleHealth/ggforestplot")


##
## Load libraries

library(data.table)
library(purrr)
library(TwoSampleMR)
library(MRInstruments)
library(dplyr)
library(ggplot2)
library(ggforestplot)


##
## List available GWASs in MR-Base catalogue

ao <- available_outcomes()


##
## Exposure info file

exp.info <- readxl::read_excel("data/exp_info.xlsx", sheet="iv_exp") %>%
  mutate_at(vars(starts_with("N")), list(~as.character(.)))


############################################################################
#                                IV selection                              #
############################################################################

##
## Format data from MR-Base

# MR-Base IDs
mrbase.ids <- exp.info %>%
  filter(source == "mr-base") %>%
  pull(id) %>% 
  as.character

# MR-Base SNP-exposure summary data
mrbase.expdat <- extract_instruments(outcomes = mrbase.ids)


##
## Format data from original paper

# Number of children ever born (NEB)

neb.file <- filter(exp.info, abbrev == "NEB") %>%
  pull(file) %>%
  as.character

neb.expdat <- list.files(path = "./data/NEB/", pattern = neb.file, full.names = T) %>%
  read.table(., header = T, sep = "\t") %>%
  filter(Strata == "NEB-Sex Combined" | Strata == "NEB-Female") %>%
  mutate(Pos = gsub(",", "", Pos)) %>%
  mutate(Pos = as.numeric(Pos)) %>%
  format_data(., 
              type="exposure",
              snp_col = "SNP",
              beta_col = "Female.NEB.Beta",
              se_col = "Female.NEB.SE",
              eaf_col = "Inc_Freq",
              effect_allele_col = "Effect_Allele",
              other_allele_col = "Other_Allele", 
              pval_col = "Female.NEB.Pvalue",
              chr_col = "Chr",
              pos_col = "Pos"
  ) %>%
  mutate(., exposure = "NEB", units.exposure = "n", id.exposure = "neb-st7", chr.exposure = as.character(chr.exposure), samplesize.exposure = 478624) # N females in discovery analyses

# 32 NEB hits = 28 sex-combined + 4 female-specific  
# SNP-NEB effect estimates from female metanalysis

# Age at menarche - 2017 - 389 SNPs - CAREFUL AS ALLELES FOR DELETIONS/INSERTIONS ARE CODED AS d/i

aam.file <- filter(exp.info, abbrev == "AAM" & ukb_included == "yes") %>%
  pull(file) %>%
  as.character

aam.id <- filter(exp.info, abbrev == "AAM" & ukb_included == "yes") %>%
  pull(id) %>%
  as.character

aam.expdat <- list.files(path = "./data/AAM/", pattern = aam.file, full.names = T) %>%
  readxl::read_excel(., sheet="S2 - Genome-wide Signals")  %>%
  format_data(., 
              type = "exposure",
              snp_col = "SNP",
              beta_col = "BETA (y/allele)",
              se_col = "SE",
              eaf_col = "EAF",
              effect_allele_col = "Effect allele",
              other_allele_col = "Other allele", 
              pval_col = "P",
              chr_col = "Chr.",
              pos_col = "Position (b37)"
  ) %>%
  mutate(., exposure = "Age at menarche", id.exposure = aam.id, chr.exposure = as.character(chr.exposure), samplesize.exposure = 329345)

# Age at natural menopause - 2021 - 290 SNPs (I used effect estimates from discovery+replication except for rs4800141 and rs6490269 as effect estimates were only available from discovery sample)

anm.file <- filter(exp.info, abbrev == "ANM" & ukb_included == "yes") %>%
  pull(file) %>%
  as.character

anm.id <- filter(exp.info, abbrev == "ANM" & ukb_included == "yes") %>%
  pull(id) %>%
  as.character

anm.expdat <- list.files(path = "./data/ANM/", pattern = anm.file, full.names = T) %>%
  read.table(., header = T, sep = "\t")  %>%
  format_data(., 
              type = "exposure",
              snp_col = "SNP",
              beta_col = "Effect",
              se_col = "SE",
              eaf_col = "EAF",
              effect_allele_col = "Effect.Allele",
              other_allele_col = "Other.Allele", 
              pval_col = "P.value",
              chr_col = "Chr",
              pos_col = "Pos..b37."
  ) %>%
  mutate(., exposure = "Age at natural menopause", id.exposure = anm.id, chr.exposure = as.character(chr.exposure), samplesize.exposure = 496151) 
# P.S.: As per ST2, all effects are on age at natural menopause in years


############################################################################
#                      Combine all SNP-exposure data                       #
############################################################################

##
## Combine all SNP-exposure data

expdat <- bind_rows(mrbase.expdat, neb.expdat, aam.expdat, anm.expdat) %>%
  merge(., exp.info, by.x = "id.exposure", by.y = "id") 

write.table(expdat, file = "data/expdat.txt")			


##
## Save list of unique rsIDs across exposures

write.table((unique(expdat$SNP)), file = "data/rsids.txt", quote = F, row.names = F, col.names = F)


##
## SNP-exposure data

expdat2 <- expdat %>%
  rename_at(vars(ends_with(".exposure")), ~sub("[.]exposure", "", .)) %>%														 
  select(id, trait, abbrev, SNP, chr, pos, effect_allele, other_allele, eaf, beta, se, pval, units, samplesize) %>%
  mutate(., f = ((beta/se)^2), r2 = f/(samplesize-2+f))


##
## Summary information on each exposure

iv_strength <- expdat2 %>%
  group_by(id) %>%
  summarise(.,
            id       = unique(id),
            Trait    = unique(trait),
            Nsnps    = n(), 
            total_R2 = round(sum(r2), 3), 
            mean_F   = round(mean(f), 0),
            N        = unique(samplesize)
  ) 


##
## Study information on each exposure																	

iv_source <- expdat %>%
  group_by(id.exposure) %>%
  count %>%
  rename(nSNPs = n) %>%
  merge(exp.info, ., by.x = "id", by.y = "id.exposure") %>%
  select(trait, abbrev, pmid, doi, unit, source, id, nSNPs, N, N_cases, N_ctrls, N_ukbb, N_ukbb_cases, N_ukbb_ctrls) %>%
  mutate_at(., vars(starts_with("N_")),~as.numeric(.)) 


##
## Save metadata

writexl::write_xlsx(list(
  snpexp = expdat2, 
  ivstrength = iv_strength, 
  ivsource = iv_source
), 
path = "./output/SupplTables.xlsx")



##
## Check number of SNPs per exposure

print(table(expdat$exposure))

q("no")
