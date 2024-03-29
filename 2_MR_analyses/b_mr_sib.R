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

library(dplyr)
library(data.table)
library(purrr)
library(TwoSampleMR)
#library(MRInstruments)
#library(xlsx)
library(readxl)
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

sib.info <- read_excel("data/sibship/SibGWAS_traits.xlsx", sheet = 1)
 

############################################################################
#                       Prepare and harmonise data                         #
############################################################################

##
## Exposures IDs

expid <- c(AAM = "ST3_41588_2017_BFng3841_MOESM301_ESM", NEB = "neb-st7", ANM = "ST2_41586_2021_3779_MOESM3_ESM")


##
## Selected SNPs 

ivdat <- read.table("data/expdat.txt") %>% 
    			mutate_at(vars(id.exposure, SNP), ~as.character(.)) %>%
    			filter(id.exposure %in% expid)
			
snps <-	map(expid, ~filter(ivdat, id.exposure %in% .)) %>%
		      map(., ~pull(., SNP))
			

##
## Create SNP link file (ie chr:position linked to rsID)

table(is.na(ivdat$chr.exposure)) # no missing

table(is.na(ivdat$pos.exposure)) # no missing

snp_link <- mutate(ivdat, snpid = paste0("chr", chr.exposure, ":", pos.exposure, ":SNP")) %>%
        			select(SNP, snpid) %>%
        			distinct(SNP, .keep_all = T)
			
			
##
## Function to read raw data

read_rawdat <- function(x) {
	
		f <- filter(sib.info, abbrev == x) %>% pull(file)
			
		df <- paste0("data/sibship/", f, ".summary.gz") %>%
				fread(.) %>%
				filter(SNP %in% snp_link$snpid) %>%
				mutate(file = f) %>%
				merge(snp_link, ., by.x = "snpid", by.y = "SNP") %>%
				merge(., sib.info, by = "file")
}
#SNP = Marker name in CHR:BP format
#CHR = Chromosome
#BP = Base pair position
#A1 = Effect allele
#A2 = Non-effect allele
#WS = within-sibship estimates
#Pop = population estimates
#Eff_N = effective N based on GWAS standard errors


##
## Exposure data

expdat <- map(names(snps), ~read_rawdat(.)) %>%
			      map2(., snps, ~filter(.x, SNP %in% .y))


##
## Outcome data

out <- filter(sib.info, ! abbrev %in% names(expid)) %>% pull(abbrev)

outdat <- map(out, ~read_rawdat(.)) 


##
## Function to format raw data 

tide_outdat <- function(x, y, z) {

		df <- format_data(
					  x,
					  type = y,
					  phenotype_col = "trait",
					  snp_col = "SNP",
					  effect_allele_col = "A1",
					  other_allele_col = "A2",
					  units_col = "unit",
					  beta_col = paste0("BETA_", z),
					  se_col = paste0("SE_BETA_", z),
					  pval_col = paste0("P_BETA_", z),
					  samplesize_col = paste0("Eff_N_", z)
					) 
}

##
## Clean and harmonise within-sib data

sib_expdat <- map(expdat, ~tide_outdat(x = ., y = "exposure", z = "WS")) %>% bind_rows
sib_outdat <- map(outdat, ~tide_outdat(x = ., y = "outcome", z = "WS")) %>% bind_rows
sib_dat <- harmonise_data(sib_expdat, sib_outdat) %>% mutate(subsamp = "Within-siblings")


##
## Clean and harmonise unrelateds data

pop_expdat <- map(expdat, ~tide_outdat(x = ., y = "exposure", z = "Pop")) %>% bind_rows
pop_outdat <- map(outdat, ~tide_outdat(x = ., y = "outcome", z = "Pop")) %>% bind_rows
pop_dat <- harmonise_data(pop_expdat, pop_outdat) %>% mutate(subsamp = "Unrelateds")


############################################################################
#                                Run MR                                    #
############################################################################

# within-sib MR results
sib_mr <- mr(sib_dat, method_list=c("mr_ivw")) %>% mutate(subsamp = "Within-siblings")

# Unrelateds MR results
pop_mr <- mr(pop_dat, method_list=c("mr_ivw")) %>% mutate(subsamp = "Unrelateds")

# Combine MR results
mr <- bind_rows(sib_mr, pop_mr)


############################################################################
#                              Plots                                       #
############################################################################

##
## Format MR results

mr_res <- mutate(mr, 
					exposure = factor(exposure, levels=c("Age at menarche", "Parity", "Age at menopause")),
					outcome = factor(outcome, levels=c("Glycated haemoglobin", "C-reactive protein", "HDL-cholesterol", "LDL-cholesterol", "Triglycerides"))
				)


##
## Forest plot - loci vs CVD risk factors

p.sib <- ggforestplot::forestplot(
                                  df = mr_res,
                                  name = outcome,
                                  estimate = b,
                                  se = se,
                                  pvalue = pval,
                                  psignif = 0.05,
                                  xlab = "SD units",
                                  #title = "",
                                  colour = subsamp,
                                  shape = subsamp,
                                  logodds = F	
                                ) +							
								theme(legend.title = element_blank(),
									  axis.title.x = element_text(size = 14),
									  axis.text.x = element_text(size = 14),
									  axis.text.y = element_text(size = 14),
									  legend.text = element_text(size = 14),
									  strip.text.x = element_text(size = 16, hjust = 0.5)							
								) +
								scale_color_manual(values=c("darkorange1", "deepskyblue2")) +
								ggforce::facet_row(
									facets = ~exposure,
									scales = "free_x",
									space = "fixed"
								  )

										
##
## Save plot
	
ggsave(paste0("results/plots/forest_repro_vs_riskfactors_SIBSHIP.png"), plot = p.sib, width = 30,  height = 20, units = "cm")

q("no")



		
 
