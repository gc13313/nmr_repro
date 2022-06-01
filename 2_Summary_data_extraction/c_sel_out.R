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

#install.packages(c("data.table", "purrr", "dplyr"))


##
## Load libraries

library(data.table)
library(purrr)
library(dplyr)


############################################################################
#               Extract SNP-metabolites summary data                       #
############################################################################

##
## Import and format data extract using bash script ("run_mr_setup.sh" - lines 65-154)

import_dat <- function(type) {

		fpath <- list.files(path = paste0(Sys.getenv('HOME'), "/tmp_nmr_ukb"), pattern = paste0("^", type), full.names = T)

		fname <- list.files(path = paste0(Sys.getenv('HOME'), "/tmp_nmr_ukb"), pattern = paste0("^", type), full.names = F)

		df <- map(fpath, ~fread(.)) %>%
				map2(., fname, ~mutate(.x, id = .y)) %>%
				map(., ~mutate_all(., ~as.character(.))) %>%
				bind_rows

}

# Extract NMR data (females)
nmr_dat_female <- import_dat(type = "nmr_dat_female")

# Save extracted NMR data
write.table(nmr_dat_female, "data/nmr/nmr_female_sumdat_ivs.txt")

q("no")
