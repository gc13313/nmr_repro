#!/bin/bash

#SBATCH --partition=mrcieu 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

##
## Record some potentially useful details about the job:
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo "Extract metabolites GWAS summary data ..."


##
## Scripts directory
scripts_dir="${SCRATCH}repro_health/scripts/"
cd $scripts_dir


############################################################################
#                   Extract SNP-outcome estimates                          #
############################################################################

##
## Set a temporary working directory
export WORK_DIR="${HOME}/tmp_nmr_ukb"
mkdir -p "${HOME}/tmp_nmr_ukb"
cd $WORK_DIR


##
## List of rsids
scp mb16066@bc4login3:${SCRATCH}repro_health/data/rsids.txt .
			
##
## Function to extract SNP-outcome data for selected IVs

extrdat() {

	# Path to input files
	dir_file=${RDSF_IEU2}/p6/119/working/data/$1
	echo ${dir_file}
	
	# List of files
	files=$(ssh mb16066@bc4login3 ls ${dir_file}/*imputed.txt.gz)
	echo ${files}
		
		for fpath in ${files[@]}
		do
		
				# Strip directory from filename
				f=`basename "${fpath}"`
				echo ${f}
				
				# Copy the input file to WD 
				scp mb16066@bc4login3:${fpath} .
								
				# Create name for output file
				outfile=$1\_${f}_TMPDAT.txt
				echo ${outfile}			
				
				# Copy header to output file
				zcat ${f} | awk NR==1 > ${outfile} 
				
				# Extract data for selected IVs and copy to output file
				zcat ${f} | grep -w -F -f rsids.txt >> ${outfile}
				# -w tells grep to match whole words only
				# -F search for fixed strings (plain text) rather than regular expressions
				# -f "data/rsids.txt" read search patterns from the file
				
				# Remove copy of input file
				rm ${f}
		done
}


# NMR metabolites 
extrdat nmr_dat_female
