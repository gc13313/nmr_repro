# nmr_repro
A project titled " From menarche to menopause: the impact of reproductive factors on the metabolic profile of over 65,000 women". This repository includes the analysis plan and all code related to the manuscript (currently on Medrxiv - doi: https://doi.org/10.1101/2022.04.17.22273947).

This repository contains:

**1) Directory: 1_MV_analyses**

Scripts:

•	a_gen_rel_var.do: Generate relevant variables (i.e each reproductive traits, confounders)

•	b_gen_statins.do: Generate statin variables

•	c_gen_statins_ds.do: Generate statin dataset

•	d_merge_outcome_ds.do: Merge biomarker outcomes dataset with exposure dataset

•	e_mvr_allexposures_allmods.do: Run all MVR models

•	f_pnc_parity_allmods.do: Run PNC models

•	g_output_plots.do: Merge MVR and PNC results

•	h_output_sumtable.do: Read in MR results and MVR/PNC results and output to supplementary excel table

**2) Directory: 2_Summary_data_extraction**

Scripts:

•	a_sel_iv.R: Selects SNPs as instruments for reproductive traits and extracts corresponding exposure GWAS summary data

•	b_extr_outdat.sh: Extracts outcome GWAS summary data for selected SNPs

•	c_sel_out.R: Formats outcome GWAS summary data for selected SNPs

**3)	Directory: 3_MR_analyses**

Scripts:

•	a_mr.R: Performs MR analyses (IVW, weighted median and MR-Egger)

•	b_mr_sib.R: Performs MR using within-siblings summary GWAS data

•	c_mvmr.R: Performs multivariable MR of age at menarche on metabolic traits accounting for BMI

•	d_mr_negcontrl.R: Performs MR analyses using negative outcome control 

