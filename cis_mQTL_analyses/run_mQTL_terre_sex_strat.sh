#!/bin/bash
/usr/local/R-3.6.2/bin/Rscript run_matrixEQTL.R --data_dir="terre_data/" \
	--SNP_fname="female_all_imputed_matrixeQTL.txt" \
	--cov_file="female_covariates_9_methy_PC.txt" \
	--cis_outfile="female_cis_all_impute_mQTL_results_9_methy_PC.txt" \
	--snp_pos="snp_pos.txt" \
	--methylation_fname="female_methylation_combat.txt" \
	--probe_pos="probe_pos.txt"

gzip terre_data/female_cis_all_impute_mQTL_results_9_methy_PC.txt

/usr/local/R-3.6.2/bin/Rscript run_matrixEQTL.R --data_dir="terre_data/" \
	--SNP_fname="male_all_imputed_matrixeQTL.txt" \
	--cov_file="male_covariates_9_methy_PC.txt" \
	--cis_outfile="male_cis_all_impute_mQTL_results_9_methy_PC.txt" \
	--snp_pos="snp_pos.txt" \
	--methylation_fname="male_methylation_combat.txt" \
	--probe_pos="probe_pos.txt"

gzip terre_data/male_cis_all_impute_mQTL_results_9_methy_PC.txt
