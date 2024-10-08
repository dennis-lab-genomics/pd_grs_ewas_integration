#!/bin/bash
/usr/local/R-3.6.2/bin/Rscript run_matrixEQTL.R --data_dir="digpd_data/" \
	--SNP_fname="all_imputed_matrixeQTL_no_mut.txt" \
	--cov_file="covariates_9_methy_PC_no_mut.txt" \
	--cis_outfile="cis_all_impute_mQTL_results_9_methy_PC_no_mut.txt" \
	--snp_pos="snp_pos.txt" \
	--methylation_fname="methylation_combat_no_mut.txt" \
	--probe_pos="probe_pos.txt"
