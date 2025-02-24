#!/bin/bash

#Updated by Samantha Schaffner July 26, 2024
#Changed R to version 4.1.0 (from 3.6.2)
#Added newly generated DNAm and covariate data from July 26, 2024

/usr/local/R-4.1.0/bin/Rscript /home1/NEURO/schaffner/PD_GRS/mQTL/run_matrixEQTL_SS.R --data_dir="/home1/NEURO/schaffner/PD_GRS/mQTL/" \
	--SNP_fname="all_imputed_matrixeQTL.txt" \
	--cov_file="covariates_CTP_PD.txt" \
	--cis_outfile="cis_all_impute_mQTL_results.txt" \
	--snp_pos="snp_pos.txt" \
	--methylation_fname="methylation.txt" \
	--probe_pos="probe_pos.txt"

pigz /home1/NEURO/schaffner/PD_GRS/mQTL/cis_all_impute_mQTL_results.txt
