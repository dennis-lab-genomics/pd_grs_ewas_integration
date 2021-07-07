#!/bin/bash
/usr/local/R-3.6.2/bin/Rscript run_matrixEQTL.R --data_dir="digpd_data/" \
  --SNP_fname="digpd_data/all_imputed_matrixeQTL_CTP.txt"\
  --cov_file="covariates_CTP.txt"\
  --cis_outfile="cis_all_impute_mQTL_results_CTP.txt"\
  --snp_pos="snp_pos_CTP.txt"\
  --methylation_fname="methylation_combat_CTP.txt"\
  --probe_pos="probe_pos_CTP.txt"

