#!/bin/bash
run_pc () {
  /usr/local/R-3.6.2/bin/Rscript run_matrixEQTL.R --data_dir="terre_data/" \
    --SNP_fname="all_imputed_matrixeQTL_chr21.txt"\
    --cov_file="covariates_${1}_methy_PC.txt"\
    --cis_outfile="cis_all_impute_mQTL_results_${1}_methy_PC_chr21.txt"\
    --snp_pos="snp_pos_chr21.txt"\
    --methylation_fname="methylation_combat_chr21.txt"\
    --probe_pos="probe_pos_chr21.txt"
}

export -f run_pc
~/parallel-20201022/src/parallel run_pc {} ::: {0..20}

/usr/local/R-3.6.2/bin/Rscript run_matrixEQTL.R --data_dir="terre_data/" \
  --SNP_fname="all_imputed_matrixeQTL_chr21.txt"\
  --cov_file="covariates_CTP.txt"\
  --cis_outfile="cis_all_impute_mQTL_results_CTP_chr21.txt"\
  --snp_pos="snp_pos_chr21.txt"\
  --methylation_fname="methylation_combat_chr21.txt"\
  --probe_pos="probe_pos_chr21.txt"