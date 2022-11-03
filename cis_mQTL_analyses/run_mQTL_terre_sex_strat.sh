#!/bin/bash


for SEX in  "female" "all"; do
    /usr/local/R-3.6.2/bin/Rscript run_matrixEQTL.R --data_dir="$PWD/" \
      --SNP_fname="${SEX}_geno_2022.txt" \
      --cov_file="${SEX}_meta_2022_prs.txt" \
      --cis_outfile="${SEX}_2022_prs_mQTL.txt" \
      --snp_pos="terre_data/snp_pos.txt" \
      --methylation_fname="${SEX}_methy_2022.txt" \
      --probe_pos="terre_data/probe_pos.txt"

    gzip "${SEX}_2022_mQTL.txt"
done


