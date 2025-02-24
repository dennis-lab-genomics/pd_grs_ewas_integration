#!/bin/bash

#Created by Samantha Schaffner Jan 20, 2025
#this takes the list of all QC'ed SNPs and the list of GRS SNPs as input, and outputs SNPs in LD (r2 > 0.2) with the GRS

#first install plink2
#downloaded alpha 6.7 (20Jan2025) from https://www.cog-genomics.org/plink/2.0/ and uploaded to home directory
#(instead of wget, as cannot connect to web and INSERM server at once)

chmod +x /home1/NEURO/schaffner/plink2

#unphased: does not consider haplotype phase
#standard for GWAS and population genetics studies
#unsure which is reference allele; ref-unknown option treats first as reference

~/plink2 --gen /home1/NEURO/schaffner/SHARE_DECIPHER/pd_grs_ewas_integration/TERRE_QC/raw_data.imputed.r2_30.maf_mismatch.gen ref-unknown --sample /home1/NEURO/schaffner/SHARE_DECIPHER/pd_grs_ewas_integration/TERRE_QC/raw_data.imputed.r2_30.maf_mismatch.sample --r2-unphased --ld-snp-list /home1/NEURO/schaffner/PD_GRS/mQTL/prs_snps.txt --ld-window-r2 0.2 --out /home1/NEURO/schaffner/PD_GRS/mQTL/grs_snps_ld
