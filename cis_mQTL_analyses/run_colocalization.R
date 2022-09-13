'
Run colocalization between GWAS loci and mQTL.
Usage:
  run_colocalization.R <mqtl> <gwas> <loci> <out> -n SAMPLES [options]
  run_colocalization.R (-h | --help)

Options:
  -h --help     Show this screen.
  -n SAMPLES    Number of Samples in mQTL study.
  --qc_fmt      Use QC/plink format GWAS summary statistics
  --plink_loci  Use Plink .clumped format for loci
  --brain_mqtl  Use format for ROSMAP brain mQTL
' -> doc
library(coloc)
library(data.table)
library(docopt)
library(glue)
library(parallel)
arguments <- docopt(doc)
print(arguments)
# Grab CpG sites within 1MB of risk loci
probes_1mb <- function(row) {
  probe_pos[chr == paste0("chr", row$CHR) & s1 >= row$POS - 1e6 & s1 <= row$POS + 1e6]$geneid
}

localize_per_locus <- function(
  i,
  mQTL_results = mQTL_results,
  snp_pos = snp_pos,
  has_methy = TRUE,
  N1 = 245,
  sum_stats=sum_stats,
  risk_scores=risk_scores # Must contain CHR,POS, SNP, and Locus Number column
){
  all_probes <- unlist(probes_1mb(risk_scores[i, ]))
  probe_res <- mclapply(
    1:length(all_probes),
    function(j) {
      snp_set <- mQTL_results[all_probes[j], on = "gene"][!duplicated(SNP)]
      coords <- snp_pos[snp_set$SNP, on = "SNP"][!duplicated(SNP)]
      tmp_stats <- sum_stats[paste0(coords$CHR, ":", coords$POS), on = "SNP"]
      if (has_methy) {
        mqtl_pval <- data.frame(
          beta = snp_set$beta,
          varbeta = snp_set$beta / snp_set$`t-stat`,
          N = N1,
          snp = paste0(coords$CHR, ":", coords$POS),
          sdY = sd(methy[all_probes[j], -c(1), on = "cpg"]),
          type = "quant",
          stringsAsFactors = F
        )
      } else {
        mqtl_pval <- data.frame(
          beta = snp_set$beta,
          N = N1,
          snp = paste0(coords$CHR, ":", coords$POS),
          MAF = tmp_stats$freq,
          pvalues = snp_set$`p-value`,
          type = "quant",
          stringsAsFactors = F
        )
      }
      gwas_pval <- data.frame(
        pvalues = tmp_stats$p,
        N = tmp_stats$N_cases + tmp_stats$N_controls,
        s = tmp_stats$N_cases / (tmp_stats$N_cases + tmp_stats$N_controls),
        snp = tmp_stats$SNP,
        type = "cc",
        MAF = tmp_stats$freq,
        stringsAsFactors = F
      )
      mqtl_pval <- as.list(mqtl_pval[!is.na(gwas_pval$pvalues), ])
      mqtl_pval$N <- unique(mqtl_pval$N)[1]
      mqtl_pval$type <- unique(mqtl_pval$type)[1]
      gwas_pval <- as.list(gwas_pval[!is.na(gwas_pval$pvalues), ])
      gwas_pval$N <- unique(gwas_pval$N)[1]
      gwas_pval$type <- unique(gwas_pval$type)[1]
      print(gwas_pval)
      print(mqtl_pval)
      if (length(gwas_pval$snp) != 0 | length(mqtl_pval$snp != 0)) {
        res <- coloc.abf(mqtl_pval, gwas_pval)
        return(
          list(
            summary=cbind(t(data.frame(res$summary)), data.frame(probe = all_probes[j], locus = risk_scores$`Locus Number`[i], locus_snp = risk_scores$SNP[i])),
            result = cbind(res$results, data.frame(probe = all_probes[j], locus = risk_scores$`Locus Number`[i], locus_snp = risk_scores$SNP[i]))
          )
        )
      } else {
        return(data.frame())
      }
    },
    mc.cores = 4
  )
  return(rbindlist(probe_res))
}


sum_stats <- fread(arguments$gwas)
mQTL_results <- fread(arguments$mqtl, key="gene")
loci <- fread(arguments$loci)

if(arguments$qc_fmt)
  sum_stats[,`:=`(SNP=paste0("chr",CHR,":",BP),freq=MAF)]

if(arguments$plink_loci)
  loci[,`:=`(POS=BP, `Locus Number`=.I)]

if(arguments$brain_mqtl){
  mQTLs_results <- mQTLs_results[, .(
    SNP = SNPid,
    gene = featureName,
    beta = SpearmanRho,
    `p-value` = pValue
  )]
  setkey(mQTLs_brain, "gene")
  snp_pos <- fread("pos_brain.txt.gz")

}else{
  snp_pos <- fread("terre_data/snp_pos.txt", key = "SNP")
}
probe_pos <- fread("terre_data/probe_pos.txt")

system.time(
  result <- rbindlist(
    mclapply(
      1:nrow(loci),
      function(i) localize_per_locus(
        i,
        mQTL_results = mQTL_results,
        snp_pos = snp_pos,
        has_methy = FALSE,
        N1 = arguments$SAMPLES
      ),
      mc.cores = 4
    )
  )
)
fwrite(result$summary, glue("{out}_pd_snp_colocalization_ph4.txt.gz"), sep = "\t", row.names = F, quote = F)
fwrite(result$result, glue("{out}_pd_snp_colocalization_per_snp.txt.gz"), sep = "\t", row.names = F, quote = F)
