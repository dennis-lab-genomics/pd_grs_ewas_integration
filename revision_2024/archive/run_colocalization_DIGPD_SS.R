# Created by Will Casazza
# Last updated Dec 29, 2024 by Samantha Schaffner

"
Run colocalization between GWAS loci and mQTL.
Usage:
  run_colocalization.R <mqtl> <gwas> <loci> <out> -n SAMPLES [options]
  run_colocalization.R (-h | --help)

Options:
  -h --help     Show this screen.
  -n SAMPLES    Number of Samples in mQTL study.
  --qc_fmt      Use QC/plink format GWAS summary statistics
  --plink_loci  Use Plink .clumped format for loci
  --smr_mqtl  Use smr query format for mQTL
" -> doc
library(coloc)
library(data.table)
library(docopt)
library(glue)
library(parallel)
#setDTthreads(8)
#arguments <- docopt(doc)
#setDTthreads(1)
 arguments <- list(
 n = "222",
 qc_fmt = TRUE,
 plink_loci = TRUE,
 #brain_mqtl = FALSE, #not used
 smr_mqtl = FALSE,
 mqtl = "~/PD_GRS/mQTL/cis_mQTL_out_DIGPD.txt",
 gwas = "~/SHARE_DECIPHER/nalls_PD.QC.gz",
 out = "cross_mqtl_cross_sumstats_digpd",
 loci = "~/SHARE_DECIPHER/nalls_PD.clumped"
 )
 #print(arguments)
 
# Grab CpG sites within 1MB of risk loci
probes_1mb <- function(row) {
  probe_pos[chr == paste0("chr", row$CHR) & s1 >= row$POS - 1e6 & s1 <= row$POS + 1e6]$geneid
}

localize_per_locus <- function(i,
                               mQTL_results = mQTL_results,
                               snp_pos = snp_pos,
                               has_methy = TRUE,
                               N1 = 222, #SS: what is N1? Should it be the same as "n" (218)?
                               sum_stats = sum_stats,
                               risk_scores = risk_scores # Must contain CHR,POS, SNP, and Locus Number column
) {
  all_probes <- unlist(probes_1mb(risk_scores[i, ]))
  probe_res <- mclapply(
    1:length(all_probes),
    function(j) {
      snp_set <- mQTL_results[all_probes[j], on = "gene"][!duplicated(SNP)]
      coords <- snp_pos[snp_set$SNP, on = "SNP"][!duplicated(SNP)]
      #tmp_stats <- na.omit(sum_stats[paste0(coords$CHR, ":", coords$POS), on = "SNP"])
      ## SS update 12/29/2024: format does not match "SNP" column, empty data table was returned here
      #sum_stats$SNP_pos <- paste0("chr",sum_stats$CHR, ":", sum_stats$POS)
      tmp_stats <- na.omit(sum_stats[coords$SNP, on = "SNP"])
      
      # match up data
      coords <- coords[match(tmp_stats$SNP, coords$SNP),] #SS 12/29/2024: changed "SNP" to "SNP_pos" here
      snp_set <- snp_set[match(coords$SNP, SNP)]
      if (has_methy) {
        mqtl_pval <- data.frame(
          beta = snp_set$beta,
          
          #varbeta = snp_set$beta / snp_set$`t-stat`, #SS 12/29/24: this should be squared to calculate standard error
          #see https://cran.r-project.org/web/packages/coloc/vignettes/a02_data.html
          #"varbeta is the standard error of beta, squared"
          
          varbeta = (snp_set$beta / snp_set$`t-stat`)^2,
          
          N = N1,
          
          #snp = paste0(coords$CHR, ":", coords$POS),
          snp=snp_set$SNP, #SS 12/29 updated
          
          #sdY = sd(methy[all_probes[j], -c(1), on = "cpg"]),
          sdY = sd(snp_set[all_probes[j], -c(1), on = "gene"]), #SS 12/29 updated to match df and column names
          
          type = "quant",
          stringsAsFactors = F
        )
      } else {
        mqtl_pval <- data.frame(
          beta = snp_set$beta,
          N = N1,
          #snp = paste0(coords$CHR, ":", coords$POS),
          snp = snp_set$SNP, #SS 12/29: changed to SNP rsid so format matches "gwas_pval"
          MAF = tmp_stats$MAF, #SS 12/29: changed from "freq" to "MAF" to match naming of tmp_stats column
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
        MAF = tmp_stats$MAF, #SS 12/29: changed from "freq" to "MAF" to match naming of tmp_stats column
        stringsAsFactors = F
      )
      mqtl_pval <- as.list(mqtl_pval[!is.na(gwas_pval$pvalues), ])
      mqtl_pval$N <- unique(mqtl_pval$N)[1]
      mqtl_pval$type <- unique(mqtl_pval$type)[1]
      gwas_pval <- as.list(gwas_pval[!is.na(gwas_pval$pvalues), ])
      gwas_pval$N <- unique(gwas_pval$N)[1]
      gwas_pval$type <- unique(gwas_pval$type)[1]
      if (length(gwas_pval$snp) != 0 | length(mqtl_pval$snp != 0)) {
        sink()
        res <- coloc.abf(mqtl_pval, gwas_pval)
        #minimum p value is: 0.02655
        #If this is what you expected, this is not a problem.
        #If this is not as small as you expected, please check you supplied var(beta) and not sd(beta) for the varbeta argument.
        
        sink()
        print("coloc completed")
        return(
          list(
            summary = cbind(t(data.frame(res$summary)), data.frame(probe = all_probes[j], locus = risk_scores$`Locus Number`[i], locus_snp = risk_scores$SNP[i])),
            result = cbind(res$results, data.frame(probe = all_probes[j], locus = risk_scores$`Locus Number`[i], locus_snp = risk_scores$SNP[i]))
          )
        )
      } else {
        print("returning blank data")
        return(list(summary = data.frame(), result = data.frame()))
      }
    },
    mc.cores = 16 #changed from 16 for testing
  )
  return(
    list(
      summary = rbindlist(mclapply(probe_res, function(res) res$summary, mc.cores = 4)),
      result = rbindlist(mclapply(probe_res, function(res) res$result, mc.cores = 4))
    )
  )
}


sum_stats <- fread(arguments$gwas)
mQTL_results <- fread(arguments$mqtl)
loci <- fread(arguments$loci)

if (arguments$qc_fmt) {
  if ("BP" %in% colnames(sum_stats)) {
    sum_stats[, `:=`(SNP = paste0("chr", CHR, ":", BP), freq = MAF)]
  } else {
    sum_stats[, `:=`(SNP = paste0("chr", CHR, ":", POS), freq = MAF)]
  }
}

if (arguments$plink_loci) {
  loci[, `:=`(POS = BP, `Locus Number` = .I)]
}

if (arguments$smr_mqtl) {
  snp_pos <- unique(mQTL_results[, .(
    SNP = SNP,
    CHR = paste0("chr", Chr),
    POS = BP
  )])


  mQTL_results <- mQTL_results[, .(
    SNP = SNP,
    gene = Probe,
    beta = b,
    `p-value` = sapply(p, max, .Machine$double.xmin)
  )]
  print("SNP pos table")
  print(head(snp_pos))
  print("mQTL table")
  print(head(mQTL_results))
} else {
  snp_pos <- fread("~/PD_GRS/mQTL/snp_pos_DIGPD.txt", key = "SNP")
}

setkey(mQTL_results, "gene")
probe_pos <- fread("~/PD_GRS/mQTL/probe_pos_DIGPD.txt")
possible_probes <- mQTL_results[, .(gene, minP = min(`p-value`)), by = "gene"][minP < 5e-8]$gene
probe_pos <- probe_pos[possible_probes, on = "geneid"]

#runtime approx 1.5 hours
system.time(
  result <- mclapply(
    1:nrow(loci), 
    function(i) {
      res <- tryCatch(localize_per_locus(
        i,
        mQTL_results = mQTL_results[possible_probes, on = "gene"],
        snp_pos = snp_pos,
        has_methy = FALSE,
        risk_scores = loci,
        sum_stats = sum_stats,
        N1 = as.numeric(arguments$n)
      ),
      error = function(e) e
      )
      if (inherits(res, "error")) {
        return(list(summary = data.table(), result = data.table()))
      } else {
        return(res)
      }
    },
    mc.cores = 4
  )
)

#from coloc.abf documentation:
#summary is a vector giving the number of SNPs analysed, and the posterior probabilities of:
#H0 (no causal variant), 
#H1 (causal variant for trait 1 only), 
#H2 (causal variant for trait 2 only), 
#H3 (two distinct causal variants) 
#and H4 (one common causal variant)
fwrite(
  rbindlist(mclapply(result, function(res) res$summary, mc.cores = 4)),
  glue("{arguments$out}_digpd_pd_snp_colocalization_ph4.txt.gz"),
  sep = "\t", row.names = F, quote = F
)

#from coloc.abf documentation:
#results is an annotated version of the input data containing 
#log Approximate Bayes Factors and intermediate calculations, and the 
#posterior probability SNP.PP.H4 of the SNP being causal for the shared signal if H4 is true. 
#This is only relevant if the posterior support for H4 in summary is convincing.
fwrite(
  rbindlist(mclapply(result, function(res) res$result, mc.cores = 4)),
  glue("{arguments$out}_digpd_pd_snp_colocalization_per_snp.txt.gz"),
  sep = "\t", row.names = F, quote = F
)
