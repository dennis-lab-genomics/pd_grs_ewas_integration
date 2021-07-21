library(mice)
library(data.table)
library(tidyverse)
library(limma)

load("/home1/NEURO/SHARE_DECIPHER/processed_DNAm_data/TERRE/TERRE_betas_combat.RData")# betas_sub

terre_prs <- fread("terre_prs.sscore")
terre_pcs <- fread("~/genotype_qc/TERRE_QC/raw_data.geno.maf.mind.sex_check.het_filter.ibd_filter.eigenvec")
terre_metadata <- read.csv("/home1/NEURO/SHARE_DECIPHER/terre_meta_master.csv") %>%
  mutate(
    FID = gsub("(PAE_[0-9]*_[1-9]*)_.*","\\1",FID),
    IID = gsub("(PAE_[0-9]*_[1-9]*)_.*","\\1",IID)
  )
all_data <- terre_prs %>% left_join(terre_metadata, by=c("IID")) %>% left_join(terre_pcs,by=c("IID"= "V1"))
reg_all_data <- all_data[na.omit(match(colnames(betas_sub),all_data$patient)),]
reg_all_data$sex <- ifelse(reg_all_data$men,"Male","Female")
reg_all_data$PD <- ifelse(reg_all_data$PD,"CASE","CONTROL")
betas_sub <- betas_sub[,colnames(betas_sub) %in% reg_all_data$patient]


pest_missing <- read.csv("/home1/NEURO/SHARE_DECIPHER/TERRE_pesticides/pesticides.csv")
pest_imputed <- read.csv("/home1/NEURO/SHARE_DECIPHER/TERRE_pesticides/pesticides_imputed.csv")


pre_mids <-rbind(cbind(data.frame(X_imputation_=0),pest_missing)[,colnames(pest_imputed)],pest_imputed)
pre_mids <- pre_mids[pre_mids$num %in% colnames(betas_sub),]
pest_mids <- as.mids(pre_mids, .imp="X_imputation_",.id="num")
terre_pest_meta <- reg_all_data[match(rownames(pest_mids$data),reg_all_data$patient),]
betas_pest <- betas_sub[,rownames(pest_mids$data)]
exposures <- colnames(pest_mids$data)[colMeans(pest_mids$data,na.rm = T) > 0.025]

for(exposure in exposures){
  design_call <- substitute(
    with(pest_mids,
    model.matrix.lm(
      ~ E*SCORE1_AVG+V3+V4+V5+age+men,
      data = terre_pest_meta,
      na.action='na.pass')
    ),
    list(E=as.name(exposure)))
  design <- eval(design_call,list(E=as.name(exposure)))
  design <- design$analyses
  res <- parallel::mclapply(
    design,
    function(x){
      na_row <- is.na(rowSums(x))
      eBayes(lmFit(betas_pest[,!na_row],x[!na_row,]))
    },
    mc.cores = 10)
  coefs <- lapply(res,coef)
  vars <- lapply(res,function(x)(x$stdev.unscaled * sqrt(x$s2.post))^2)
  coefs_per_cpg <- lapply(
    1:nrow(betas_pest),
    function(i){
      unlist(lapply(coefs,function(mat) mat[i,9]))
    }
  )
  vars_per_cpg <- lapply(
    1:nrow(betas_pest),
    function(i){
      unlist(lapply(vars,function(mat) mat[i,9]))
    }
  )
  pooled_res <- parallel::mclapply(
    1:nrow(betas_pest),
    function(i){
      pooled <- pool.scalar(coefs_per_cpg[[i]],vars_per_cpg[[i]],ncol(betas_pest),k =9)
      out <- data.frame(pooled[-c(2,3)])
      out$P <- 2*pt(abs(out$t),out$df,lower.tail=F)
      return(out)
    },
    mc.cores=20
  )
  names(pooled_res) <- rownames(betas_pest)
  sum_stats <- rbindlist(pooled_res,idcol="probe")
  fwrite(sum_stats,paste0(exposure,"_PRS_interaction.txt.gz"),sep="\t",row.names = F,quote=F)
}
