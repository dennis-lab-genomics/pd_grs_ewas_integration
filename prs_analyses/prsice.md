PRSice2 Development of risk score
================

# Step 0: Prepare covariates and input files

``` r
IDs <- fread("~/genotype_qc/TERRE_QC/all_imputed_r2_30_rsid_hard_call.fam")[, .(FID = V1, IID = V2)]
covariate <- fread("../cis_mQTL_analyses/terre_data/covariates_10_methy_PC.txt")[id != "head_trauma_loc"]
covariate <- cbind(IDs, covariate %>% transpose(make.names = "id"))

PD <- fread("../cis_mQTL_analyses/terre_data/covariates_CTP_PD.txt")[id == "PD"]
PD <- cbind(IDs, PD %>% transpose(make.names = "id"))

head(covariate)
```

    ##          FID       IID          PC1       PC2        PC3        PC4        PC5
    ## 1: PAE_015_1 PAE_015_1  -7.03264917  8.057218 -1.1892932   7.680491 -4.3100504
    ## 2: PAE_015_2 PAE_015_2  10.58721134  8.624483 -2.4599118  -0.673395  0.1088861
    ## 3: PAE_015_3 PAE_015_3   8.92736676  7.675730 -1.8317991   2.280433 -1.8181081
    ## 4: PAE_037_1 PAE_037_1  -0.06600985 -4.239611  0.4591169   2.324907 -6.0032639
    ## 5: PAE_037_2 PAE_037_2   2.24229994 -3.128131  6.2919567   6.159199  5.4985993
    ## 6: PAE_037_3 PAE_037_3 -21.97379292  4.993849 -5.0696918 -15.341385 -0.5031651
    ##           PC6       PC7        PC8        PC9        PC10          V3
    ## 1:  3.7019364 -1.705740  1.9124255 -0.5864368 -0.41608154 -0.04530380
    ## 2: -1.8585109 -1.320930 -1.0983885 -0.8055798  1.22520130 -0.07638310
    ## 3:  4.5074378 -2.019942  0.1664676  2.1014179 -2.13724595  0.00291674
    ## 4:  5.8307433  1.913521  1.1002500 -1.4405261 -0.10823954  0.02952700
    ## 5: -2.0491117 -2.030263 -2.0831077 -0.3149775 -0.19753391  0.02969210
    ## 6: -0.4946656 -2.511126 -2.3577727  2.2627508 -0.05729025  0.01579970
    ##             V4          V5      age men
    ## 1:  0.00875784 -0.00251219 56.22998   1
    ## 2:  0.01037420 -0.03894880 60.13689   1
    ## 3:  0.02180150  0.01303290 58.09446   1
    ## 4: -0.00408671  0.00811965 73.64271   0
    ## 5: -0.02416570 -0.02666510 75.90965   0
    ## 6:  0.04695520  0.05616860 75.84668   0

``` r
head(PD)
```

    ##          FID       IID PD
    ## 1: PAE_015_1 PAE_015_1  0
    ## 2: PAE_015_2 PAE_015_2  0
    ## 3: PAE_015_3 PAE_015_3  1
    ## 4: PAE_037_1 PAE_037_1  1
    ## 5: PAE_037_2 PAE_037_2  0
    ## 6: PAE_037_3 PAE_037_3  0

``` r
fwrite(PD, "TERRE.pheno", sep = "\t")
fwrite(covariate, "TERRE.covariate", sep = "\t")
```

``` r
covariate
```

    ##            FID       IID          PC1        PC2        PC3        PC4
    ##   1: PAE_015_1 PAE_015_1  -7.03264917   8.057218 -1.1892932  7.6804906
    ##   2: PAE_015_2 PAE_015_2  10.58721134   8.624483 -2.4599118 -0.6733950
    ##   3: PAE_015_3 PAE_015_3   8.92736676   7.675730 -1.8317991  2.2804332
    ##   4: PAE_037_1 PAE_037_1  -0.06600985  -4.239611  0.4591169  2.3249073
    ##   5: PAE_037_2 PAE_037_2   2.24229994  -3.128131  6.2919567  6.1591991
    ##  ---
    ## 241: PAE_604_3 PAE_604_3  10.92477832   8.135720  5.4413797 -6.9671935
    ## 242: PAE_604_4 PAE_604_4   8.99908327   5.422101  5.4763081 -1.6292101
    ## 243: PAE_623_2 PAE_623_2 -11.56325217 -15.926213  4.5826028 -4.9899042
    ## 244: PAE_652_2 PAE_652_2  -7.20731283 -16.961343 -9.4614013 -5.3316739
    ## 245: PAE_652_3 PAE_652_3 -13.08511696 -18.737425  0.8786771 -0.8225876
    ##             PC5        PC6        PC7        PC8        PC9       PC10
    ##   1: -4.3100504  3.7019364 -1.7057397  1.9124255 -0.5864368 -0.4160815
    ##   2:  0.1088861 -1.8585109 -1.3209301 -1.0983885 -0.8055798  1.2252013
    ##   3: -1.8181081  4.5074378 -2.0199416  0.1664676  2.1014179 -2.1372460
    ##   4: -6.0032639  5.8307433  1.9135206  1.1002500 -1.4405261 -0.1082395
    ##   5:  5.4985993 -2.0491117 -2.0302626 -2.0831077 -0.3149775 -0.1975339
    ##  ---
    ## 241:  9.1301257  5.2281240 -0.2522033  1.5687913 -0.2783727  1.3299875
    ## 242: -5.5858751 -4.2163950  1.1466702 -0.2100942 -1.4269395 -0.5337795
    ## 243: -2.9969255 -0.9507375 -0.6454093  1.9120981  1.5514008  0.5228948
    ## 244:  2.0355945  1.6939029 -1.3403652 -0.3418679  0.6723842  0.2808460
    ## 245:  6.7523450 -4.6575242  1.1036263 -0.2747632 -3.7139110  1.3510379
    ##               V3          V4          V5      age men
    ##   1: -0.04530380  0.00875784 -0.00251219 56.22998   1
    ##   2: -0.07638310  0.01037420 -0.03894880 60.13689   1
    ##   3:  0.00291674  0.02180150  0.01303290 58.09446   1
    ##   4:  0.02952700 -0.00408671  0.00811965 73.64271   0
    ##   5:  0.02969210 -0.02416570 -0.02666510 75.90965   0
    ##  ---
    ## 241:  0.05126360 -0.02965190 -0.06904540 74.59548   1
    ## 242:  0.13972400 -0.11876200 -0.14161800 75.55921   1
    ## 243: -0.03615590 -0.02229510 -0.04420550 50.67214   0
    ## 244:  0.02219500 -0.02890760 -0.00347380 71.79192   0
    ## 245: -0.04783330 -0.01125360  0.01105620 69.76591   0

# Step 1: Run PRSice-2 on Nalls et al 2019 Sumstats

``` bash
Rscript /home1/NEURO/casazza/PRSice.R \
    --prsice /home1/NEURO/casazza/PRSice_linux\
    --base /home1/NEURO/casazza/nalls_PD.QC.gz\
    --base-info INFO:0.8 \
    --base-maf MAF:0.01 \
    --cov TERRE.covariate \
    --binary-target T\
    --beta  \
    --out TERRE_PRSice \
    -q 5\
    --all-score\
    --pheno TERRE.pheno \
    --snp SNP \
    --stat b \
    --pvalue p\
    --target /home1/NEURO/casazza/genotype_qc/TERRE_QC/all_imputed_r2_30_rsid_hard_call \
    --thread 32

Rscript /home1/NEURO/casazza/PRSice.R \
    --prsice /home1/NEURO/casazza/PRSice_linux\
    --base /home1/NEURO/casazza/nalls_PD.QC.gz\
    --base-info INFO:0.8 \
    --base-maf MAF:0.01 \
    --cov TERRE.covariate \
    --binary-target T\
    --beta  \
    --out TERRE_PRSice_nalls_male \
    -q 5\
    --all-score\
    --pheno TERRE.pheno \
    --snp SNP \
    --stat b \
    --pvalue p\
    --target /home1/NEURO/casazza/genotype_qc/TERRE_QC/all_imputed_r2_30_rsid_hard_call_male \
    --thread 32

Rscript /home1/NEURO/casazza/PRSice.R \
    --prsice /home1/NEURO/casazza/PRSice_linux\
    --base /home1/NEURO/casazza/nalls_PD.QC.gz\
    --base-info INFO:0.8 \
    --base-maf MAF:0.01 \
    --cov TERRE.covariate \
    --binary-target T\
    --beta  \
    --out TERRE_PRSice_nalls_female \
    -q 5\
    --all-score\
    --pheno TERRE.pheno \
    --snp SNP \
    --stat b \
    --pvalue p\
    --target /home1/NEURO/casazza/genotype_qc/TERRE_QC/all_imputed_r2_30_rsid_hard_call_female \
    --thread 32
```

# Step 2: Evaluate output

``` r
include_graphics("prsice_images/TERRE_PRSice_BARPLOT_2021-11-02.png")
```

<img src="prsice_images/TERRE_PRSice_BARPLOT_2021-11-02.png" width="400px" />

``` r
include_graphics("prsice_images/TERRE_PRSice_HIGH-RES_PLOT_2021-11-02.png")
```

<img src="prsice_images/TERRE_PRSice_HIGH-RES_PLOT_2021-11-02.png" width="400px" />

``` r
include_graphics("prsice_images/TERRE_PRSice_QUANTILES_PLOT_2021-11-02.png")
```

<img src="prsice_images/TERRE_PRSice_QUANTILES_PLOT_2021-11-02.png" width="400px" />

``` r
include_graphics("prsice_images/TERRE_PRSice_nalls_male_BARPLOT_2022-01-13.png")
```

<img src="prsice_images/TERRE_PRSice_nalls_male_BARPLOT_2022-01-13.png" width="400px" />

``` r
include_graphics("prsice_images/TERRE_PRSice_nalls_male_HIGH-RES_PLOT_2022-01-13.png")
```

<img src="prsice_images/TERRE_PRSice_nalls_male_HIGH-RES_PLOT_2022-01-13.png" width="400px" />

``` r
include_graphics("prsice_images/TERRE_PRSice_nalls_male_QUANTILES_PLOT_2022-01-13.png")
```

<img src="prsice_images/TERRE_PRSice_nalls_male_QUANTILES_PLOT_2022-01-13.png" width="400px" />

``` r
include_graphics("prsice_images/TERRE_PRSice_nalls_female_BARPLOT_2022-01-13.png")
```

<img src="prsice_images/TERRE_PRSice_nalls_female_BARPLOT_2022-01-13.png" width="400px" />

``` r
include_graphics("prsice_images/TERRE_PRSice_nalls_female_HIGH-RES_PLOT_2022-01-13.png")
```

<img src="prsice_images/TERRE_PRSice_nalls_female_HIGH-RES_PLOT_2022-01-13.png" width="400px" />

``` r
include_graphics("prsice_images/TERRE_PRSice_nalls_female_QUANTILES_PLOT_2022-01-13.png")
```

<img src="prsice_images/TERRE_PRSice_nalls_female_QUANTILES_PLOT_2022-01-13.png" width="400px" />
\#\# Plotting PRSice Data on my own

``` r
library(ggnewscale)
prsice_male_meta <- fread("prsice_nalls_male_data/TERRE_PRSice_nalls_male.prsice")
ggplot(prsice_male_meta[Threshold <= 0.5], aes(Threshold, Num_SNP, color = -log10(P))) +
  geom_point() +
  scale_y_continuous(breaks = c(seq(0, 1e5, 2.5e4), seq(2e5, 6e5, 1e5))) +
  theme_minimal()
```

![](prsice_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
prsice_female_meta <- fread("prsice_nalls_female_data/TERRE_PRSice_nalls_female.prsice")
ggplot(prsice_female_meta[Threshold <= 0.5], aes(Threshold, Num_SNP, color = -log10(P))) +
  geom_point() +
  scale_y_continuous(breaks = c(seq(0, 1e5, 2.5e4), seq(2e5, 6e5, 1e5))) +
  theme_minimal()
```

![](prsice_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
prsice_meta <- fread("prsice_data/TERRE_PRSice.prsice")
ggplot(mapping = aes(Threshold, R2, color = -log10(P))) +
  geom_point(data = prsice_female_meta, size = 1) +
  scale_color_gradient(low = "darkred", high = "pink") +
  new_scale_color() +
  geom_point(data = prsice_male_meta, size = 1, aes(color = -log10(P))) +
  scale_color_gradient(low = "darkblue", high = "lightblue2") +
  labs(y = bquote("R"^2), x = "GWAS P-Value Threshold", color = bquote("-log"["10"] ~ "(P)")) +
  theme_minimal()
```

![](prsice_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
ggplot(mapping = aes(Threshold, R2, color = -log10(P))) +
  geom_point(data = prsice_female_meta, size = 1) +
  scale_color_gradient(low = "darkred", high = "pink") +
  new_scale_color() +
  geom_point(data = prsice_male_meta, size = 1, aes(color = -log10(P))) +
  scale_color_gradient(low = "darkblue", high = "lightblue2") +
  new_scale_color() +
  geom_point(data = prsice_meta, size = 1, aes(color = -log10(P))) +
  labs(y = bquote("R"^2), x = "GWAS P-Value Threshold", color = bquote("-log"["10"] ~ "(P)")) +
  theme_minimal()
```

![](prsice_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

# Step 3 run linear model at different thresholds for SNP inclusion

``` r
load("/home1/NEURO/SHARE_DECIPHER/processed_DNAm_data/TERRE_processed_2021/betas_combat.RData") # betas_sub
# Assign genotyping ID to data
original_covars <- fread("/home1/NEURO/SHARE_DECIPHER/terre_meta_master.csv")[, .(patient, IID = gsub("_PAE.*", "", IID))]
colnames(betas_combat) <- original_covars$IID[match(colnames(betas_combat), original_covars$patient)]
betas_combat <- betas_combat[, colnames(betas_combat) %in% covariate$IID]
```

Let’s check how the data looks for the first 5 subjects:

``` r
ggplot(betas_combat[, 1:5] %>% as.data.table(keep.rownames = T) %>% melt(id.vars = "rn", value.name = "betas", variable.name = "subject"), aes(betas, color = subject)) +
  geom_density()
```

![](prsice_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Match DNA, PRS, and metadata

``` r
prsice_all <- fread("prsice_data/TERRE_PRSice.all_score")[match(colnames(betas_combat), IID), .(FID, IID, `Pt_0.0219001`, `Pt_5e-08`, `Pt_5.005e-05`, `Pt_0.00010005`, `Pt_0.00100005`, `Pt_0.0101501`, `Pt_0.1`, `Pt_0.2`, `Pt_0.3`, `Pt_0.4`, `Pt_0.5`, `Pt_1`)]
covariate <- covariate[match(colnames(betas_combat), IID)]
all(covariate$IID == colnames(betas_combat))
```

    ## [1] TRUE

``` r
all(covariate$IID == prsice_all$IID)
```

    ## [1] TRUE

``` r
covariate_male <- covariate[men == 1] %>% select(-men)
betas_male <- betas_combat[, covariate_male$IID]
prsice_male_all <- prsice_all[match(colnames(betas_male), IID), .(FID, IID, `Pt_0.0219001`, `Pt_5e-08`, `Pt_5.005e-05`, `Pt_0.00010005`, `Pt_0.00100005`, `Pt_0.0101501`, `Pt_0.1`, `Pt_0.2`, `Pt_0.3`, `Pt_0.4`, `Pt_0.5`, `Pt_1`)]
covariate_male <- covariate_male[match(colnames(betas_male), IID)]
all(covariate_male$IID == colnames(betas_male))
```

    ## [1] TRUE

``` r
all(covariate_male$IID == prsice_male_all$IID)
```

    ## [1] TRUE

``` r
covariate_female <- covariate[men == 0] %>% select(-men)
betas_female <- betas_combat[, covariate_female$IID]
prsice_female_all <- prsice_all[match(colnames(betas_female), IID), .(FID, IID, `Pt_0.0219001`, `Pt_5e-08`, `Pt_5.005e-05`, `Pt_0.00010005`, `Pt_0.00100005`, `Pt_0.0101501`, `Pt_0.1`, `Pt_0.2`, `Pt_0.3`, `Pt_0.4`, `Pt_0.5`, `Pt_1`)]
covariate_female <- covariate_female[match(colnames(betas_female), IID)]
all(covariate_female$IID == colnames(betas_female))
```

    ## [1] TRUE

``` r
all(covariate_female$IID == prsice_female_all$IID)
```

    ## [1] TRUE

``` r
covariate_male <- covariate[men == 1] %>% select(-men)
betas_male <- betas_combat[, covariate_male$IID]
prsice_male_nalls_all <- fread("prsice_nalls_male_data/TERRE_PRSice_nalls_male.all_score")[match(colnames(betas_male), IID), .(FID, IID, `Pt_0.0219001`, `Pt_5e-08`, `Pt_5.005e-05`, `Pt_0.00010005`, `Pt_0.00100005`, `Pt_0.0101501`, `Pt_0.1`, `Pt_0.2`, `Pt_0.3`, `Pt_0.4`, `Pt_0.5`, `Pt_1`)]
covariate_male <- covariate_male[match(colnames(betas_male), IID)]
all(covariate_male$IID == colnames(betas_male))
```

    ## [1] TRUE

``` r
all(covariate_male$IID == prsice_male_nalls_all$IID)
```

    ## [1] TRUE

``` r
covariate_female <- covariate[men == 0] %>% select(-men)
betas_female <- betas_combat[, covariate_female$IID]
prsice_female_nalls_all <- fread("prsice_nalls_female_data/TERRE_PRSice_nalls_female.all_score")[match(colnames(betas_female), IID), .(FID, IID, `Pt_0.0219001`, `Pt_5e-08`, `Pt_5.005e-05`, `Pt_0.00010005`, `Pt_0.00100005`, `Pt_0.0101501`, `Pt_0.1`, `Pt_0.2`, `Pt_0.3`, `Pt_0.4`, `Pt_0.5`, `Pt_1`)]
covariate_female <- covariate_female[match(colnames(betas_female), IID)]
all(covariate_female$IID == colnames(betas_female))
```

    ## [1] TRUE

``` r
all(covariate_female$IID == prsice_female_nalls_all$IID)
```

    ## [1] TRUE

### Run limma

``` r
mvalues <- lumi::beta2m(betas_combat)
```

    ## Setting options('download.file.method.GEOquery'='auto')

    ## Setting options('GEOquery.inmemory.gpl'=FALSE)

    ## No methods found in package 'RSQLite' for request: 'dbListFields' when loading 'lumi'

``` r
prs_mat <- prsice_all[, -c(1, 2)]
cov_mat <- covariate[, -c(1, 2)]

mvalues_male <- lumi::beta2m(betas_male)
prs_mat_male <- prsice_male_all[, -c(1, 2)]
cov_mat_male <- covariate_male[, -c(1, 2)]

mvalues_female <- lumi::beta2m(betas_female)
prs_mat_female <- prsice_female_all[, -c(1, 2)]
cov_mat_female <- covariate_female[, -c(1, 2)]

prs_mat_nalls_male <- prsice_male_nalls_all[, -c(1, 2)]
prs_mat_nalls_female <- prsice_female_nalls_all[, -c(1, 2)]
```

``` r
registerDoParallel(ncol(prs_mat) / 2)
hits <- foreach(prs_thresh = colnames(prs_mat)) %dopar% {
  design_prs <- model.matrix(~., data = cbind(prs_mat[, ..prs_thresh], cov_mat))
  prs_fit <- lmFit(mvalues, design_prs)
  prs_fit <- eBayes(prs_fit)
  topTable(prs_fit, coef = 2, adjust.method = "bonf", p.value = 0.05, number = Inf, genelist = rownames(mvalues))
}
names(hits) <- colnames(prs_mat)
hits_by_thresh_bonf <- rbindlist(hits, idcol = "threshold", fill = TRUE)

registerDoParallel(ncol(prs_mat_male) / 2)
hits_male <- foreach(prs_thresh = colnames(prs_mat_male)) %dopar% {
  design_prs_male <- model.matrix(~., data = cbind(prs_mat_male[, ..prs_thresh], cov_mat_male))
  prs_fit_male <- lmFit(mvalues_male, design_prs_male)
  prs_fit_male <- eBayes(prs_fit_male)
  topTable(prs_fit_male, coef = 2, adjust.method = "bonf", p.value = 0.05, number = Inf, genelist = rownames(mvalues_male))
}
names(hits_male) <- colnames(prs_mat_male)
hits_by_thresh_bonf_male <- rbindlist(hits_male, idcol = "threshold", fill = TRUE)

registerDoParallel(ncol(prs_mat_female) / 2)
hits_female <- foreach(prs_thresh = colnames(prs_mat_female)) %dopar% {
  design_prs_female <- model.matrix(~., data = cbind(prs_mat_female[, ..prs_thresh], cov_mat_female))
  prs_fit_female <- lmFit(mvalues_female, design_prs_female)
  prs_fit_female <- eBayes(prs_fit_female)
  topTable(prs_fit_female, coef = 2, adjust.method = "bonf", p.value = 0.05, number = Inf, genelist = rownames(mvalues_female))
}
names(hits_female) <- colnames(prs_mat_female)
hits_by_thresh_bonf_female <- rbindlist(hits_female, idcol = "threshold", fill = TRUE)

registerDoParallel(ncol(prs_mat_male) / 2)
hits_nalls_male <- foreach(prs_thresh = colnames(prs_mat_nalls_male)) %dopar% {
  design_prs_nalls_male <- model.matrix(~., data = cbind(prs_mat_nalls_male[, ..prs_thresh], cov_mat_male))
  prs_fit_nalls_male <- lmFit(mvalues_male, design_prs_nalls_male)
  prs_fit_nalls_male <- eBayes(prs_fit_nalls_male)
  topTable(prs_fit_nalls_male, coef = 2, adjust.method = "bonf", p.value = 0.05, number = Inf, genelist = rownames(mvalues_male))
}
names(hits_nalls_male) <- colnames(prs_mat_nalls_male)
hits_by_thresh_bonf_nalls_male <- rbindlist(hits_nalls_male, idcol = "threshold", fill = TRUE)

registerDoParallel(ncol(prs_mat_female) / 2)
hits_nalls_female <- foreach(prs_thresh = colnames(prs_mat_nalls_female)) %dopar% {
  design_prs_nalls_female <- model.matrix(~., data = cbind(prs_mat_nalls_female[, ..prs_thresh], cov_mat_female))
  prs_fit_nalls_female <- lmFit(mvalues_female, design_prs_nalls_female)
  prs_fit_nalls_female <- eBayes(prs_fit_nalls_female)
  topTable(prs_fit_nalls_female, coef = 2, adjust.method = "bonf", p.value = 0.05, number = Inf, genelist = rownames(mvalues_female))
}
names(hits_nalls_female) <- colnames(prs_mat_nalls_female)
hits_by_thresh_bonf_nalls_female <- rbindlist(hits_nalls_female, idcol = "threshold", fill = TRUE)
```

``` r
ggplot(hits_by_thresh_bonf[, .(hits = .N), by = threshold] %>% mutate(threshold = recode_factor(threshold, `Pt_0.0219001` = "0.0219", `Pt_5e-08` = "5e-8", `Pt_5.005e-05` = "5e-5", `Pt_0.00010005` = "1e-4", `Pt_0.00100005` = "1e-3", `Pt_0.0101501` = "1e-2", `Pt_0.1` = "0.1", `Pt_0.2` = "0.2", `Pt_0.3` = "0.3", `Pt_0.4` = "0.4", `Pt_0.5` = "0.5", `Pt_1` = "1.0")), aes(threshold, hits, label = hits)) +
  geom_text(vjust = -0.25) +
  geom_col() +
  labs(x = "GWAS P Value Threshold", y = "EWAS Hits") +
  theme_minimal()
```

![](prsice_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
ggplot(hits_by_thresh_bonf_male[, .(hits = .N), by = threshold] %>% mutate(threshold = recode_factor(threshold, `Pt_0.0219001` = "0.0219", `Pt_5e-08` = "5e-8", `Pt_5.005e-05` = "5e-5", `Pt_0.00010005` = "1e-4", `Pt_0.00100005` = "1e-3", `Pt_0.0101501` = "1e-2", `Pt_0.1` = "0.1", `Pt_0.2` = "0.2", `Pt_0.3` = "0.3", `Pt_0.4` = "0.4", `Pt_0.5` = "0.5", `Pt_1` = "1.0")), aes(threshold, hits, label = hits)) +
  geom_text(vjust = -0.25) +
  geom_col() +
  labs(x = "GWAS P Value Threshold", y = "EWAS Hits") +
  theme_minimal()
```

![](prsice_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
ggplot(hits_by_thresh_bonf_female[, .(hits = .N), by = threshold] %>% mutate(threshold = recode_factor(threshold, `Pt_0.0219001` = "0.0219", `Pt_5e-08` = "5e-8", `Pt_5.005e-05` = "5e-5", `Pt_0.00010005` = "1e-4", `Pt_0.00100005` = "1e-3", `Pt_0.0101501` = "1e-2", `Pt_0.1` = "0.1", `Pt_0.2` = "0.2", `Pt_0.3` = "0.3", `Pt_0.4` = "0.4", `Pt_0.5` = "0.5", `Pt_1` = "1.0")), aes(threshold, hits, label = hits)) +
  geom_text(vjust = -0.25) +
  geom_col() +
  labs(x = "GWAS P Value Threshold", y = "EWAS Hits") +
  theme_minimal()
```

![](prsice_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

``` r
hits_by_thresh_bonf_nalls_male[, .(hits = .N), by = threshold]
```

    ##        threshold hits
    ## 1:      Pt_5e-08   34
    ## 2:  Pt_5.005e-05   10
    ## 3: Pt_0.00010005    6

``` r
ggplot(hits_by_thresh_bonf_nalls_male[, .(hits = .N), by = threshold] %>% mutate(threshold = recode_factor(threshold, `Pt_0.0219001` = "0.0219", `Pt_5e-08` = "5e-8", `Pt_5.005e-05` = "5e-5", `Pt_0.00010005` = "1e-4", `Pt_0.00100005` = "1e-3", `Pt_0.0101501` = "1e-2", `Pt_0.1` = "0.1", `Pt_0.2` = "0.2", `Pt_0.3` = "0.3", `Pt_0.4` = "0.4", `Pt_0.5` = "0.5", `Pt_1` = "1.0")), aes(threshold, hits, label = hits)) +
  geom_text(vjust = -0.25) +
  geom_col() +
  labs(x = "GWAS P Value Threshold", y = "EWAS Hits") +
  theme_minimal()
```

![](prsice_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->

``` r
hits_by_thresh_bonf_nalls_female[, .(hits = .N), by = threshold]
```

    ##    threshold hits
    ## 1:  Pt_5e-08   18

``` r
ggplot(hits_by_thresh_bonf_nalls_female[, .(hits = .N), by = threshold] %>% mutate(threshold = recode_factor(threshold, `Pt_0.0219001` = "0.0219", `Pt_5e-08` = "5e-8", `Pt_5.005e-05` = "5e-5", `Pt_0.00010005` = "1e-4", `Pt_0.00100005` = "1e-3", `Pt_0.0101501` = "1e-2", `Pt_0.1` = "0.1", `Pt_0.2` = "0.2", `Pt_0.3` = "0.3", `Pt_0.4` = "0.4", `Pt_0.5` = "0.5", `Pt_1` = "1.0")), aes(threshold, hits, label = hits)) +
  geom_text(vjust = -0.25) +
  geom_col() +
  labs(x = "GWAS P Value Threshold", y = "EWAS Hits") +
  theme_minimal()
```

![](prsice_files/figure-gfm/unnamed-chunk-11-5.png)<!-- -->

``` r
display_venn <- function(x, ...) {
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
display_venn(list(both = hits_by_thresh_bonf[threshold == "Pt_5e-08"]$ID, male = hits_by_thresh_bonf_male[threshold == "Pt_5e-08"]$ID, female = hits_by_thresh_bonf_female[threshold == "Pt_5e-08"]$ID), fill = c("gray", "blue", "pink"))
```

    ## Loading required package: grid

    ## Loading required package: futile.logger

![](prsice_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
display_venn(list(both = hits_by_thresh_bonf[threshold == "Pt_5e-08"]$ID, male = hits_by_thresh_bonf_nalls_male[threshold == "Pt_5e-08"]$ID, female = hits_by_thresh_bonf_nalls_female[threshold == "Pt_5e-08"]$ID), fill = c("gray", "blue", "pink"))
```

![](prsice_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

``` r
library(missMethyl)
gometh(hits_by_thresh_bonf[threshold == "Pt_5e-08"]$ID, array.type = "EPIC", all.cpg = rownames(mvalues))
```

``` r
top_design_prs <- model.matrix(~., data = cbind(prs_mat[, `Pt_5e-08`], cov_mat))
top_prs_fit <- lmFit(mvalues, top_design_prs)
top_prs_fit <- eBayes(top_prs_fit)
top_prs_hits <- topTable(top_prs_fit, coef = 2, adjust.method = "bonf", number = Inf, genelist = rownames(mvalues))

top_design_prs_gxe <- model.matrix(~., data = cbind(prs_mat[, `Pt_5e-08`], cov_mat))
top_prs_fit_gxe <- lmFit(mvalues, top_design_prs_gxe)
top_prs_fit_gxe <- eBayes(top_prs_fit_gxe)
top_prs_hits_gxe <- topTable(top_prs_fit_gxe, coef = 2, adjust.method = "bonf", number = Inf, genelist = rownames(mvalues))
```

``` r
manifest <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other %>%
  as.data.frame() %>%
  rownames_to_column(var = "name")
prs_annot <- data.table(top_prs_hits)[manifest, gene := gsub(";.*", "", UCSC_RefGene_Name), on = c(ID = "name")]
ggplot(prs_annot, aes(logFC, -log10(P.Value))) +
  geom_point() +
  geom_point(
    data = subset(prs_annot, adj.P.Val < 0.05 & abs(logFC) > 0.03),
    color = "green",
    mapping = aes(logFC, -log10(P.Value))
  ) +
  geom_point(
    data = subset(prs_annot, ID %in% hits_by_thresh_bonf_nalls_male[threshold == "Pt_5e-08"]$ID),
    color = "lightblue2",
    mapping = aes(logFC, -log10(P.Value))
  ) +
  geom_point(
    data = subset(prs_annot, ID %in% hits_by_thresh_bonf_nalls_female[threshold == "Pt_5e-08"]$ID),
    color = "pink",
    mapping = aes(logFC, -log10(P.Value))
  ) +
  geom_hline(
    linetype = "dashed",
    yintercept = min(-log10(prs_annot$P.Value[prs_annot$adj.P.Val < 0.05]))
  ) +
  geom_vline(linetype = "dashed", xintercept = 0.03) +
  geom_vline(linetype = "dashed", xintercept = -0.03) +
  geom_text_repel(
    data = prs_annot %>% filter(abs(logFC) > 0.03 & adj.P.Val < 0.05, P.Value < 1e-10),
    color = "dodgerblue",
    mapping = aes(logFC, -log10(P.Value), label = ifelse(gene != "", gene, ID)),
    size = 3,
    max.overlaps = 20
  ) +
  labs(y = bquote("log"[10] ~ "(P)"), x = quote(Delta ~ beta ~ Methylation)) +
  theme_minimal()
```

![](prsice_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

# Manhattan plot vs GWAS manhattan plot

``` r
library(qqman)
cpg_pos <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations %>%
  as.data.frame() %>%
  rownames_to_column(var = "name")
copy_annot <- prs_annot[cpg_pos, on = c(ID = "name")] %>%
  filter(chr != "chrX|chrY") %>%
  mutate(chr = as.numeric(gsub("chr", "", chr)))
manhattan(copy_annot[, .(SNP = gene, CHR = chr, BP = pos, P = P.Value)],
  annotatePval = max(prs_annot$P.Value[prs_annot$adj.P.Val < 0.05]),
  annotateTop = TRUE,
  genomewideline = min(-log10(prs_annot$P.Value[prs_annot$adj.P.Val < 0.05])),
  suggestiveline = FALSE
)
manhattan(copy_annot[chr == 17, .(SNP = gene, CHR = chr, BP = pos, P = P.Value)],
  annotatePval = max(prs_annot$P.Value[prs_annot$adj.P.Val < 0.05]),
  annotateTop = FALSE,
  genomewideline = min(-log10(prs_annot$P.Value[prs_annot$adj.P.Val < 0.05])),
  suggestiveline = FALSE
)
```

``` r
GWAS <- fread("~/nalls_PD.QC.gz")
```

``` r
manhattan(GWAS,
  annotatePval = 5e-8,
  bp = "POS",
  p = "p",
  suggestiveline = FALSE
)
```
