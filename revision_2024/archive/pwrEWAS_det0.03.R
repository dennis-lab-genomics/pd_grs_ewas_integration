########### Run the rest of the simulations via command line (lengthy)
library(pwrEWAS)
ExperimentHub::setExperimentHubOption('cache', "/home1/NEURO/schaffner/.local/share/R/ExperimentHub")

# All QC'ed single CpGs
# 800k probes, 27 target DM-CpGs, sample size between 30-210, DB 0.03, FDR 0.05
pwr_800k_27 <- pwrEWAS(minTotSampleSize=30, maxTotSampleSize=210, SampleSizeSteps=30, NcntPer=0.67,
                       targetDelta = 0.03, J = 8e+05, targetDmCpGs=27,
                       tissueType = "Blood adult",
                       detectionLimit = 0.03, DMmethod = "limma",
                       FDRcritVal = 0.05, core = 1, sims = 50)
pwr_800k_27$meanPower
#0.03
#30  0.1461637
#60  0.3377294
#90  0.4325097
#120 0.4698856
#150 0.4996136
#180 0.5413532
#210 0.53667021

# Variability-filtered single CpGs
# 440k probes, 27 target DM-CpGs, sample size between 30-210, DB 0.03, FDR 0.05
pwr_440k_27 <- pwrEWAS(minTotSampleSize=30, maxTotSampleSize=210, SampleSizeSteps=30, NcntPer=0.67,
                       targetDelta = 0.03, J = 4.4e+05, targetDmCpGs=27,
                       tissueType = "Blood adult",
                       detectionLimit = 0.03, DMmethod = "limma",
                       FDRcritVal = 0.05, core = 1, sims = 50)
pwr_440k_27$meanPower
#0.03
#30  0.1716477
#60  0.3666470
#90  0.4419913
#120 0.4689176
#150 0.5067490
#180 0.5157613
#210 0.5557704

# 200k probes, 27 target DM-CpGs, sample size between 30-210, DB 0.03, FDR 0.05
pwr_200k_27 <- pwrEWAS(minTotSampleSize=30, maxTotSampleSize=210, SampleSizeSteps=30, NcntPer=0.67,
                       targetDelta = 0.03, J = 2e+05, targetDmCpGs=27,
                       tissueType = "Blood adult",
                       detectionLimit = 0.03, DMmethod = "limma",
                       FDRcritVal = 0.05, core = 1, sims = 50)
pwr_200k_27$meanPower
#0.03
#30  0.2095202
#60  0.4027408
#90  0.4839355
#120 0.5177796
#150 0.5437318
#180 0.5731346
#210 0.5917916

# 100k probes, 27 target DM-CpGs, sample size between 30-210, DB 0.03, FDR 0.05
pwr_100k_27 <- pwrEWAS(minTotSampleSize=30, maxTotSampleSize=210, SampleSizeSteps=30, NcntPer=0.67,
                       targetDelta = 0.03, J = 1e+05, targetDmCpGs=27,
                       tissueType = "Blood adult",
                       detectionLimit = 0.03, DMmethod = "limma",
                       FDRcritVal = 0.05, core = 1, sims = 50)
pwr_100k_27$meanPower
#0.03
#30  0.2836289
#60  0.5115510
#90  0.6370082
#120 0.6751151
#150 0.7422825
#180 0.7691940
#210 0.7665549

# 50k probes, 27 target DM-CpGs, sample size between 30-210, DB 0.03, FDR 0.05
pwr_50k_27 <- pwrEWAS(minTotSampleSize=30, maxTotSampleSize=210, SampleSizeSteps=30, NcntPer=0.67,
                      targetDelta = 0.03, J = 5e+04, targetDmCpGs=27,
                      tissueType = "Blood adult",
                      detectionLimit = 0.03, DMmethod = "limma",
                      FDRcritVal = 0.05, core = 1, sims = 50)
pwr_50k_27$meanPower
#0.03
#30  0.3659739
#60  0.6008366
#90  0.7607395
#120 0.8139900
#150 0.8140195
#180 0.8234715
#210 0.8598860

