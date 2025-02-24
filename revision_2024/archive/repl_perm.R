repl_perm <- function(sig.cpgs, test.cpgs, n_repl, permutation_number=1000){
   
   #sig.cpgs = data frame of EWAS results
   #discovery_betas = data frame or matrix of betas tested for differential methylation discovery (CpGs as rows, samples as columns)
   #discovery_DBs = named vector of delta betas, corresponding to rownames of discovery_betas
   #validation_betas = data frame or matrix of betas tested for differential methylation validation (CpGs as rows, samples as columns)
   #validation_DBs = named vector of delta betas, corresponding to rownames of validation_betas
   #n_validation = number of CMRs that had the same effect direction in validation dataset
 
   bootstrap_effect <-lapply(1:permutation_number, function(x){
    set.seed(x)
   
    # pull a random set of CpGs from TERRE, the same size as the number of posthoc CpGs
    rnd_cpgs <- sig.cpgs[sample(1:nrow(sig.cpgs), nrow(test.cpgs)),]
    return(nrow(rnd_cpgs[rnd_cpgs$P.Value<=0.05 & abs(rnd_cpgs$delta_M)>1.5,]))
    
    })
   
   bootstrap_effect <- do.call(rbind, bootstrap_effect)
   
   print("Permutation P values for enrichment and depletion")
   #how many iterations of bootstrapping found MORE genes coexpressed than the real data? Divide this by the number of permutations to get a p-value
   enrich_p<- length(which(bootstrap_effect>=n_repl))/permutation_number
   #how many iterations of bootstrapping found MORE genes coexpressed than the real data? Divide this by the number of permutations to get a p-value
   depletion_p<- length(which(bootstrap_effect<=n_repl))/permutation_number
   print(paste("Enrichment: ", enrich_p, "; Depletion ", depletion_p, sep="")) }
