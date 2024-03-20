args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  print(length(args))
  stop("At least one argument must be supplied (input file).\n", call. = FALSE)
} else {
  print(paste0("Arg input:  ", args))
  print(length(args))
}



setwd("~/family_variance/")
library(hypred)



# best_parents <- read.table("select_best_parents/best_parents.txt")

marker_list <- readRDS("create_marker_list/phased_marker_list.rds")
marker_info_cM_parent <- readRDS("create_marker_list/marker_info_cM_parent.rds")
indiv_names <- readRDS("create_marker_list/indiv_names.rds")
chrom_names <- readRDS("create_marker_list/chrom_names.rds")

effective_marker_sizes <- c(4, 16, 64, 256, 512, 1024)
h2s <- c(0.9, 0.7, 0.5, 0.3, 0.1)

remove_marker_indices <- readRDS("simulate_phenotypes/remove_marker_indices.rds")
effective_marker_indices <- readRDS("simulate_phenotypes/effective_marker_indices.rds")
lmms <- readRDS("simulate_phenotypes/lmms.rds")
bayesC <- readRDS("simulate_phenotypes/bayesC.rds")
alphas <- readRDS("simulate_phenotypes/alphas.rds")



effective_loci_number <- as.numeric(args[1])
replic <- as.numeric(args[2])
parents_names_all <- as.character(args[3:22])

# effective_loci_number <- as.numeric(best_parents[1, 1])
# replic <- as.numeric(best_parents[1, 2])
# parents_names_all <- as.character(best_parents[1, 3:22])



parents_names <- data.frame(parent1=NA, parent2=NA)
for (i in 1:19){
  for (j in (i+1):20){
    parents_names = rbind(parents_names, 
                          c(parents_names_all[i], parents_names_all[j]))
  }
}
parents_names <- parents_names[-1, ]

set.seed(1)
for (g in 1:nrow(parents_names)){
  offsprings <- matrix(NA, ncol=200)
  offspring_list = vector(mode="list", length=(length(chrom_names) + 1))
  offspring_list[[1]] = parents_names[g, ]
  for (h in 1:length(chrom_names)){
    marker_info = marker_info_cM_parent[marker_info_cM_parent$CHROM==chrom_names[h], c("CHROM", "cM")]
    chrom_length = (max(marker_info$cM) - min(marker_info$cM)) / 100
    marker_numb = nrow(marker_info)
    genome = hypredGenome(num.chr = 1, len.chr = chrom_length, num.snp = marker_numb)
    genome = hypredNewMap(genome, new.map = ((marker_info$cM - min(marker_info$cM)) / 100))
    
    parent1A = marker_list[[h]][, 1, which(indiv_names==parents_names[g, 1])]
    parent1B = marker_list[[h]][, 2, which(indiv_names==parents_names[g, 1])]
    parent2A = marker_list[[h]][, 1, which(indiv_names==parents_names[g, 2])]
    parent2B = marker_list[[h]][, 2, which(indiv_names==parents_names[g, 2])]
    
    gametes1 = matrix(nrow=marker_numb, ncol=200)
    gametes2 = matrix(nrow=marker_numb, ncol=200)
    offspring = matrix(nrow=marker_numb, ncol=200)
    
    for (k in 1:200){
      gametes1[, k] = hypredRecombine(genome, genomeA = parent1A, genomeB = parent1B, mutate = F, 
                                      block = FALSE)
      gametes2[, k] = hypredRecombine(genome, genomeA = parent2A, genomeB = parent2B, mutate = F, 
                                      block = FALSE)
      offspring[, k] = gametes1[, k] + gametes2[, k]
      
      offspring_phased = paste(as.character(gametes1[, k]), as.character(gametes2[, k]), sep="|")
    }
    
    offspring = offspring - 1
    offsprings = rbind(offsprings, offspring)
    
    offspring_list[[h+1]] = offspring_phased
  }
  
  saveRDS(offspring_list,
          paste("simulate_crosses_best_parents/offspring_genotypes/offspring_list_",
                parents_names[g, 1],
                "_",
                parents_names[g, 2],
                ".rds",
                sep=""))
  
  offsprings = offsprings[-1, ]
  offsprings = t(offsprings)
  
  # get family BV, mean, variance
  offsprings_Z = offsprings[, !remove_marker_indices]
  
  offsprings_BV = offsprings_Z %*% 
    alphas[[which(effective_loci_number==effective_marker_sizes)]][[replic]]
  offsprings_BV_fammean = mean(offsprings_BV[, 1])
  offsprings_BV_famvar = var(offsprings_BV[, 1])
  
  offsprings_obj = list(offsprings_BV, offsprings_BV_fammean, offsprings_BV_famvar)
  
  saveRDS(offsprings_obj, paste("simulate_crosses_best_parents/offspring_family_info_true/", 
                                paste(effective_loci_number, replic, 
                                      parents_names[g, 1], parents_names[g, 2], 
                                      sep="_"), 
                                "_.rds", 
                                sep=""))
  
  
  offsprings_Z2 = 
    offsprings_Z[, -effective_marker_indices[[
      which(effective_loci_number==effective_marker_sizes)]][, replic]]
  offsprings_predY_RR = matrix(NA, nrow=200, ncol=5)
  offsprings_predY = matrix(NA, nrow=200, ncol=5)
  for (i in 1:length(lmms)){
    offsprings_predY_RR[, i] = offsprings_Z2 %*% 
      lmms[[i]][[which(effective_loci_number==effective_marker_sizes)]][[replic]]$u
    offsprings_predY[, i] = offsprings_Z2 %*% 
      bayesC[[i]][[which(effective_loci_number==effective_marker_sizes)]][[replic]]$ETA[[1]]$b
  }
  
  offsprings_predY_RR_fammean = sapply(as.data.frame(offsprings_predY_RR), mean)
  offsprings_predY_RR_famvar = sapply(as.data.frame(offsprings_predY_RR), var)
  
  offsprings_predY_fammean = sapply(as.data.frame(offsprings_predY), mean)
  offsprings_predY_famvar = sapply(as.data.frame(offsprings_predY), var)
  
  offsprings_obj = list(offsprings_predY_RR, offsprings_predY_RR_fammean, offsprings_predY_RR_famvar, 
                        offsprings_predY, offsprings_predY_fammean, offsprings_predY_famvar)
  
  saveRDS(offsprings_obj, paste("simulate_crosses_best_parents/offspring_family_info_pred/", 
                                paste(effective_loci_number, replic, 
                                      parents_names[g, 1], parents_names[g, 2], 
                                      sep="_"), 
                                "_.rds", 
                                sep=""))
}
