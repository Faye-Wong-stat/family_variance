setwd("~/family_variance/")
library(hypred)



# load the marker list
marker_list <- readRDS("create_marker_list/phased_marker_list.rds")
marker_info_cM_parent <- readRDS("create_marker_list/marker_info_cM_parent.rds")
indiv_names <- readRDS("create_marker_list/indiv_names.rds")
chrom_names <- readRDS("create_marker_list/chrom_names.rds")

dim(marker_info_cM_parent)
# [1] 34911    10
length(indiv_names)
# [1] 1007

set.seed(1)
for (i in 1:length(indiv_names)){
  if(i %% 50 == 0){
    print(i)
  }
  
  gamete_list = vector(mode="list", length=length(chrom_names)+1)
  gamete_list[[1]] = indiv_names[i]
  for (j in 1:length(chrom_names)){
    marker_info = marker_info_cM_parent[marker_info_cM_parent$CHROM==chrom_names[j], c("CHROM", "cM")]
    chrom_length = (max(marker_info$cM) - min(marker_info$cM)) / 100
    marker_numb = nrow(marker_info)
    genome = hypredGenome(num.chr = 1, len.chr = chrom_length, num.snp = marker_numb)
    genome = hypredNewMap(genome, new.map = ((marker_info$cM - min(marker_info$cM)) / 100))
    
    parentA = marker_list[[j]][, 1, i]
    parentB = marker_list[[j]][, 2, i]
    gametes = matrix(nrow=marker_numb, ncol=200)
    
    for (k in 1:200){
      gametes[, k] = hypredRecombine(genome, genomeA = parentA, genomeB = parentB, mutate = F,
                                     block = FALSE)
    }
    
    gamete_list[[j+1]] = gametes
  }
  
  saveRDS(gamete_list, paste("simulate_gametes/", indiv_names[i], ".rds", sep=""))
}



file_names <- list.files("simulate_gametes/", pattern=".rds")
file_names <- paste("simulate_gametes/", file_names, sep="")
write(file_names, "simulate_gametes/file_names.txt")


