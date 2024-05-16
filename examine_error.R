setwd("~/family_variance/")
source("codes/phase_correct_funcs.R")



correct_phases <- readRDS("get_phase_error/correct_phases.rds")
geno <- matrix(c(0,0, 0,1, 1,0, 1,1), ncol=4)
off_genotypes <- array(c(0,0, 0,1, 0,2, 1,0, 1,1, 1,2, 2,0, 2,1, 2,2), dim=c(2, 1, 9))
one_par_geno <- array(NA, dim=c(2, 2, 16))
for (i in 1:4){
  for (j in 1:4){
    one_par_geno[, , (i-1)*4 + j] = array(c(geno[, i], geno[, j]), dim=c(2, 2))
  }
  
}



phasing_error <- c(0.005, 0.01, 0.025, 0.05, 0.1, 0.15)

marker_list <- readRDS("create_marker_list/phased_marker_list.rds")
all_ind_names <- readRDS("create_marker_list/indiv_names.rds")
chrom_names <- readRDS("create_marker_list/chrom_names.rds")
# pedigree <- readRDS("create_marker_list_real_data2/pedigree.rds")
cM_dist_markers <- readRDS("get_markers_dist/cM_dist_markers.rds")
marker_info_cM_parent <- readRDS("get_markers_cM/marker_info_cM_parent.rds")



marker_list_error_names <- list.files("introduce_error/", pattern="marker_list_error")
marker_list_error <- vector(mode="list", length=7)
marker_list_error[[1]] <- marker_list
marker_list_error[2:7] <- lapply(marker_list_error_names, FUN=function(x){
  readRDS(paste("introduce_error/", x, sep=""))
})



offspring_names <- list.files("simulate_crosses/", pattern="offspring_list_")
set.seed(1)
offspring_names <- sample(offspring_names, 500)
offspring_temp <- readRDS(paste("simulate_crosses/", offspring_names[1], sep=""))

offspring <- vector(mode="list", length=length(chrom_names))
for (j in 1:length(chrom_names)){
  offspring[[j]] = matrix(NA, nrow=nrow(offspring_temp[[j+1]]), ncol=length(offspring_names))
}

offspring_ped <- data.frame(parent1 = rep(NA, length(offspring_names)), 
                             parent2 = rep(NA, length(offspring_names)))

# set.seed(1)
offspring_random <- sample(1:200, length(offspring_names), replace=T)
head(offspring_random)
# [1]  97 177 156  12  61  56
for (i in 1:length(offspring_names)){
  offspring_temp = readRDS(paste("simulate_crosses/", offspring_names[i], sep=""))
  offspring_ped[i, ] = offspring_temp[[1]]
  for (j in 1:length(chrom_names)){
    offspring[[j]][, i] = sapply(offspring_temp[[j+1]][, offspring_random[i]], FUN=function(x){
      sum(as.numeric(strsplit(x, split="\\|")[[1]]))
    })
    print(c(i, j))
  }
}



par1_names <- offspring_ped$parent1
par2_names <- offspring_ped$parent2

marker_info_cM_parent2 <- marker_info_cM_parent[!is.na(marker_info_cM_parent$cM), ]

for (j in 1:length(cM_dist_markers)){
  marker_names2 = rownames(marker_info_cM_parent2[marker_info_cM_parent2$CHROM == chrom_names[j], ])
  marker_names = rownames(marker_info_cM_parent[marker_info_cM_parent$CHROM == chrom_names[j], ])
  
  for (h in 1:length(cM_dist_markers[[j]])){
    cM_dist_markers[[j]][[h]] = apply(cM_dist_markers[[j]][[h]], 2, FUN=function(x){
      return(match(marker_names2[x], marker_names))
    })
  }
}



score <- vector(mode="list", length=length(marker_list_error))

for (h in 1:length(marker_list_error)){
  score[[h]] = array(NA, dim=c(4, 1000, 5, length(chrom_names), length(offspring_names)))
  
  for (j in 1:length(chrom_names)){
    for (i in 1:length(offspring_names)){
      par1A = marker_list_error[[h]][[j]][, 1, which(all_ind_names==par1_names[i])]
      par1B = marker_list_error[[h]][[j]][, 2, which(all_ind_names==par1_names[i])]
      par2A = marker_list_error[[h]][[j]][, 1, which(all_ind_names==par2_names[i])]
      par2B = marker_list_error[[h]][[j]][, 2, which(all_ind_names==par2_names[i])]
      off = offspring[[j]][, i]
      
      for (l in 1:5){
        dist_markers = cM_dist_markers[[j]][[l]]
        dist_markers = matrix(dist_markers, ncol=2)
        
        if(nrow(dist_markers)==0){next}
        
        for (k in 1:nrow(dist_markers)){
          marker_index = dist_markers[k, ]
          chrom_matx = matrix(c(par1A[marker_index], 
                                par1B[marker_index], 
                                par2A[marker_index], 
                                par2B[marker_index], 
                                off[marker_index]), ncol=5)
          
          keep_marker = apply(chrom_matx, 1, FUN=get_index_cor_off)
          
          if (any(is.na(keep_marker)) | (!all(keep_marker))){next}
          
          which_par1 = which(apply(one_par_geno, c(3), function(x){
            all(x==chrom_matx[, 1:2])}))
          which_par2 = which(apply(one_par_geno, c(3), function(x){
            all(x==chrom_matx[, 3:4])}))
          which_off = which(apply(off_genotypes, c(3), function(x){
            all(x==chrom_matx[, 5])}))
          
          if (all(c(which_par1, which_par2, which_off))){
            score[[h]][1:2, k, l, j, i] = marker_index
            score[[h]][3:4, k, l, j, i] = 
              correct_phases[[which_par1]][[which_par2]][[which_off]]
          }
          
        } # end of k loop
      } # end of l loop
    } # end of i loop
  }
}



saveRDS(score, "examine_error/score.rds")












