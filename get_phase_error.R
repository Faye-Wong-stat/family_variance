setwd("~/family_variance/")
source("codes/phase_correct_funcs.R")
library(vcfR)
library(abind)



geno <- matrix(c(0,0, 0,1, 1,0, 1,1), ncol=4)
off_geno <- array(c(0,0, 0,1, 0,2, 1,0, 1,1, 1,2, 2,0, 2,1, 2,2), dim=c(2, 1, 9))

one_par_geno <- c()
for (i in 1:4){
  for (j in 1:4){
    one_par_geno[[(i-1)*4 + j]] = array(c(geno[, i], geno[, j]), dim=c(2, 2))
  }
}
names(one_par_geno) <- 1:16
two_pars_geno <- c()
for (i in 1:16){
  for (j in 1:16){
    two_pars_geno[[(i-1)*16+j]] = cbind(one_par_geno[[i]], one_par_geno[[j]])
  }
}
names(two_pars_geno) <- 1:256
pars_off_corrt_arry <- c()
for (i in 1:256){
  four_off = gen_off_geno(two_pars_geno[[i]])
  parent = array(two_pars_geno[[i]], dim=c(2, 4, dim(four_off)[3]))
  par_off = abind(parent, four_off, along=2)
  pars_off_corrt_arry = abind(pars_off_corrt_arry, par_off, along=3)
}

correct_phases <- c()
for (i in 1:16){
  correct_phases[[i]] = vector(mode="list", length=16)
  for (j in 1:16){
    correct_phases[[i]][[j]] = vector(mode="list", length=9)
    for (k in 1:9){
      correct_phases[[i]][[j]][[k]] = c(NA, NA)
    }
    names(correct_phases[[i]][[j]]) = 1:9
  }
  names(correct_phases[[i]]) = 1:16
}
names(correct_phases) = 1:16
for (i in 1:16){
  for (j in 1:16){
    for (k in 1:9){
      par1_unphased = convert_geno_mtx(one_par_geno[[i]])
      par2_unphased = convert_geno_mtx(one_par_geno[[j]])
      
      check_par2 = T
      
      if (is_11(par1_unphased)){
        par1_phased = one_par_geno[[i]]
        par2_phased = one_par_geno[[j]]
        offspring_unphased = off_geno[, , k]
        par1_off = cbind(par1_phased, offspring_unphased)
        
        par1_tilde = twist_11(par1_phased)
        par1_tilde_off = cbind(par1_tilde, offspring_unphased)
        
        par1_off_logc = apply(pars_off_corrt_arry[, c(1:2, 5), ], c(3), function(x){
          all(x==par1_off)})
        par1_tilde_off_logc = apply(pars_off_corrt_arry[, c(1:2, 5), ], c(3), function(x){
          all(x==par1_tilde_off)})
        
        if (any(par1_off_logc) & 
            !any(par1_tilde_off_logc)){
          # score correct
          correct_phases[[i]][[j]][[k]][1] = T
        } else if (!any(par1_off_logc) & 
                   any(par1_tilde_off_logc)){
          # score incorrect
          correct_phases[[i]][[j]][[k]][1] = F
        } else if (any(par1_off_logc) & 
                   any(par1_tilde_off_logc)){
          # look at the second parent 
          
          par1_par2_off = cbind(par1_phased, par2_phased, offspring_unphased)
          par1_tilde_par2_off = cbind(par1_tilde, par2_phased, offspring_unphased)
          par1_par2_off_logc = apply(pars_off_corrt_arry, c(3), function(x){
            all(x==par1_par2_off)})
          par1_tilde_par2_off_logc = apply(pars_off_corrt_arry, c(3), function(x){
            all(x==par1_tilde_par2_off)})
          if (any(par1_par2_off_logc) & 
              !any(par1_tilde_par2_off_logc)){
            # score correct
            correct_phases[[i]][[j]][[k]][1] = T
          } else if (!any(par1_par2_off_logc) & 
                     any(par1_tilde_par2_off_logc)){
            # score incorrect
            correct_phases[[i]][[j]][[k]][1] = F
          } else if (!any(par1_par2_off_logc) & 
                     !any(par1_tilde_par2_off_logc)){
            print(c(i, j, k))
            print("either type1 or type2 1/1 parent1 along with parent2 works, check genotype")
          }
          check_par2 = F
          
        } else {
          print(c(i, j, k))
          print("either type1 or type2 1/1 parent1 works, check genotype")
        }
      }
      
      if (is_11(par2_unphased) & check_par2){
        par1_phased = one_par_geno[[i]]
        par2_phased = one_par_geno[[j]]
        offspring_unphased = off_geno[, , k]
        par2_off = cbind(par2_phased, offspring_unphased)
        
        par2_tilde = twist_11(par2_phased)
        par2_tilde_off = cbind(par2_tilde, offspring_unphased)
        
        par2_off_logc = apply(pars_off_corrt_arry[, 3:5, ], c(3), function(x){
          all(x==par2_off)})
        par2_tilde_off_logc = apply(pars_off_corrt_arry[, 3:5, ], c(3), function(x){
          all(x==par2_tilde_off)})
        
        if (any(par2_off_logc) & 
            !any(par2_tilde_off_logc)){
          # score correct 
          correct_phases[[i]][[j]][[k]][2] = T
        } else if (!any(par2_off_logc) & 
                   any(par2_tilde_off_logc)){
          # score incorrect
          correct_phases[[i]][[j]][[k]][2] = F
        } else if (any(par2_off_logc) & 
                   any(par2_tilde_off_logc)){
          # look at the second parent 
          
          par1_par2_off = cbind(par1_phased, par2_phased, offspring_unphased)
          par1_par2_tilde_off = cbind(par1_phased, par2_tilde, offspring_unphased)
          par1_par2_off_logc = apply(pars_off_corrt_arry, c(3), function(x){
            all(x==par1_par2_off)})
          par1_par2_tilde_off_logc = apply(pars_off_corrt_arry, c(3), function(x){
            all(x==par1_par2_tilde_off)})
          if (any(par1_par2_off_logc) & 
              !any(par1_par2_tilde_off_logc)){
            # score correct
            correct_phases[[i]][[j]][[k]][2] = T
          } else if (!any(par1_par2_off_logc) & 
                     any(par1_par2_tilde_off_logc)){
            # score incorrect
            correct_phases[[i]][[j]][[k]][2] = F
          } else if (!any(par1_par2_off_logc) & 
                     !any(par1_par2_tilde_off_logc)){
            print(c(i, j, k))
            print("either type1 or type2 1/1 parent2 along with parent1 works, check genotype")
          }
          
        } else {
          print(c(i, j, k))
          print("either type1 or type2 1/1 parent2 works, check genotype")
        }
      }
      
      
    }
  }
}

saveRDS(correct_phases, "get_phase_error/correct_phases.rds")



marker_matrix <- readRDS("generate_vcffiles/marker_matrix.rds")
marker_info_cM_parent <- readRDS("get_markers_cM/marker_info_cM_parent.rds")
marker_info_cM_valid <- readRDS("get_markers_cM/marker_info_cM_valid.rds")
# marker_info_cM_parent <- read.csv("cleaned_data/marker_info_cM_parent.csv", row.names=1)
# marker_info_cM_valid <- read.csv("cleaned_data/marker_info_cM_valid.csv", row.names=1)
pedigree_valid <- readRDS("generate_vcffiles/pedigree_valid.rds")
# index_low_geno_error <- readRDS("get_errors/index_off_low_error.rds")
# marker_matrix <- read.csv("cleaned_data/marker_matrix.csv", row.names=1, check.names=F)

phased_marker_parent <- read.vcfR("phased_data/phased_parent.vcf")
phased_marker_info_parent <- as.data.frame(phased_marker_parent@fix)
phased_marker_matrix_parent <- as.data.frame(phased_marker_parent@gt)
phased_marker_info_parent$FORMAT <- phased_marker_matrix_parent$FORMAT
phased_marker_matrix_parent <- phased_marker_matrix_parent[, -1]
rownames(phased_marker_matrix_parent) <- phased_marker_info_parent$ID

cM_dist_markers <- readRDS("get_markers_dist/cM_dist_markers.rds")

par1_names <- pedigree_valid$P1
par2_names <- pedigree_valid$P2
off_names <- pedigree_valid$ID
chrom_names <- unique(marker_info_cM_parent[,"CHROM"])

# par1_names_low_geno_error <- par1_names[which(index_low_geno_error)]
# par2_names_low_geno_error <- par2_names[which(index_low_geno_error)]
# off_names_low_geno_error <- off_names[which(index_low_geno_error)]

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



geno <- matrix(c(0,0, 0,1, 1,0, 1,1), ncol=4)
off_genotypes <- array(c(0,0, 0,1, 0,2, 1,0, 1,1, 1,2, 2,0, 2,1, 2,2), dim=c(2, 1, 9))
one_par_geno <- array(NA, dim=c(2, 2, 16))
for (i in 1:4){
  for (j in 1:4){
    one_par_geno[, , (i-1)*4 + j] = array(c(geno[, i], geno[, j]), dim=c(2, 2))
  }
}

score <- array(NA, dim=c(4, 1000, 5, length(chrom_names), length(off_names)))
for (i in 1:length(off_names)){
  par1_geno = as.character(phased_marker_matrix_parent[, colnames(phased_marker_matrix_parent)==
                                                         par1_names[i]])
  par2_geno = as.character(phased_marker_matrix_parent[, colnames(phased_marker_matrix_parent)==
                                                         par2_names[i]])
  off_geno = as.character(marker_matrix[, colnames(marker_matrix)==off_names[i]])
  
  par1A = gsub("\\|[01]", "", par1_geno)
  par1B = gsub("[01]\\|", "", par1_geno)
  par2A = gsub("\\|[01]", "", par2_geno)
  par2B = gsub("[01]\\|", "", par2_geno)
  par1A = as.numeric(par1A)
  par1B = as.numeric(par1B)
  par2A = as.numeric(par2A)
  par2B = as.numeric(par2B)
  
  # off_geno[off_geno=="0/0"] = 0
  # off_geno[off_geno=="0/1"] = 1
  # off_geno[off_geno=="1/1"] = 2
  off_geno = as.numeric(off_geno)
  
  for (j in 1:length(chrom_names)){
    chrom_pos = which(marker_info_cM_parent$CHROM == chrom_names[j])
    off_geno_chrom = off_geno[chrom_pos]
    par1A_chrom = par1A[chrom_pos]
    par1B_chrom = par1B[chrom_pos]
    par2A_chrom = par2A[chrom_pos]
    par2B_chrom = par2B[chrom_pos]
    
    for (h in 1:5){
      dist_markers = cM_dist_markers[[j]][[h]]
      dist_markers = matrix(dist_markers, ncol=2)
      
      if (nrow(dist_markers)!=0){
        
        for (k in 1:nrow(dist_markers)){
          marker_index = dist_markers[k, ]
          chrom_matx = matrix(c(par1A_chrom[marker_index], 
                                par1B_chrom[marker_index], 
                                par2A_chrom[marker_index], 
                                par2B_chrom[marker_index], 
                                off_geno_chrom[marker_index]), ncol=5)
          
          keep_marker = apply(chrom_matx, 1, FUN=get_index_cor_off)
          
          if (!any(is.na(keep_marker))){
            if (all(keep_marker)){
              
              which_par1 = which(apply(one_par_geno, c(3), function(x){
                all(x==chrom_matx[, 1:2])}))
              which_par2 = which(apply(one_par_geno, c(3), function(x){
                all(x==chrom_matx[, 3:4])}))
              which_off = which(apply(off_genotypes, c(3), function(x){
                all(x==chrom_matx[, 5])}))
              
              if (all(c(which_par1, which_par2, which_off))){
                score[1:2, k, h, j, i] = marker_index
                score[3:4, k, h, j, i] = correct_phases[[which_par1]][[which_par2]][[which_off]]
              }
              
            }
          }
          
        } # end of k loop
        
      }
      
    } # end of h loop
  }
}

saveRDS(score, "get_phase_error/score.rds")




# cM_dist_markers <- readRDS("codes/get_markers_dist2_cM_dist_markers.rds")
# for (j in 1:length(cM_dist_markers)){
#   for (h in 1:length(cM_dist_markers[[j]])){
#     rownames(cM_dist_markers[[j]][[h]]) = NULL
#     cM_dist_markers[[j]][[h]] = cM_dist_markers[[j]][[h]][-1, ]
#   }
# }
# for (j in 1:length(cM_dist_markers)){
#   marker_names2 = rownames(marker_info_cM_parent2[marker_info_cM_parent2$CHROM == chrom_names[j], ])
#   marker_names = rownames(marker_info_cM_parent[marker_info_cM_parent$CHROM == chrom_names[j], ])
# 
#   for (h in 1:length(cM_dist_markers[[j]])){
#     cM_dist_markers[[j]][[h]] = apply(cM_dist_markers[[j]][[h]], 2, FUN=function(x){
#       return(match(marker_names2[x], marker_names))
#     })
#   }
# }
# cM_dist_markers2 <- vector(mode="list", length=length(cM_dist_markers))
# for (j in 1:length(cM_dist_markers2)){
#   cM_dist_markers2[[j]] = vector(mode="list", length=length(cM_dist_markers[[j]]))
#   for (h in 1:length(cM_dist_markers2[[j]])){
#     cM_dist_markers2[[j]][[h]] = vector(mode="list", length=2)
#     cM_dist_markers2[[j]][[h]][[1]] = matrix(NA,
#                                              nrow=nrow(cM_dist_markers[[j]][[h]]),
#                                              ncol=ncol(cM_dist_markers[[j]][[h]]))
#     cM_dist_markers2[[j]][[h]][[2]] = matrix(NA,
#                                              nrow=nrow(cM_dist_markers[[j]][[h]]),
#                                              ncol=ncol(cM_dist_markers[[j]][[h]]))
#   }
# }
# for (j in 1:length(cM_dist_markers2)){
#   for (h in 1:length(cM_dist_markers2[[j]])){
#     if_consecutive = apply(cM_dist_markers[[j]][[h]], 1, function(x){
#       ifelse(x[2] - x[1] == 1, T, F)
#     })
#     cM_dist_markers2[[j]][[h]][[1]] = cM_dist_markers[[j]][[h]][if_consecutive, ]
#     cM_dist_markers2[[j]][[h]][[2]] = cM_dist_markers[[j]][[h]][!if_consecutive, ]
#   }
# }
# for (j in 1:length(cM_dist_markers2)){
#   for (h in 1:length(cM_dist_markers2[[j]])){
#     for (g in 1:2){
#       cM_dist_markers2[[j]][[h]][[g]] =
#         cM_dist_markers2[[j]][[h]][[g]][complete.cases(cM_dist_markers2[[j]][[h]][[g]]), ]
#     }
#   }
# }
# for (j in 1:length(cM_dist_markers2)){
#   for (h in 1:length(cM_dist_markers2[[j]])){
#     for (g in 1:2){
#       print(c(j, h, g))
#       print(dim(cM_dist_markers2[[j]][[h]][[g]]))
#     }
#   }
# }
# saveRDS(cM_dist_markers2, "codes/get_phase_error6_cM_dist_markers2.rds")


