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



marker_matrix <- read.csv("generate_vcffiles/marker_matrix.csv", row.names=1, check.names=F)
marker_info_cM_parent <- read.csv("get_markers_cM/marker_info_cM_parent.csv", row.names=1)
marker_info_cM_valid <- read.csv("get_markers_cM/marker_info_cM_valid.csv", row.names=1)
# marker_info_cM_parent <- read.csv("cleaned_data/marker_info_cM_parent.csv", row.names=1)
# marker_info_cM_valid <- read.csv("cleaned_data/marker_info_cM_valid.csv", row.names=1)

pedigree_valid <- read.csv("generate_vcffiles/pedigree_valid.csv")
index_low_geno_error <- readRDS("get_errors/index_off_low_error.rds")
# marker_matrix <- read.csv("cleaned_data/marker_matrix.csv", row.names=1, check.names=F)

phased_marker_parent <- read.vcfR("phased_data/phased_parent.vcf")
phased_marker_matrix_parent <- as.data.frame(phased_marker_parent@gt)
phased_marker_matrix_parent <- phased_marker_matrix_parent[, -1]

par1_names <- pedigree_valid$Integrated_P1
par2_names <- pedigree_valid$Integrated_P2
off_names <- pedigree_valid$Accession_ID
chrom_names <- unique(marker_info_cM_parent[,"CHROM"])

par1_names_low_geno_error <- par1_names[which(index_low_geno_error)]
par2_names_low_geno_error <- par2_names[which(index_low_geno_error)]
off_names_low_geno_error <- off_names[which(index_low_geno_error)]

marker_info_cM_parent2 <- marker_info_cM_parent[!is.na(marker_info_cM_parent$cM), ]


cM_dist_markers <- readRDS("get_markers_dist/cM_dist_markers.rds")
# cM_dist_markers <- readRDS("codes/get_markers_dist2_cM_dist_markers.rds")
for (j in 1:length(cM_dist_markers)){
  for (h in 1:length(cM_dist_markers[[j]])){
    rownames(cM_dist_markers[[j]][[h]]) = NULL
    cM_dist_markers[[j]][[h]] = cM_dist_markers[[j]][[h]][-1, ]
  }
}
for (j in 1:length(cM_dist_markers)){
  marker_names2 = rownames(marker_info_cM_parent2[marker_info_cM_parent2$CHROM == chrom_names[j], ])
  marker_names = rownames(marker_info_cM_parent[marker_info_cM_parent$CHROM == chrom_names[j], ])
  
  for (h in 1:length(cM_dist_markers[[j]])){
    cM_dist_markers[[j]][[h]] = apply(cM_dist_markers[[j]][[h]], 2, FUN=function(x){
      return(match(marker_names2[x], marker_names))
    })
  }
}
cM_dist_markers2 <- vector(mode="list", length=length(cM_dist_markers))
for (j in 1:length(cM_dist_markers2)){
  cM_dist_markers2[[j]] = vector(mode="list", length=length(cM_dist_markers[[j]]))
  for (h in 1:length(cM_dist_markers2[[j]])){
    cM_dist_markers2[[j]][[h]] = vector(mode="list", length=2)
    cM_dist_markers2[[j]][[h]][[1]] = matrix(NA, 
                                             nrow=nrow(cM_dist_markers[[j]][[h]]), 
                                             ncol=ncol(cM_dist_markers[[j]][[h]])) 
    cM_dist_markers2[[j]][[h]][[2]] = matrix(NA, 
                                             nrow=nrow(cM_dist_markers[[j]][[h]]), 
                                             ncol=ncol(cM_dist_markers[[j]][[h]])) 
  }
}
for (j in 1:length(cM_dist_markers2)){
  for (h in 1:length(cM_dist_markers2[[j]])){
    if_consecutive = apply(cM_dist_markers[[j]][[h]], 1, function(x){
      ifelse(x[2] - x[1] == 1, T, F)
    })
    cM_dist_markers2[[j]][[h]][[1]] = cM_dist_markers[[j]][[h]][if_consecutive, ]
    cM_dist_markers2[[j]][[h]][[2]] = cM_dist_markers[[j]][[h]][!if_consecutive, ]
  }
}
for (j in 1:length(cM_dist_markers2)){
  for (h in 1:length(cM_dist_markers2[[j]])){
    for (g in 1:2){
      cM_dist_markers2[[j]][[h]][[g]] = 
        cM_dist_markers2[[j]][[h]][[g]][complete.cases(cM_dist_markers2[[j]][[h]][[g]]), ]
    }
  }
}
for (j in 1:length(cM_dist_markers2)){
  for (h in 1:length(cM_dist_markers2[[j]])){
    for (g in 1:2){
      print(c(j, h, g))
      print(dim(cM_dist_markers2[[j]][[h]][[g]]))
    }
  }
}



saveRDS(correct_phases, "codes/get_phase_error6_correct_phases.rds")
saveRDS(cM_dist_markers2, "codes/get_phase_error6_cM_dist_markers2.rds")


