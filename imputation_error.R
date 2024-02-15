setwd("family_variance/")
library(vcfR)



gen_off_geno <- function(vec){
  # vec = c(p1, p2)
  # p1, p2 take values of 0, 1, 2
  if (setequal(vec, c(0, 0))){
    poss_off_geno = 0
  } else if (setequal(vec, c(0, 1))){
    poss_off_geno = c(0, 1)
  } else if (setequal(vec, c(0, 2))){
    poss_off_geno = 1
  } else if (setequal(vec, c(1, 1))){
    poss_off_geno = c(0, 1, 2)
  } else if (setequal(vec, c(1, 2))){
    poss_off_geno = c(1, 2)
  } else if (setequal(vec, c(2, 2))){
    poss_off_geno = 2
  } else {
    poss_off_geno = NA
    print("unvalid genotype input")
  }
  return(poss_off_geno)
}
if_genotyped <- function(vec){
  # vec = c(p1, p2, off)
  if (!any(is.na(vec))){
    return(T)
  } else {
    return(F)
  }
}



{
marker_info_valid <- read.csv("generate_vcffiles/marker_info_valid.csv", row.names=1, check.names=F)
marker_matrix_valid <- read.csv("generate_vcffiles/marker_matrix_valid.csv", row.names=1, check.names=F)
pedigree_valid <- read.csv("generate_vcffiles/pedigree_valid.csv")
phased_marker_parent <- read.vcfR("phased_data/phased_parent.vcf")
phased_marker_info_parent <- as.data.frame(phased_marker_parent@fix)
phased_marker_matrix_parent <- as.data.frame(phased_marker_parent@gt)

phased_marker_info_parent$FORMAT <- phased_marker_matrix_parent$FORMAT
phased_marker_matrix_parent <- phased_marker_matrix_parent[, -1]

dim(marker_info_valid)
# [1] 37285     9
dim(marker_matrix_valid)
# [1] 37285   436
dim(pedigree_valid)
# [1] 436   3
dim(phased_marker_info_parent)
# [1] 37285     9
dim(phased_marker_matrix_parent)
# [1] 37285   563



marker_matrix_valid2 <- matrix(NA, nrow=nrow(marker_matrix_valid), ncol=ncol(marker_matrix_valid))
marker_matrix_valid2[marker_matrix_valid=="0/0"] <- 0
marker_matrix_valid2[marker_matrix_valid=="0/1"] <- 1
marker_matrix_valid2[marker_matrix_valid=="1/1"] <- 2
colnames(marker_matrix_valid2) <- colnames(marker_matrix_valid)

phased_marker_matrix_parent2 <- matrix(NA, 
                                       nrow=nrow(phased_marker_matrix_parent), 
                                       ncol=ncol(phased_marker_matrix_parent))
phased_marker_matrix_parent2[phased_marker_matrix_parent=="0|0"] <- 0
phased_marker_matrix_parent2[phased_marker_matrix_parent=="0|1"] <- 1
phased_marker_matrix_parent2[phased_marker_matrix_parent=="1|0"] <- 1
phased_marker_matrix_parent2[phased_marker_matrix_parent=="1|1"] <- 2
colnames(phased_marker_matrix_parent2) <- colnames(phased_marker_matrix_parent)



off_names <- pedigree_valid$Accession_ID
par1_names <- pedigree_valid$Integrated_P1
par2_names <- pedigree_valid$Integrated_P2
chrom_names <- unique(phased_marker_info_parent[,"CHROM"])
}



imputation_error <- array(NA, dim=c(4, 2000, length(chrom_names), length(off_names)))
for (i in 1){
  off_geno = marker_matrix_valid2[, colnames(marker_matrix_valid2)==off_names[i]]
  par1_geno = phased_marker_matrix_parent2[, colnames(phased_marker_matrix_parent2)==par1_names[i]]
  par2_geno = phased_marker_matrix_parent2[, colnames(phased_marker_matrix_parent2)==par2_names[i]]
  
  for (j in 1){
    chrom_pos = which(marker_info_valid$CHROM == chrom_names[j])
    off_chrom = off_geno[chrom_pos]
    par1_chrom = par1_geno[chrom_pos]
    par2_chrom = par1_geno[chrom_pos]
    
    trio_chrom = matrix(c(par1_chrom, par2_chrom, off_chrom), nrow=3, byrow=T)
    genotyped_pos = which(apply(trio_chrom, 2, FUN=if_genotyped))
    
    if (length(genotyped_pos)!=0){
      error_pos = c()
      for (k in 1:length(genotyped_pos)){
        if (!(off_chrom[genotyped_pos[k]] %in% gen_off_geno(trio_chrom[1:2, genotyped_pos[k]]))){
          error_pos = c(error_pos, genotyped_pos[k])
        }
      }
      if (length(error_pos)!=0){
        imputation_error[1, 1:length(error_pos), j, i] = error_pos
        imputation_error[2, 1:length(error_pos), j, i] = off_chrom[error_pos]
        imputation_error[3:4, 1:length(error_pos), j, i] = trio_chrom[1:2, error_pos]
      }
    }
  }
}



saveRDS(imputation_error, file="imputation_error/imputation_error.rds")

