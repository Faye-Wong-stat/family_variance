setwd("family_variance/")



marker_info_valid <- read.csv("generate_vcffiles/marker_info_valid.csv", row.names=1, check.names=F)
marker_matrix_valid <- read.csv("generate_vcffiles/marker_matrix_valid.csv", row.names=1, check.names=F)
pedigree_valid <- read.csv("generate_vcffiles/pedigree_valid.csv")



off_names <- pedigree_valid$Accession_ID
par1_names <- pedigree_valid$Integrated_P1
par2_names <- pedigree_valid$Integrated_P2
chrom_names <- unique(marker_info_valid[,"CHROM"])



genotype_error <- readRDS("genotype_error/genotype_error.rds")
imputation_error <- readRDS("imputation_error/imputation_error.rds")
# low_poss <- readRDS("codes/view_prob_two_loci_low_poss.rds")



# error <- array("", dim=c(12, 2000, 28, 436))
error <- array("", dim=c(8, 2000, 28, 436))
for (i in 1:436){
  off_geno = as.character(marker_matrix_valid[, colnames(marker_matrix_valid)==off_names[i]])
  
  for (j in 1:28){
    chrom_pos = which(marker_info_valid$CHROM == chrom_names[j])
    off_chrom = off_geno[chrom_pos]
    
    index_geno = genotype_error[1, , j, i]
    index_impu = imputation_error[1, , j, i]
    # index_phas = as.numeric(low_poss[2, , j, i])
    index_geno = index_geno[!is.na(index_geno)]
    index_impu = index_impu[!is.na(index_impu)]
    index_impu = index_impu[!index_impu %in% index_geno]
    # index_phas = index_phas[!is.na(index_phas)]

    # index_un = union(union(index_geno, index_impu), index_phas)
    index_un = union(index_geno, index_impu)
    index_un = index_un[order(index_un)]
    
    if (length(index_un)!=0){
      error[1, 1:length(index_un), j, i] = as.character(index_un)
      
      error[2, 1:length(index_un), j, i] = off_chrom[index_un]
      
      error[3, 1:length(index_un), j, i] = as.character(index_geno[match(index_un, index_geno)])
      for (k in 1:2){
        error[k+3, 1:length(index_un), j, i] = 
          as.character(genotype_error[k+2, , j, i][match(index_un, index_geno)])
      }
      
      error[6, 1:length(index_un), j, i] = as.character(index_impu[match(index_un, index_impu)])
      for (k in 1:2){
        error[k+6, 1:length(index_un), j, i] = 
          as.character(imputation_error[k+2, , j, i][match(index_un, index_impu)])
      }
      
      # for (k in 1:2){
      #   error[k+8, 1:length(index_un), j, i] = 
      #     low_poss[k, , j, i][match(index_un, index_phas)]
      # }
      # for (k in 1:2){
      #   error[k+10, 1:length(index_un), j, i] = 
      #     low_poss[k+3, , j, i][match(index_un, index_phas)]
      # }
    }
  }
}

saveRDS(error, file="codes/view_errors.rds")






