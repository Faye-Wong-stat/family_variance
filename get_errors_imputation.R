setwd("~/family_variance/")
library(vcfR)



marker_info <- readRDS("generate_vcffiles/marker_info.rds")
marker_matrix <- readRDS("generate_vcffiles/marker_matrix.rds")
pedigree_valid <- readRDS("generate_vcffiles/pedigree_valid.rds")
phased_marker_parent <- read.vcfR("phased_data/phased_parent.vcf")
phased_marker_info_parent <- as.data.frame(phased_marker_parent@fix)
phased_marker_matrix_parent <- as.data.frame(phased_marker_parent@gt)
phased_marker_info_parent$FORMAT <- phased_marker_matrix_parent$FORMAT
phased_marker_matrix_parent <- phased_marker_matrix_parent[, -1]

marker_matrix2 <- matrix(NA, nrow=nrow(marker_matrix), ncol=ncol(marker_matrix))
marker_matrix2[marker_matrix=="0/0"] <- 0
marker_matrix2[marker_matrix=="0/1"] <- 1
marker_matrix2[marker_matrix=="1/1"] <- 2
rownames(marker_matrix2) <- rownames(marker_matrix)
colnames(marker_matrix2) <- colnames(marker_matrix)

phased_marker_matrix_parent2 <- 
  matrix(NA, nrow=nrow(phased_marker_matrix_parent), ncol=ncol(phased_marker_matrix_parent))
phased_marker_matrix_parent2[phased_marker_matrix_parent=="0|0"] <- 0
phased_marker_matrix_parent2[phased_marker_matrix_parent=="0|1"] <- 1
phased_marker_matrix_parent2[phased_marker_matrix_parent=="1|0"] <- 1
phased_marker_matrix_parent2[phased_marker_matrix_parent=="1|1"] <- 2
rownames(phased_marker_matrix_parent2) <- phased_marker_info_parent$ID
colnames(phased_marker_matrix_parent2) <- colnames(phased_marker_matrix_parent)

par1_names <- pedigree_valid$P1
par2_names <- pedigree_valid$P2
off_names <- pedigree_valid$ID



trio_geno <- array(NA, dim=c(nrow(marker_matrix2), 5, length(off_names)))
for (i in 1:length(off_names)){
  off_geno = marker_matrix2[, colnames(marker_matrix2)==off_names[i]]
  par1_geno = marker_matrix2[, colnames(marker_matrix2)==par1_names[i]]
  par2_geno = marker_matrix2[, colnames(marker_matrix2)==par2_names[i]]
  par1_impu = par1_impu2 = 
    phased_marker_matrix_parent2[, colnames(phased_marker_matrix_parent2)==par1_names[i]]
  par2_impu = par2_impu2 =
    phased_marker_matrix_parent2[, colnames(phased_marker_matrix_parent2)==par2_names[i]]
  
  index_par1_missing = is.na(par1_geno)
  index_par2_missing = is.na(par2_geno)
  index_missing_1par = index_par1_missing | index_par2_missing
  index_missing_2par = index_par1_missing & index_par2_missing
  index_missing_1par = ifelse(index_missing_2par, F, index_missing_1par)
  
  par1_impu[!index_missing_1par] = NA
  par2_impu[!index_missing_1par] = NA
  par1_impu2[!index_missing_2par] = NA
  par2_impu2[!index_missing_2par] = NA
  # print(sum(!is.na(par1_impu)))
  # print(sum(!is.na(par1_impu2)))
  trio_geno[, , i] = cbind(par1_impu, par2_impu, par1_impu2, par2_impu2, off_geno)
}



error_rate_missing1 <- apply(trio_geno[, c(1, 2, 5), ], c(1, 3), FUN=function(x){
  if (!any(is.na(x))) {
    
    if (x[1]==x[2] & x[1]!=1){
      if (x[2]==x[3]){
        return(T)
      } else {
        return(F)
      }
    } else if (setequal(x[1:2], c(0, 2))){
      if (x[3]==1){
        return(T)
      } else {
        return(F)
      }
    } else {
      return(NA) 
    }
    
  } else {
    return(NA)
  }
})

error_rate_missing2 <- apply(trio_geno[, c(3, 4, 5), ], c(1, 3), FUN=function(x){
  if (!any(is.na(x))) {
    
    if (x[1]==x[2] & x[1]!=1){
      if (x[2]==x[3]){
        return(T)
      } else {
        return(F)
      }
    } else if (setequal(x[1:2], c(0, 2))){
      if (x[3]==1){
        return(T)
      } else {
        return(F)
      }
    } else {
      return(NA) 
    }
    
  } else {
    return(NA)
  }
})



saveRDS(trio_geno, "get_errors_imputation/trio_geno.rds")
saveRDS(error_rate_missing1, "get_errors_imputation/error_rate_missing1.rds")
saveRDS(error_rate_missing2, "get_errors_imputation/error_rate_missing2.rds")



mean(error_rate_missing1, na.rm=T)
# [1] 0.9388043
mean(error_rate_missing2, na.rm=T)
# [1] 0.8741175

error_rate_missing1_mean <- apply(error_rate_missing1, 2, mean, na.rm=T)
error_rate_missing2_mean <- apply(error_rate_missing2, 2, mean, na.rm=T)

pdf("get_errors_imputation/plots/p1.pdf")
hist(error_rate_missing1_mean, xlab="percentage of correct offspring imputation per trio", 
     main="when one parent is missing")
dev.off()
pdf("get_errors_imputation/plots/p2.pdf")
hist(error_rate_missing2_mean, xlab="percentage of correct offspring imputation per trio", 
     main="when two parents are missing")
dev.off()














