setwd("~/family_variance/")
library(vcfR)



get_index_homoz_par_sort <- function(vec){
  # vec is a vector with length 3
  # vec = c(p1, p2, off)
  # it the the genotype (0, 1, or 2) of two parents and their offspring
  if (!any(is.na(vec))){
    if (vec[1]==0 & vec[2]==0){
      return(c(T, F, F))
    } else if ( setequal(vec[1:2], c(0, 2)) ){
      return(c(F, T, F))
    } else if (vec[1]==2 & vec[2]==2){
      return(c(F, F, T))
    } else{
      return(c(F, F, F))
    }
  } else{
    return(c(NA, NA, NA))
  }
}
get_index_homoz_par_cor_off <- function(vec){
  # vec is a vector with length 3
  # vec = c(p1, p2, off)
  # it the the genotype (0, 1, or 2) of two parents and their offspring
  if (!any(is.na(vec))){
    if (vec[3]==1){
      return(T)
    } else{
      return(F)
    }
  } else{
    return(NA)
  }
}

get_index_non_two_het_par <- function(vec){
  # vec is a vector with length 3
  # vec = c(p1, p2, off)
  # it the the genotype (0, 1, or 2) of two parents and their offspring
  if (!any(is.na(vec))){
    if (! (vec[1]==vec[2] & vec[1]==1) ){
      return(T)
    } else {
      return(F)
    }
  } else{
    return(NA)
  }
}
get_index_non_two_het_par_cor_off <- function(vec){
  # vec is a vector with length 3
  # vec = c(p1, p2, off)
  # it the the genotype (0, 1, or 2) of two parents and their offspring
  if (!any(is.na(vec))){
    if (vec[1]==vec[2] & vec[1]==0){
      ifelse(vec[3]==0, T, F)
    } else if (vec[1]==vec[2] & vec[1]==2){
      ifelse(vec[3]==2, T, F)
    } else if (setequal(vec[1:2], c(0, 2))){
      ifelse(vec[3]==1, T, F)
    } else if (setequal(vec[1:2], c(0, 1))){
      ifelse(vec[3]!=2, T, F)
    } else if (setequal(vec[1:2], c(1, 2))){
      ifelse(vec[3]!=0, T, F)
    } else{
      return(NA)
    }
  } else{
    return(NA)
  }
}



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

index_low_geno_error <- readRDS("get_errors/index_off_low_error.rds")
par1_names_low_geno_error <- par1_names[which(index_low_geno_error)]
par2_names_low_geno_error <- par2_names[which(index_low_geno_error)]
off_names_low_geno_error <- off_names[which(index_low_geno_error)]



# find the imputed parental genotypes, set up missing one vs two parental genotypes scenarios 
trio_geno <- array(NA, dim=c(nrow(marker_matrix2), 5, length(off_names_low_geno_error)))
for (i in 1:length(off_names_low_geno_error)){
  off_geno = marker_matrix2[, colnames(marker_matrix2)==off_names_low_geno_error[i]]
  par1_geno = marker_matrix2[, colnames(marker_matrix2)==par1_names_low_geno_error[i]]
  par2_geno = marker_matrix2[, colnames(marker_matrix2)==par2_names_low_geno_error[i]]
  par1_impu = par1_impu2 = 
    phased_marker_matrix_parent2[, colnames(phased_marker_matrix_parent2)==par1_names_low_geno_error[i]]
  par2_impu = par2_impu2 =
    phased_marker_matrix_parent2[, colnames(phased_marker_matrix_parent2)==par2_names_low_geno_error[i]]
  
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

saveRDS(trio_geno, "get_imputated_errors/trio_geno.rds")
# trio_geno <- readRDS("get_imputated_errors/trio_geno.rds")

# look at imputed parental genotypes with 0, 2, when only one missing genotyped parent
index_homoz_par_missing_1par_sort <- 
  apply(trio_geno[, c(1:2, 5), ], c(1, 3), FUN=get_index_homoz_par_sort)
index_homoz_par_missing_1par_02 <- index_homoz_par_missing_1par_sort[2, , ]
dim(index_homoz_par_missing_1par_02) <- 
  c(dim(index_homoz_par_missing_1par_02)[1], 1, dim(index_homoz_par_missing_1par_02)[2])
index_homoz_par_missing_1par_02 <- index_homoz_par_missing_1par_02[, rep(1,3), ]

trio_geno_02 <- array(NA, dim=c(nrow(marker_matrix2), 3, length(off_names_low_geno_error)))
trio_geno_02[which(index_homoz_par_missing_1par_02)] <- 
  trio_geno[, c(1:2, 5), ][which(index_homoz_par_missing_1par_02)]

error_rate_02 <- apply(trio_geno_02, c(1, 3), FUN=get_index_homoz_par_cor_off)
error_rate_02_means <- apply(error_rate_02, 2, mean, na.rm=T)
length(error_rate_02_means)
# [1] 362
sum(is.na(error_rate_02_means))
# [1] 0
mean(error_rate_02_means)
# [1] 0.9061018
1 - 0.9061018
# [1] 0.0938982
mean(error_rate_02, na.rm=T)
# [1] 0.9016224
1 - 0.9016224
# [1] 0.0983776
# imputation error for parental genotypes 0, 2 with one missing parental genotype
pdf("get_imputated_errors/plots/p1.pdf")
hist(error_rate_02_means, 
     xlab="correct imputation rate", 
     main="histogram for average correct imputations across 294 trios \n missing one parental genotype")
plot(density(na.omit(error_rate_02_means)), 
     xlab="correct imputation rate", 
     main="density for average correct imputations across 294 trios")
dev.off()

# index_off_low_error_missing_1par <- apply(error_rate_02, 2, function(x){
#   ifelse(mean(x, na.rm=T) > 0.95, T, F)
# })
# length(which(index_off_low_error_missing_1par))
# # [1] 159
# 
# saveRDS(index_off_low_error_missing_1par, 
#         "codes/get_imputated_errors_index_off_low_error_missing_1par.rds")

# do the same thing for missing both parental genotypes
# not a lot of markers where both parents were not genotyped
index_homoz_par_missing_2par_sort <- 
  apply(trio_geno[, 3:5, ], c(1, 3), FUN=get_index_homoz_par_sort)
index_homoz_par_missing_2par_02 <- index_homoz_par_missing_2par_sort[2, , ]
dim(index_homoz_par_missing_2par_02) <- 
  c(dim(index_homoz_par_missing_2par_02)[1], 1, dim(index_homoz_par_missing_2par_02)[2])
index_homoz_par_missing_2par_02 <- index_homoz_par_missing_2par_02[, rep(1,3), ]

trio_geno_missing_2par_02 <- array(NA, dim=c(nrow(marker_matrix2), 3, length(off_names_low_geno_error)))
trio_geno_missing_2par_02[which(index_homoz_par_missing_2par_02)] <- 
  trio_geno[, 3:5, ][which(index_homoz_par_missing_2par_02)]

error_rate_missing_2par_02 <- apply(trio_geno_missing_2par_02, c(1, 3), FUN=get_index_homoz_par_cor_off)
error_rate_missing_2par_02_means <- apply(error_rate_missing_2par_02, 2, mean, na.rm=T)
length(error_rate_missing_2par_02_means)
# [1] 362
sum(is.na(error_rate_missing_2par_02_means))
# [1] 155
sum(!is.na(error_rate_missing_2par_02_means))
# [1] 207
mean(error_rate_missing_2par_02_means, na.rm=T)
# [1] 0.8756039
1 - 0.8756039
# [1] 0.1243961
mean(error_rate_missing_2par_02, na.rm=T)
# [1] 0.8983957
1 - 0.8983957
# [1] 0.1016043
# imputation error for parental genotypes 0, 2 with both missing parental genotypes
pdf("get_imputated_errors/plots/p2.pdf")
hist(error_rate_missing_2par_02_means, 
     xlab="correct imputation rate", 
     main="histogram for average correct imputations across 126 trios \n missing two parental genotype")
plot(density(na.omit(error_rate_missing_2par_02_means)), 
     xlab="correct imputation rate", 
     main="density for average correct imputations across 126 trios")
dev.off()



# error rate when at least one parent isn't heterozygous
index_non_two_het_par_missing_1par <- 
  apply(trio_geno[, c(1:2,5), ], c(1, 3), FUN=get_index_non_two_het_par)
dim(index_non_two_het_par_missing_1par) <- 
  c(dim(index_non_two_het_par_missing_1par)[1], 1, dim(index_non_two_het_par_missing_1par)[2])
index_non_two_het_par_missing_1par <- index_non_two_het_par_missing_1par[, rep(1,3), ]

trio_geno_missing_1par_non_two_het_par <- 
  array(NA, dim=c(nrow(marker_matrix2), 3, length(off_names_low_geno_error)))
trio_geno_missing_1par_non_two_het_par[which(index_non_two_het_par_missing_1par)] <- 
  trio_geno[, c(1:2,5), ][which(index_non_two_het_par_missing_1par)]

error_rate_missing_1par_non_two_het_par <- 
  apply(trio_geno_missing_1par_non_two_het_par, c(1, 3), FUN=get_index_non_two_het_par_cor_off)
error_rate_missing_1par_non_two_het_par_means <- 
  apply(error_rate_missing_1par_non_two_het_par, 2, mean, na.rm=T)
mean(error_rate_missing_1par_non_two_het_par_means, na.rm=T)
# [1] 0.9595174
1 - 0.9595174
# [1] 0.0404826
pdf("get_imputated_errors/plots/p3.pdf")
hist((1 - error_rate_missing_1par_non_two_het_par_means), 
     main="histogram for average imputations error across 295 trios \n missing one parental genotype")
plot(density(na.omit(1 - error_rate_missing_1par_non_two_het_par_means)), 
     main="density for average imputations error across 295 trios")
dev.off()



index_non_two_het_par_missing_2par <- 
  apply(trio_geno[, 3:5, ], c(1, 3), FUN=get_index_non_two_het_par)
dim(index_non_two_het_par_missing_2par) <- 
  c(dim(index_non_two_het_par_missing_2par)[1], 1, dim(index_non_two_het_par_missing_2par)[2])
index_non_two_het_par_missing_2par <- index_non_two_het_par_missing_2par[, rep(1,3), ]

trio_geno_missing_2par_non_two_het_par <- 
  array(NA, dim=c(nrow(marker_matrix2), 3, length(off_names_low_geno_error)))
trio_geno_missing_2par_non_two_het_par[which(index_non_two_het_par_missing_2par)] <- 
  trio_geno[, 3:5, ][which(index_non_two_het_par_missing_2par)]

error_rate_missing_2par_non_two_het_par <- 
  apply(trio_geno_missing_2par_non_two_het_par, c(1, 3), FUN=get_index_non_two_het_par_cor_off)
error_rate_missing_2par_non_two_het_par_means <- 
  apply(error_rate_missing_2par_non_two_het_par, 2, mean, na.rm=T)
mean(error_rate_missing_2par_non_two_het_par_means, na.rm=T)
# [1] 0.8615724
1 - 0.8615724
# [1] 0.1384276
pdf("get_imputated_errors/plots/p4.pdf")
hist((1 - error_rate_missing_2par_non_two_het_par_means), 
     main="histogram for average imputations error across 295 trios \n missing two parental genotype")
plot(density(na.omit(1 - error_rate_missing_2par_non_two_het_par_means)), 
     main="density for average imputations error across 295 trios")
dev.off()



# head(error_rate_missing_2par_02_means)
# # [1] NaN   1   1 NaN NaN NaN
# sum(is.na(error_rate_missing_2par_02[, 1]))
# # [1] 37254
# sum(is.na(trio_geno_missing_2par_02[, , 1])) / 3
# # [1] 37254
# sum(index_homoz_par_missing_2par_sort[2, , 1], na.rm=T)
# # [1] 0
# 
# sum(!is.na(trio_geno[, 3, 1]))
# # [1] 15
# trio_geno[!is.na(trio_geno[, 3, 1]), 3:5, 1]
# #       [,1] [,2] [,3]
# #  [1,]    0    0    0
# #  [2,]    0    0   NA
# #  [3,]    0    1    1
# #  [4,]    2    1    1
# #  [5,]    0    0    0
# #  [6,]    1    2    2
# #  [7,]    1    1    1
# #  [8,]    0    1    1
# #  [9,]    0    0   NA
# # [10,]    1    1    1
# # [11,]    1    1    2
# # [12,]    1    1    1
# # [13,]    1    1    1
# # [14,]    1    1    2
# # [15,]    1    1    2




# index_off_low_error_missing_2par <- apply(error_rate_missing_2par_02, 2, function(x){
#   ifelse(mean(x, na.rm=T) > 0.95, T, F)
# })
# length(which(index_off_low_error_missing_2par))
# # [1] 128
# 
# saveRDS(index_off_low_error_missing_2par, 
#         "codes/get_imputated_errors_index_off_low_error_missing_2par.rds")
