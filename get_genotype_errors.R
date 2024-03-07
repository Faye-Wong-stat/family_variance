setwd("~/family_variance/")
# library(vcfR)



get_missing_genotype <- function(vec){
  # vec is the genotype ("0/0", etc) of two parents
  # vec = c(p1, p2)
  if (sum(is.na(vec))==1){
    return(1)
  } else if (sum(is.na(vec))==2){
    return(2)
  } else if (sum(is.na(vec))==0 & !any(vec=="")){
    return(0)
  } else {
    return(NA)
  }
}
# get the index of homozygous parents 
# add if(length(vec)==3)
get_index_homo_par <- function(vec){
  # vec is the genotype (0, 1, 2) of the offspring and two parents
  # vec = c(off, p1, p2)
  if ((!any(is.na(vec)))){
    if (setequal(vec[2:3], c(0, 2)) | (vec[2]==0 & vec[3]==0) | (vec[2]==2 & vec[3]==2)){
      return(T)
    } else{
      return(F)
    }
  } else{
    return(F)
  }
}
get_index_homo_par_cor_off <- function(vec){
  # vec is the genotype (0, 1, 2) of the offspring and two parents
  # vec = c(off, p1, p2)
  if (vec[2]==vec[3] & vec[2]==vec[1]){
    return(T)
  } else if (setequal(vec[2:3], c(0, 2)) & vec[1]==1){
    return(T)
  } else{
    return(F)
  }
}
get_index_homo_par_sort <- function(vec){
  # vec is the genotype (0, 1, 2) of the two parents
  # vec = c(p1, p2)
  if ((!any(is.na(vec)))){
    if ((vec[1]==0 & vec[2]==0)){
      return(c(T, F, F))
    } else if (vec[1]==2 & vec[2]==2){
      return(c(F, F, T))
    } else if (setequal(vec[1:2], c(0, 2))){
      return(c(F, T, F))
    } else{
      return(c(F, F, F))
    }
  } else{
    return(c(F, F, F))
  }
}



marker_info <- read.csv("generate_vcffiles/marker_info.csv", row.names=1, check.names=F)
marker_matrix <- read.csv("generate_vcffiles/marker_matrix.csv", row.names=1, check.names=F)
pedigree_valid <- read.csv("generate_vcffiles/pedigree_valid.csv")

# marker_info_valid <- read.csv("cleaned_data/marker_info_valid.csv", row.names=1, check.names=F)
# marker_matrix_valid <- read.csv("cleaned_data/marker_matrix_valid.csv", row.names=1, check.names=F)
# marker_info_parent <- read.csv("cleaned_data/marker_info_parent.csv", row.names=1, check.names=F)
# marker_matrix_parent <- read.csv("cleaned_data/marker_matrix_parent.csv", row.names=1, check.names=F)
# phased_marker_parent <- read.vcfR("phased_data/phased_parent.vcf")
# phased_marker_info_parent <- as.data.frame(phased_marker_parent@fix)
# phased_marker_matrix_parent <- as.data.frame(phased_marker_parent@gt)
# pedigree_valid <- read.csv("cleaned_data/pedigree_end_progenies.csv")

# phased_marker_info_parent$FORMAT <- phased_marker_matrix_parent$FORMAT
# phased_marker_matrix_parent <- phased_marker_matrix_parent[, -1]

marker_matrix2 <- matrix(NA, nrow=nrow(marker_matrix), ncol=ncol(marker_matrix))
marker_matrix2[marker_matrix=="0/0"] <- 0
marker_matrix2[marker_matrix=="0/1"] <- 1
marker_matrix2[marker_matrix=="1/1"] <- 2
colnames(marker_matrix2) <- colnames(marker_matrix)
rownames(marker_matrix2) <- rownames(marker_matrix)

# marker_matrix_parent2 <- matrix(NA, nrow=nrow(marker_matrix_parent), ncol=ncol(marker_matrix_parent))
# marker_matrix_parent2[marker_matrix_parent=="0/0"] <- 0
# marker_matrix_parent2[marker_matrix_parent=="0/1"] <- 1
# marker_matrix_parent2[marker_matrix_parent=="1/1"] <- 2
# colnames(marker_matrix_parent2) <- colnames(marker_matrix_parent)
# 
# phased_marker_matrix_parent2 <- matrix(NA,
#                                        nrow=nrow(phased_marker_matrix_parent),
#                                        ncol=ncol(phased_marker_matrix_parent))
# phased_marker_matrix_parent2[phased_marker_matrix_parent=="0|0"] <- 0
# phased_marker_matrix_parent2[phased_marker_matrix_parent=="0|1"] <- 1
# phased_marker_matrix_parent2[phased_marker_matrix_parent=="1|0"] <- 1
# phased_marker_matrix_parent2[phased_marker_matrix_parent=="1|1"] <- 2
# colnames(phased_marker_matrix_parent2) <- colnames(phased_marker_matrix_parent)

off_names <- pedigree_valid$Accession_ID
par1_names <- pedigree_valid$Integrated_P1
par2_names <- pedigree_valid$Integrated_P2
chrom_names <- unique(marker_info[,"CHROM"])
set.seed(1)
off_names_shuff <- sample(off_names, replace=F)
head(off_names_shuff)
# [1] 16C005P028 16C536P008 16C050P005 16C017P035 16C061P028 16C048P012


{
# length(off_names)
# length(colnames(marker_matrix_valid))
# sum(!off_names %in% colnames(marker_matrix_valid))
# sum(!colnames(marker_matrix_valid) %in% off_names)
# off_names[!off_names %in% colnames(marker_matrix_valid)]
# colnames(marker_matrix_valid)[!colnames(marker_matrix_valid) %in% off_names]
}


# how many parents genotypes are imputated? (missing)
# trio_genotype in format of "0/0", "0/1", "1/1"
trio_genotype <- array("", dim=c(3, 2000, 28, 436))
for (i in 1:436){
  off_geno = as.character(marker_matrix[, colnames(marker_matrix)==off_names[i]])
  par1_geno = as.character(marker_matrix[, colnames(marker_matrix)==par1_names[i]])
  par2_geno = as.character(marker_matrix[, colnames(marker_matrix)==par2_names[i]])
  
  for (j in 1:28){
    chrom_pos = marker_info$CHROM == chrom_names[j]
    off_chrom = off_geno[chrom_pos]
    par1_chrom = par1_geno[chrom_pos]
    par2_chrom = par1_geno[chrom_pos]
    
    trio_chrom = matrix(c(off_chrom, par1_chrom, par2_chrom), nrow=3, byrow=T)
    trio_genotype[, 1:dim(trio_chrom)[2], j, i] = trio_chrom
  }
}

saveRDS(trio_genotype, file="get_genotype_errors/trio_genotype.rds")
# trio_genotype <- readRDS("get_genotype_errors/trio_genotype.rds")



missing_genotype <- apply(trio_genotype[2:3, , , ], c(2, 3, 4), FUN=get_missing_genotype)
sum(missing_genotype==0, na.rm=T)
# [1] 15,914,324
sum(missing_genotype==1, na.rm=T)
# [1] 0
sum(missing_genotype==2, na.rm=T)
# [1] 341,936
15914324 + 341936 == 436 * 37285
# TRUE
15914324 + 341936
# [1] 16,256,260
# total number of markers

missing_0_parents_off <- 0
for (i in 1:436){
  for (j in 1:28){
    off_geno = trio_genotype[1, which(missing_genotype[, j, i]==0), j, i]
    missing_0_parents_off = missing_0_parents_off + sum(is.na(off_geno))
  }
}
# out of the 15,914,324 markers that both parental markers are genotyped, 
# 140,911 of the offspring markers are not genotyped 
missing_2_parents_off <- 0
for (i in 1:436){
  for (j in 1:28){
    off_geno = trio_genotype[1, which(missing_genotype[, j, i]==2), j, i]
    missing_2_parents_off = missing_2_parents_off + sum(is.na(off_geno))
  }
}
# out of the 341,292 markers that neither parental markers is genotyped, 
# 12,959 of the offspring markers are not genotyped 



# trio_genotype2 in format of 0, 1, 2
trio_genotype2 <- array(NA, dim=c(3, 2000, 28, 436))
for (i in 1:436){
  off_geno = marker_matrix2[, colnames(marker_matrix2)==off_names[i]]
  par1_geno = marker_matrix2[, colnames(marker_matrix2)==par1_names[i]]
  par2_geno = marker_matrix2[, colnames(marker_matrix2)==par2_names[i]]
  
  for (j in 1:28){
    chrom_pos = marker_info$CHROM == chrom_names[j]
    off_chrom = off_geno[chrom_pos]
    par1_chrom = par1_geno[chrom_pos]
    par2_chrom = par1_geno[chrom_pos]
    
    trio_chrom = matrix(c(off_chrom, par1_chrom, par2_chrom), nrow=3, byrow=T)
    trio_genotype2[, 1:dim(trio_chrom)[2], j, i] = trio_chrom
  }
}

saveRDS(trio_genotype2, file="get_genotype_errors/trio_genotype2.rds")
# trio_genotype2 <- readRDS("get_genotype_errors/trio_genotype2.rds")

# trio_genotype3 is the trios with homozygous parents 
trio_genotype3 <- array(NA, dim=dim(trio_genotype2))
homo_par_num <- matrix(NA, nrow=28, ncol=436)
for (i in 1:436){
  for (j in 1:28){
    index_homo_par = apply(trio_genotype2[, , j, i], 2, FUN=get_index_homo_par)
    index_homo_par = which(index_homo_par)
    homo_par_num[j, i] = length(index_homo_par)
    trio_genotype3[, 1:length(index_homo_par), j, i] = trio_genotype2[, index_homo_par, j, i]
  }
}
sum(homo_par_num)
# [1] 10,846,295

homo_par_cor_off_num <- matrix(NA, nrow=28, ncol=436)
for (i in 1:436){
  for (j in 1:28){
    index_cor_off = 
      apply(trio_genotype3[, 1:homo_par_num[j, i], j, i], 2, FUN=get_index_homo_par_cor_off)
    index_incor_off = which(!index_cor_off)
    index_cor_off = which(index_cor_off)
    homo_par_cor_off_num[j, i] = length(index_cor_off)
    
  }
}
trio_genotype3[, head(index_incor_off, 20), j, i]
sum(homo_par_cor_off_num)
# [1] 8,708,064
1 - 8708064 / 10846295
# [1] 0.1971393
# genotyping error rate with homozygous parental markers and offspring marker

saveRDS(trio_genotype3, file="get_genotype_errors/trio_genotype3.rds")
saveRDS(homo_par_num, file="get_genotype_errors/homo_par_num.rds")
saveRDS(homo_par_cor_off_num, file="get_genotype_errors/homo_par_cor_off_num.rds")

# trio_genotype3 <- readRDS("get_genotype_errors/trio_genotype3.rds")
# homo_par_num <- readRDS("get_genotype_errors/homo_par_num.rds")
# homo_par_cor_off_num <- readRDS("get_genotype_errors/homo_par_cor_off_num.rds")



# trio_genotype4 in format of 0, 1, 2, but with shuffled offspring
trio_genotype4 <- array(NA, dim=c(3, 2000, 28, 436))
for (i in 1:436){
  off_geno = marker_matrix2[, colnames(marker_matrix2)==off_names_shuff[i]]
  par1_geno = marker_matrix2[, colnames(marker_matrix2)==par1_names[i]]
  par2_geno = marker_matrix2[, colnames(marker_matrix2)==par2_names[i]]
  
  for (j in 1:28){
    chrom_pos = marker_info$CHROM == chrom_names[j]
    off_chrom = off_geno[chrom_pos]
    par1_chrom = par1_geno[chrom_pos]
    par2_chrom = par1_geno[chrom_pos]
    
    trio_chrom = matrix(c(off_chrom, par1_chrom, par2_chrom), nrow=3, byrow=T)
    trio_genotype4[, 1:dim(trio_chrom)[2], j, i] = trio_chrom
  }
}

# trio_genotype5 is the trios with homozygous parents and shuffled offsprings
trio_genotype5 <- array(NA, dim=dim(trio_genotype4))
homo_par_num_shuff <- matrix(NA, nrow=28, ncol=436)
for (i in 1:436){
  for (j in 1:28){
    index_homo_par = apply(trio_genotype4[, , j, i], 2, FUN=get_index_homo_par)
    index_homo_par = which(index_homo_par)
    homo_par_num_shuff[j, i] = length(index_homo_par)
    trio_genotype5[, 1:length(index_homo_par), j, i] = trio_genotype4[, index_homo_par, j, i]
  }
}
sum(homo_par_num_shuff)
# [1] 10,836,901
# this number is diff from sum(homo_par_num) because some trio were discarded/added due to NA

homo_par_cor_off_num_shuff <- matrix(NA, nrow=28, ncol=436)
for (i in 1:436){
  for (j in 1:28){
    index_cor_off = 
      apply(trio_genotype5[, 1:homo_par_num_shuff[j, i], j, i], 2, FUN=get_index_homo_par_cor_off)
    index_cor_off = which(index_cor_off)
    homo_par_cor_off_num_shuff[j, i] = length(index_cor_off)
  }
}
sum(homo_par_cor_off_num_shuff)
# [1] 7,626,557
1 - 7626557/10836901
# [1] 0.2962419
# genotyping error rate with homozygous parental markers and shuffled offspring

saveRDS(trio_genotype4, file="get_genotype_errors/trio_genotype4.rds")
saveRDS(trio_genotype5, file="get_genotype_errors/trio_genotype5.rds")
saveRDS(homo_par_num_shuff, file="get_genotype_errors/homo_par_num_shuff.rds")
saveRDS(homo_par_cor_off_num_shuff, file="get_genotype_errors/homo_par_cor_off_num_shuff.rds")

# trio_genotype4 <- readRDS("get_genotype_errors/trio_genotype4.rds")
# trio_genotype5 <- readRDS("get_genotype_errors/trio_genotype5.rds")
# homo_par_num_shuff <- readRDS("get_genotype_errors/homo_par_num_shuff.rds")
# homo_par_cor_off_num_shuff <- readRDS("get_genotype_errors/homo_par_cor_off_num_shuff.rds")



# examine error rate with different combinations of homozygous parents (0,0; 0,2; 2,2)
# trio_genotype3 <- readRDS("codes/get_genotype_errors_trio_genotype3.rds")
trio_genotype3_00 <- array(NA, dim=dim(trio_genotype3))
trio_genotype3_02 <- array(NA, dim=dim(trio_genotype3))
trio_genotype3_22 <- array(NA, dim=dim(trio_genotype3))
homo_par_num_00 <- matrix(NA, nrow=28, ncol=436)
homo_par_num_02 <- matrix(NA, nrow=28, ncol=436)
homo_par_num_22 <- matrix(NA, nrow=28, ncol=436)

for (i in 1:436){
  # print(i)
  for (j in 1:28){
    # print(j)
    index_homo_par = apply(trio_genotype3[2:3, , j, i], 2, FUN=get_index_homo_par_sort)
    index_homo_par_00 = index_homo_par[1, ]
    index_homo_par_02 = index_homo_par[2, ]
    index_homo_par_22 = index_homo_par[3, ]
    index_homo_par_00 = which(index_homo_par_00)
    index_homo_par_02 = which(index_homo_par_02)
    index_homo_par_22 = which(index_homo_par_22)
    
    homo_par_num_00[j, i] = length(index_homo_par_00)
    homo_par_num_02[j, i] = length(index_homo_par_02)
    homo_par_num_22[j, i] = length(index_homo_par_22)
    
    if (length(index_homo_par_00)!=0){
      trio_genotype3_00[, 1:length(index_homo_par_00), j, i] = trio_genotype3[, index_homo_par_00, j, i]
    }
    if (length(index_homo_par_02)!=0){
      trio_genotype3_02[, 1:length(index_homo_par_02), j, i] = trio_genotype3[, index_homo_par_02, j, i]
    }
    if (length(index_homo_par_22)!=0){
      trio_genotype3_22[, 1:length(index_homo_par_22), j, i] = trio_genotype3[, index_homo_par_22, j, i]
    }
  }
}

homo_par_cor_off_num_00 <- matrix(NA, nrow=28, ncol=436)
homo_par_cor_off_num_02 <- matrix(NA, nrow=28, ncol=436)
homo_par_cor_off_num_22 <- matrix(NA, nrow=28, ncol=436)
for (i in 1:436){
  for (j in 1:28){
    index_cor_off = 
      apply(trio_genotype3_00[, 1:homo_par_num_00[j, i], j, i], 2, FUN=get_index_homo_par_cor_off)
    index_cor_off = which(index_cor_off)
    homo_par_cor_off_num_00[j, i] = length(index_cor_off)
  }
}
for (i in 1:436){
  for (j in 1:28){
    if (homo_par_num_02[j, i] != 0){
      index_cor_off = 
        apply(trio_genotype3_02[, 1:homo_par_num_02[j, i], j, i], 2, FUN=get_index_homo_par_cor_off)
      index_cor_off = which(index_cor_off)
      homo_par_cor_off_num_02[j, i] = length(index_cor_off)
    }
  }
}
for (i in 1:436){
  for (j in 1:28){
    index_cor_off = 
      apply(trio_genotype3_22[, 1:homo_par_num_22[j, i], j, i], 2, FUN=get_index_homo_par_cor_off)
    index_cor_off = which(index_cor_off)
    homo_par_cor_off_num_22[j, i] = length(index_cor_off)
  }
}

sum(homo_par_num_00)
# [1] 4659480
sum(homo_par_num_02)
# [1] 0
sum(homo_par_num_22)
# [1] 6186815
sum(homo_par_cor_off_num_00)
# [1] 3666338
sum(homo_par_cor_off_num_02, na.rm=T)
# [1] 0
sum(homo_par_cor_off_num_22)
# [1] 5041726
1 - 3666338 / 4659480
# [1] 0.2131444
# genotyping error rate with 0, 0 parents
1 - 5037199 / 6186815
# [1] 0.1858171
# genotyping error rate with 2, 2 parents

# there are no (0, 2) parental genotype
# for (i in 1:436){
#   for (j in 1:28){
#     print(length(which(apply(trio_genotype3[2:3, , j, i], 2, FUN=function(x){setequal(x, c(0, 2))}))))
#   }
# }

saveRDS(trio_genotype3_00, file="get_genotype_errors/trio_genotype3_00.rds")
saveRDS(trio_genotype3_02, file="get_genotype_errors/trio_genotype3_02.rds")
saveRDS(trio_genotype3_22, file="get_genotype_errors/trio_genotype3_22.rds")
saveRDS(homo_par_num_00, file="get_genotype_errors/homo_par_num_00.rds")
saveRDS(homo_par_num_02, file="get_genotype_errors/homo_par_num_02.rds")
saveRDS(homo_par_num_22, file="get_genotype_errors/homo_par_num_22.rds")
saveRDS(homo_par_cor_off_num_00, file="get_genotype_errors/homo_par_cor_off_num_00.rds")
saveRDS(homo_par_cor_off_num_02, file="get_genotype_errors/homo_par_cor_off_num_02.rds")
saveRDS(homo_par_cor_off_num_22, file="get_genotype_errors/homo_par_cor_off_num_22.rds")



# 
# trio_genotype3_00 <- readRDS("codes/get_genotype_errors_trio_genotype3_00.rds")
# trio_genotype3_02 <- readRDS("codes/get_genotype_errors_trio_genotype3_02.rds")
# trio_genotype3_22 <- readRDS("codes/get_genotype_errors_trio_genotype3_22.rds")
# 
# sum(trio_genotype3_00[1, , , ]==1, na.rm=T)
# # [1] 934246
# sum(trio_genotype3_00[1, , , ]==2, na.rm=T)
# # [1] 58030
# 934246 / (4656811 - 3664535)
# # [1] 0.9415183
# 58030 / (4656811 - 3664535)
# # [1] 0.05848171
# sum(trio_genotype3_22[1, , , ]==1, na.rm=T)
# # [1] 1083833
# sum(trio_genotype3_22[1, , , ]==0, na.rm=T)
# # [1] 60249
# 1083833 / (6181281 - 5037199)
# # [1] 0.9473386
# 60249 / (6181281 - 5037199)
# # [1] 0.05266144







