setwd("family_variance/")



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



marker_info <- read.csv("generate_vcffiles/marker_info.csv", row.names = 1)
marker_matrix <- read.csv("generate_vcffiles/marker_matrix.csv", row.names = 1)
# marker_info_valid <- read.csv("cleaned_data/marker_info_valid.csv", row.names=1, check.names=F)
# marker_matrix_valid <- read.csv("cleaned_data/marker_matrix_valid.csv", row.names=1, check.names=F)
# marker_info_parent <- read.csv("cleaned_data/marker_info_parent.csv", row.names=1, check.names=F)
# marker_matrix_parent <- read.csv("cleaned_data/marker_matrix_parent.csv", row.names=1, check.names=F)
pedigree_valid <- read.csv("generate_vcffiles/pedigree_valid.csv")
dim(marker_info)
# [1] 37285     9
dim(marker_matrix)
# [1] 37285   999
dim(pedigree_valid)
# [1] 436   3



off_names <- pedigree_valid$Accession_ID
par1_names <- pedigree_valid$Integrated_P1
par2_names <- pedigree_valid$Integrated_P2
chrom_names <- unique(marker_info[,"CHROM"])



marker_matrix2 <- matrix(NA, nrow=nrow(marker_matrix), ncol=ncol(marker_matrix))
marker_matrix2[marker_matrix=="0/0"] <- 0
marker_matrix2[marker_matrix=="0/1"] <- 1
marker_matrix2[marker_matrix=="1/1"] <- 2

colnames(marker_matrix2) <- colnames(marker_matrix)
rownames(marker_matrix2) <- rownames(marker_matrix)



genotype_error <- array(NA, dim=c(4, 2000, length(chrom_names), length(off_names)))
for (i in 1:length(off_names)){
  # print(i)
  off_geno = marker_matrix2[, colnames(marker_matrix2)==off_names[i]]
  par1_geno = marker_matrix2[, colnames(marker_matrix2)==par1_names[i]]
  par2_geno = marker_matrix2[, colnames(marker_matrix2)==par2_names[i]]
  
  for (j in 1:length(chrom_names)){
    # print(j)
    chrom_pos = which(marker_info$CHROM == chrom_names[j])
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
        genotype_error[1, 1:length(error_pos), j, i] = error_pos
        genotype_error[2, 1:length(error_pos), j, i] = off_chrom[error_pos]
        genotype_error[3:4, 1:length(error_pos), j, i] = trio_chrom[1:2, error_pos]
      }
    }
  }
}



saveRDS(genotype_error, file="genotype_error/genotype_error.rds")



{
# marker_info_valid[1:10, ]
# #       CHROM    POS ID REF ALT QUAL FILTER INFO FORMAT
# # 29785    1A   9047 NA   C   N   99   PASS   NA     GT
# # 36850    1A  10255 NA   G   N   99   PASS   NA     GT
# # 11745    1A  64241 NA   G   N   99   PASS   NA     GT
# # 579      1A  83188 NA   C   N   99   PASS   NA     GT
# # 15282    1A  88065 NA   T   N   99   PASS   NA     GT
# # 24893    1A 101233 NA   C   N   99   PASS   NA     GT
# # 28432    1A 115766 NA   A   N   99   PASS   NA     GT
# # 11353    1A 150544 NA   T   N   99   PASS   NA     GT
# # 14802    1A 153468 NA   G   N   99   PASS   NA     GT
# # 6236     1A 182814 NA   A   N   99   PASS   NA     GT
# marker_matrix_valid[1:10, 1:5]
# #              MSU73 05C092P011 05C120P009 06C132P002 08C013P605
# # AX.184660696   0/0        1/1        1/1        1/1        1/1
# # AX.184955090   1/1        1/1        1/1        1/1        0/1
# # AX.184194282   1/1        1/1        1/1        1/1        0/1
# # AX.123361821   0/0        1/1        1/1        1/1        0/1
# # AX.184279290   0/0        0/0        0/0        0/0        0/1
# # AX.184521015   1/1        1/1        1/1        1/1        0/1
# # AX.184617406   0/0        0/0        0/0        0/0        0/1
# # AX.184184503   1/1        0/0        0/0        0/0        0/1
# # AX.184268232   1/1        1/1        1/1        1/1        1/1
# # AX.184066211   0/0        0/0        0/0        0/0        0/1
# marker_matrix_valid2[1:10, 1:5]
# #              MSU73 05C092P011 05C120P009 06C132P002 08C013P605
# # AX.184660696     0          2          2          2          2
# # AX.184955090     2          2          2          2          1
# # AX.184194282     2          2          2          2          1
# # AX.123361821     0          2          2          2          1
# # AX.184279290     0          0          0          0          1
# # AX.184521015     2          2          2          2          1
# # AX.184617406     0          0          0          0          1
# # AX.184184503     2          0          0          0          1
# # AX.184268232     2          2          2          2          2
# # AX.184066211     0          0          0          0          1
# marker_info_parent[1:10, ]
# #       CHROM    POS ID REF ALT QUAL FILTER INFO FORMAT
# # 29785    1A   9047 NA   C   N   99   PASS   NA     GT
# # 36850    1A  10255 NA   G   N   99   PASS   NA     GT
# # 11745    1A  64241 NA   G   N   99   PASS   NA     GT
# # 579      1A  83188 NA   C   N   99   PASS   NA     GT
# # 15282    1A  88065 NA   T   N   99   PASS   NA     GT
# # 24893    1A 101233 NA   C   N   99   PASS   NA     GT
# # 28432    1A 115766 NA   A   N   99   PASS   NA     GT
# # 11353    1A 150544 NA   T   N   99   PASS   NA     GT
# # 14802    1A 153468 NA   G   N   99   PASS   NA     GT
# # 6236     1A 182814 NA   A   N   99   PASS   NA     GT
# marker_matrix_parent[1:10, 1:5]
# #              Akashi005 BeaverSweet Everbearing185 Hayazaki  K1
# # AX.184660696       0/1         0/1            0/1      0/0 0/0
# # AX.184955090       0/1         1/1            0/1      1/1 1/1
# # AX.184194282       0/0         0/1            0/1      1/1 1/1
# # AX.123361821       0/0         0/0            0/1      0/0 0/1
# # AX.184279290       0/1         0/0            0/1      0/0 0/0
# # AX.184521015       0/1         0/1            0/1      1/1 1/1
# # AX.184617406       0/1         0/0            0/1      0/0 0/0
# # AX.184184503       1/1         1/1            1/1      1/1 1/1
# # AX.184268232       1/1         0/1            1/1      1/1 1/1
# # AX.184066211       0/1         0/1            0/1      0/0 0/0
# marker_matrix_parent2[1:10, 1:5]
# #              Akashi005 BeaverSweet Everbearing185 Hayazaki K1
# # AX.184660696         1           1              1        0  0
# # AX.184955090         1           2              1        2  2
# # AX.184194282         0           1              1        2  2
# # AX.123361821         0           0              1        0  1
# # AX.184279290         1           0              1        0  0
# # AX.184521015         1           1              1        2  2
# # AX.184617406         1           0              1        0  0
# # AX.184184503         2           2              2        2  2
# # AX.184268232         2           1              2        2  2
# # AX.184066211         1           1              1        0  0
}




