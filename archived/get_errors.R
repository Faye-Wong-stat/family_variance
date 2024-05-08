setwd("~/family_variance/")



# get_index_homoz_par_sort <- function(vec){
#   # vec is a vector with length 3
#   # vec = c(p1, p2, off)
#   # it the the genotype (0, 1, or 2) of two parents and their offspring
#   if (!any(is.na(vec))){
#     if (vec[1]==0 & vec[2]==0){
#       return(c(T, F, F))
#     } else if ( setequal(vec[1:2], c(0, 2)) ){
#       return(c(F, T, F))
#     } else if (vec[1]==2 & vec[2]==2){
#       return(c(F, F, T))
#     } else{
#       return(c(F, F, F))
#     }
#   } else{
#     return(c(NA, NA, NA))
#   }
# }
# get_index_homoz_par_cor_off <- function(vec){
#   # vec is a vector with length 3
#   # vec = c(p1, p2, off)
#   # it the the genotype (0, 1, or 2) of two parents and their offspring
#   if (!any(is.na(vec))){
#     if (vec[3]==1){
#       return(T)
#     } else{
#       return(F)
#     }
#   } else{
#     return(NA)
#   }
# }



pedigree_valid <- readRDS("generate_vcffiles/pedigree_valid.rds")
marker_matrix <- readRDS("generate_vcffiles/marker_matrix.rds")

par1_names <- pedigree_valid$P1
par2_names <- pedigree_valid$P2
off_names <- pedigree_valid$ID

marker_matrix2 <- matrix(NA, nrow=nrow(marker_matrix), ncol=ncol(marker_matrix))
marker_matrix2[marker_matrix=="0/0"] <- 0
marker_matrix2[marker_matrix=="0/1"] <- 1
marker_matrix2[marker_matrix=="1/1"] <- 2
rownames(marker_matrix2) <- rownames(marker_matrix)
colnames(marker_matrix2) <- colnames(marker_matrix)



trio_geno <- array(NA, dim=c(nrow(marker_matrix2), 3, length(off_names)))
for (i in 1:length(off_names)){
  par1_geno = marker_matrix2[, colnames(marker_matrix2)==par1_names[i]]
  par2_geno = marker_matrix2[, colnames(marker_matrix2)==par2_names[i]]
  off_geno = marker_matrix2[, colnames(marker_matrix2)==off_names[i]]
  trio_geno[, , i] = cbind(par1_geno, par2_geno, off_geno)
}
dim(trio_geno)
# [1] 37282     3   593
saveRDS(trio_geno, "get_errors/trio_geno.rds")
# trio_geno <- readRDS("get_errors/trio_geno.rds")

# A <- trio_geno[1:5, , 1:4]
# A[3, 3, 1] <- 1

error_rate <- apply(trio_geno, c(1, 3), FUN=function(x){
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

saveRDS(error_rate, "get_errors/error_rate.rds")

# index_homoz_par_sort <- apply(trio_geno, c(1, 3), FUN=get_index_homoz_par_sort)
# # sum(index_homoz_par_sort[2, ,], na.rm=T)
# # [1] 794963
# # 794963 markers across 417 trios have homozygous parental genotype 02
# index_homoz_par_02 <- index_homoz_par_sort[2, , ]
# # dim(index_homoz_par_02)
# # [1] 37282   417
# dim(index_homoz_par_02) <- c(dim(index_homoz_par_02)[1], 1, dim(index_homoz_par_02)[2])
# index_homoz_par_02 <- index_homoz_par_02[, rep(1,3), ]
# # dim(index_homoz_par_02)
# # [1] 37282     3   417
# 
# trio_geno_02 <- array(NA, dim=dim(trio_geno))
# trio_geno_02[which(index_homoz_par_02)] <- trio_geno[which(index_homoz_par_02)]
# # sum(!is.na(trio_geno_02[, 3, ]))
# # [1] 794963
# # matches with 
# # sum(index_homoz_par_sort[2, ,], na.rm=T)
# 
# error_rate_02 <- apply(trio_geno_02, c(1, 3), FUN=get_index_homoz_par_cor_off)
# error_rate_02_means <- apply(error_rate_02, 2, mean, na.rm=T)
# pdf("get_errors/plots/p1.pdf")
# hist(error_rate_02_means, xlab="percentage of correct offspring genotype per trio")
# # plot(density(na.omit(error_rate_02_means)))
# dev.off()

mean(error_rate, na.rm=T)
# [1] 0.9764082
error_rate_means <- apply(error_rate, 2, mean, na.rm=T)
pdf("get_errors/plots/p1.pdf")
hist(error_rate_means, xlab="percentage of correct offspring genotype per trio")
# plot(density(na.omit(error_rate_02_means)))
dev.off()



index_off_low_error <- apply(error_rate, 2, function(x){
  ifelse(mean(x, na.rm=T) > 0.95, T, F)
})
sum(index_off_low_error, na.rm = T)
# [1] 541
index_off_low_error_90 <- apply(error_rate, 2, function(x){
  ifelse(mean(x, na.rm=T) > 0.9, T, F)
})
length(which(index_off_low_error_90))
# [1] 562

# saveRDS(index_off_low_error, "get_errors/index_off_low_error.rds")
# index_off_low_error <- readRDS("get_errors/index_off_low_error.rds")

# error_rate_02_exam <- apply(trio_geno_02[, , which(index_off_low_error)],
#                        c(1, 3),
#                        FUN=get_index_homoz_par_cor_off)
# error_rate_02_exam_means <- apply(error_rate_02_exam, 2, mean, na.rm=T)
# mean(error_rate_02_exam_means)
# # [1] 0.9686263
# 1 - 0.9686263
# # [1] 0.0313737
# pdf("get_errors/plots/p2.pdf")
# hist(error_rate_02_exam_means, xlab="percentage of correct offspring genotype per trio")
# dev.off()




