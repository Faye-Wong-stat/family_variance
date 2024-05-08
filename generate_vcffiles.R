setwd("~/family_variance/")
library(vcfR)



# load data
marker <- readRDS("remove_end_progenies/marker.rds")
marker_info <- readRDS("remove_end_progenies/marker_info.rds")
pedigree_parent <- readRDS("remove_end_progenies/pedigree_parents.rds")
pedigree_valid <- readRDS("remove_end_progenies/pedigree_valid.rds")
pedigree_genotyped <- readRDS("remove_end_progenies/pedigree_genotyped.rds")



# functions
# check if there are any special characters in individual names 
any_spec_char <- function(spec_char, x){
  return(any(grepl(spec_char, x)))
}
# replace "_" and "-" with "" in vectors
replace_spec_char <- function(x){
  return(gsub("-", "", gsub("_", "", x)))
}

# check if there are any special characters in individual names 
any_spec_char("_", rownames(marker))
# TRUE
any_spec_char(".", rownames(marker))
# TRUE
any_spec_char("-", rownames(marker))
# TRUE
length(unique(rownames(marker)))
# 1007
length(unique( gsub("-", "", gsub("_", "", rownames(marker))) ))
# 1007

# replace "_" and "-" with "" in marker matrix, parents pedigree, and validation pedigree
marker_clean <- marker
rownames(marker_clean) <- replace_spec_char(rownames(marker_clean))
pedigree_parent_clean <- pedigree_parent
pedigree_parent_clean <- as.data.frame(apply(pedigree_parent_clean, 2, FUN=replace_spec_char))
pedigree_valid_clean <- pedigree_valid
pedigree_valid_clean <- as.data.frame(apply(pedigree_valid_clean, 2, FUN=replace_spec_char))
pedigree_clean <- pedigree_genotyped
pedigree_clean <- as.data.frame(apply(pedigree_clean, 2, FUN=replace_spec_char))

marker_info_clean <- marker_info
# pdf("generate_vcffiles/plots/hist of marker quality.pdf")
# hist(marker_info_clean$hit_score)
# dev.off()
marker_clean <- t(marker_clean)
marker_clean <- as.data.frame(marker_clean)

sum(!marker_info_clean$probe_id %in% rownames(marker_clean))
# [1] 0
apply(pedigree_clean, 2, FUN=function(x){
  sum(! x %in% colnames(marker_clean))
})
# ID  P1  P2 
# 0 221 270
apply(pedigree_valid_clean, 2, FUN=function(x){
  sum(! x %in% colnames(marker_clean))
})
# ID P1 P2 
# 0  0  0
apply(pedigree_parent_clean, 2, FUN=function(x){
  sum(! x %in% colnames(marker_clean))
})
# ID  P1  P2 
# 0 221 270
dim(marker_clean)
# [1] 37441  1007
dim(marker_info_clean)
# [1] 37441    13
dim(pedigree_clean)
# [1] 1007   3
dim(pedigree_valid_clean)
# [1] 593   3
dim(pedigree_parent_clean)
# [1] 414   3





# generate vcf file metadata section
metaData = c(
  "##fileformat=VCFv4.3",
  "##source=StrawberryOutcross",
  "##phasing=none",
  "##FILTER=<ID=PASS,Description=PASS>",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)





# create marker info matrix for all individuals

marker_info <- 
  cbind(marker_info_clean[match(rownames(marker_clean), marker_info_clean$probe_id), 
                          c("Chromosome", "ref_site", "probe_id", "ref_nt")], 
        rep("N", nrow(marker_clean)), 
        rep(99, nrow(marker_clean)), 
        rep("PASS", nrow(marker_clean)), 
        rep(NA, nrow(marker_clean)), 
        rep("GT", nrow(marker_clean)))
colnames(marker_info)[1:9] <- 
  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
dim(marker_info)
# [1] 37441     9

# remove markers with NA position
na_position <- which(is.na(marker_info$POS))
# remove markers with "-" for ref allele
na_ref <- grep("[^ATCG]", marker_info$REF)
# remove duplicated markers 
dup <- which(duplicated(marker_info[, c("CHROM", "POS")]))
# create a union of all the markers need to be removed
marker_remove <- union(union(na_position, na_ref), dup)
length(marker_remove)
# 159
marker_info <- marker_info[-marker_remove, ]

# reorder the markers by chromosome first then position
marker_order <- order(marker_info$CHROM, marker_info$POS)
marker_info <- marker_info[marker_order, ]
rownames(marker_info) <- marker_info$ID

# create marker matrix for all individuals
marker_clean <- marker_clean[-marker_remove, ]
marker_clean <- marker_clean[marker_order, ]
sum(rownames(marker_info)!=rownames(marker_clean))
# [1] 0

marker_info$POS <- gsub(" ", "", marker_info$POS)

saveRDS(as.data.frame(marker_info), "generate_vcffiles/marker_info.rds")
saveRDS(as.data.frame(marker_clean), "generate_vcffiles/marker_matrix.rds")





# examine the genotyping error and select trios that are more likely to be real trios
par1_names <- pedigree_valid_clean$P1
par2_names <- pedigree_valid_clean$P2
off_names <- pedigree_valid_clean$ID

marker_matrix2 <- as.matrix(marker_clean)

trio_geno <- array(NA, dim=c(nrow(marker_matrix2), 3, length(off_names)))
for (i in 1:length(off_names)){
  par1_geno = marker_matrix2[, colnames(marker_matrix2)==par1_names[i]]
  par2_geno = marker_matrix2[, colnames(marker_matrix2)==par2_names[i]]
  off_geno = marker_matrix2[, colnames(marker_matrix2)==off_names[i]]
  trio_geno[, , i] = cbind(par1_geno, par2_geno, off_geno)
}
dim(trio_geno)
# [1] 37282     3   593
saveRDS(trio_geno, "generate_vcffiles/trio_geno.rds")

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

saveRDS(error_rate, "generate_vcffiles/error_rate.rds")

mean(error_rate, na.rm=T)
# [1] 0.9764082
error_rate_means <- apply(error_rate, 2, mean, na.rm=T)
pdf("generate_vcffiles/plots/p1.pdf")
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
# index_off_low_error is the logic of offspring of possible true trios

off_names <- off_names[index_off_low_error]

pedigree_valid_clean <- pedigree_clean[pedigree_clean$ID %in% off_names, ]
pedigree_parent_clean <- pedigree_clean[!pedigree_clean$ID %in% off_names, ]
dim(pedigree_valid_clean)
# [1] 541   3
dim(pedigree_parent_clean)
# [1] 466   3





marker_clean <- as.data.frame(marker_clean)
marker_clean[marker_clean==2] <- "1/1"
marker_clean[marker_clean==1] <- "0/1"
marker_clean[marker_clean==0] <- "0/0"

marker_clean <- as.matrix(marker_clean, row.names=NULL)

# create a "vcfR" object and save as vcf file
marker_all <- new("vcfR", meta=metaData, fix=as.matrix(marker_info), gt=marker_clean)
write.vcf(marker_all, "generate_vcffiles/genotype.vcf.gz")





# marker matrix for parent individuals 
# parent here means testing set
marker_clean_parent <- marker_clean[, colnames(marker_clean) %in% pedigree_parent_clean$ID]
marker_clean_parent[marker_clean_parent==2] <- "1/1"
marker_clean_parent[marker_clean_parent==1] <- "0/1"
marker_clean_parent[marker_clean_parent==0] <- "0/0"

sum(rownames(marker_info)!=rownames(marker_clean_parent))
# [1] 0

saveRDS(marker_info, "generate_vcffiles/marker_info_parent.rds")
saveRDS(marker_clean_parent, "generate_vcffiles/marker_matrix_parent.rds")
# create a "vcfR" object and save as vcf file
marker_parent <- new("vcfR", meta=metaData, fix=as.matrix(marker_info), gt=marker_clean_parent)
write.vcf(marker_parent, "generate_vcffiles/genotype_parent.vcf.gz")



# marker matrix for validation individuals 
marker_clean_valid <- marker_clean[, colnames(marker_clean) %in% pedigree_valid_clean$ID]
marker_clean_valid[marker_clean_valid==2] <- "1/1"
marker_clean_valid[marker_clean_valid==1] <- "0/1"
marker_clean_valid[marker_clean_valid==0] <- "0/0"

sum(rownames(marker_info)!=rownames(marker_clean_valid))
# [1] 0

saveRDS(marker_info, "generate_vcffiles/marker_info_valid.rds")
saveRDS(marker_clean_valid, "generate_vcffiles/marker_matrix_valid.rds")
# create a "vcfR" object and save as vcf file
marker_valid <- new("vcfR", meta=metaData, fix=as.matrix(marker_info), gt=marker_clean_valid)
write.vcf(marker_valid, "generate_vcffiles/genotype_valid.vcf.gz")





saveRDS(pedigree_clean, "generate_vcffiles/pedigree.rds")
saveRDS(pedigree_parent_clean, "generate_vcffiles/pedigree_parents.rds")
saveRDS(pedigree_valid_clean, "generate_vcffiles/pedigree_valid.rds")
