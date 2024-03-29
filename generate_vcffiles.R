setwd("~/family_variance/")
library(vcfR)



# load data
marker <- readRDS("remove_end_progenies/marker.rds")
marker_info <- readRDS("remove_end_progenies/marker_info.rds")
# genotype <- read.csv("new_provided_data/Factorial_50K_SNP_Array_Genotypes.csv")
pedigree_parent <- readRDS("remove_end_progenies/pedigree_parents.rds")
pedigree_valid <- readRDS("remove_end_progenies/pedigree_valid.rds")
pedigree_genotyped <- readRDS("remove_end_progenies/pedigree_genotyped.rds")



# function
# check if there are any special characters in individual names 
any_spec_char <- function(spec_char, x){
  return(any(grepl(spec_char, x)))
}
# replace "_" and "-" with "" in vectors
replace_spec_char <- function(x){
  return(gsub("-", "", gsub("_", "", x)))
}
{
# any_spec_char("_", marker$X)
# # TRUE
# any_spec_char("\\.", marker$X)
# # FALSE
# any_spec_char("-", marker$X)
# # TRUE
# # any_spec_char("", marker$X)
# length(unique(marker$X))
# # 999
# length(unique( gsub("-", "", gsub("_", "", marker$X)) ))
# # 999
}
# replace "_" and "-" with "" in marker matrix, parents pedigree, and validation pedigree
marker_clean <- marker
rownames(marker_clean) <- replace_spec_char(rownames(marker_clean))
pedigree_parent_clean <- pedigree_parent
pedigree_parent_clean <- as.data.frame(apply(pedigree_parent_clean, 2, FUN=replace_spec_char))
pedigree_valid_clean <- pedigree_valid
pedigree_valid_clean <- as.data.frame(apply(pedigree_valid_clean, 2, FUN=replace_spec_char))
pedigree_clean <- pedigree_genotyped
pedigree_clean <- as.data.frame(apply(pedigree_clean, 2, FUN=replace_spec_char))
# # change "-" into "." in SNP ID in marker_info_clean
marker_info_clean <- marker_info
# marker_info_clean$probe_id <- gsub("-", ".", marker_info_clean$probe_id)

sum(!marker_info_clean$probe_id %in% colnames(marker_clean))
# [1] 0
apply(pedigree_clean, 2, FUN=function(x){
  sum(! x %in% rownames(marker_clean))
})
# ID  P1  P2 
# 0 282 270
apply(pedigree_valid_clean, 2, FUN=function(x){
  sum(! x %in% rownames(marker_clean))
})
# ID P1 P2 
# 0  0  0
apply(pedigree_parent_clean, 2, FUN=function(x){
  sum(! x %in% rownames(marker_clean))
})
# ID  P1  P2 
# 0 221 270
dim(marker_clean)
# [1]   1035 37441
dim(marker_info_clean)
# [1] 37444    13
dim(pedigree_clean)
# [1] 1007   3
dim(pedigree_valid_clean)
# [1] 417   3
dim(pedigree_parent_clean)
# [1] 590   3



# pedigree_clean$Accession_ID[!pedigree_clean$Accession_ID %in% marker_clean$X]



# generate vcf file metadata section
metaData = c(
  "##fileformat=VCFv4.3",
  "##source=StrawberryOutcross",
  "##phasing=none",
  "##FILTER=<ID=PASS,Description=PASS>",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)



# marker matrix for parent individuals 
# parent here means testing set
marker_clean_parent <- marker_clean[rownames(marker_clean) %in% pedigree_parent_clean$ID, ]
marker_clean_parent <- t(marker_clean_parent)
marker_clean_parent <- as.data.frame(marker_clean_parent)
marker_clean_parent[marker_clean_parent==2] <- "1/1"
marker_clean_parent[marker_clean_parent==1] <- "0/1"
marker_clean_parent[marker_clean_parent==0] <- "0/0"

# pdf("generate_vcffiles/plots/hist of marker quality.pdf")
# hist(marker_info_clean$hit_score)
# dev.off()

# create marker info matrix for parent individuals
marker_info_parent <- 
  cbind(marker_info_clean[match(rownames(marker_clean_parent), marker_info_clean$probe_id), 
                          c("Chromosome", "ref_site", "probe_id", "ref_nt")], 
        rep("N", nrow(marker_clean_parent)), 
        rep(99, nrow(marker_clean_parent)), 
        rep("PASS", nrow(marker_clean_parent)), 
        rep(NA, nrow(marker_clean_parent)), 
        rep("GT", nrow(marker_clean_parent)))
colnames(marker_info_parent)[1:9] <- 
  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
nrow(marker_info_parent)
# 37441

# remove markers with NA position
na_position_parent <- which(is.na(marker_info_parent$POS))
# remove markers with "-" for ref allele
na_ref_parent <- grep("[^ATCG]", marker_info_parent$REF)
# remove duplicated markers 
dup <- which(duplicated(marker_info_parent[, c("CHROM", "POS")]))
# create a union of all the markers need to be removed
remove_parent <- union(union(na_position_parent, na_ref_parent), dup)
length(remove_parent)
# 159
marker_info_parent <- marker_info_parent[-remove_parent, ]

# reorder the markers by chromosome first then position
order_parent <- order(marker_info_parent$CHROM, marker_info_parent$POS)
marker_info_parent <- marker_info_parent[order_parent, ]
rownames(marker_info_parent) <- marker_info_parent$ID

marker_clean_parent <- marker_clean_parent[-remove_parent, ]
marker_clean_parent <- marker_clean_parent[order_parent, ]
# make sure two data sets are in the right order
sum(rownames(marker_info_parent)!=rownames(marker_clean_parent))
# [1] 0

saveRDS(marker_info_parent, "generate_vcffiles/marker_info_parent.rds")
saveRDS(marker_clean_parent, "generate_vcffiles/marker_matrix_parent.rds")



marker_info_parent <- as.matrix(marker_info_parent, row.names=NULL)

# remove " " empty space before the positions
marker_info_parent[,2] <- gsub(" ", "", marker_info_parent[,2])

marker_clean_parent <- as.matrix(marker_clean_parent, row.names=NULL)

# create a "vcfR" object and save as vcf file
marker_parent <- new("vcfR", meta=metaData, fix=marker_info_parent, gt=marker_clean_parent)
write.vcf(marker_parent, "generate_vcffiles/genotype_parent.vcf.gz")




# marker matrix for validation individuals 
marker_clean_valid <- marker_clean[rownames(marker_clean) %in% pedigree_valid_clean$ID, ]
marker_clean_valid <- t(marker_clean_valid)
marker_clean_valid <- as.data.frame(marker_clean_valid)
marker_clean_valid[marker_clean_valid==2] <- "1/1"
marker_clean_valid[marker_clean_valid==1] <- "0/1"
marker_clean_valid[marker_clean_valid==0] <- "0/0"

# # create marker info matrix for valid individuals
# marker_info_valid <-
#   cbind(marker_info_clean[match(rownames(marker_clean_valid), marker_info_clean$probe_id), 
#                           c("Chromosome", "ref_site", "probe_id", "ref_nt")],
#         rep("N", nrow(marker_clean_valid)),
#         rep(99, nrow(marker_clean_valid)),
#         rep("PASS", nrow(marker_clean_valid)),
#         rep(NA, nrow(marker_clean_valid)),
#         rep("GT", nrow(marker_clean_valid)))
# colnames(marker_info_valid)[1:9] <-
#   c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
# nrow(marker_info_valid)
# # 37444
# 
# # remove markers with NA position
# na_position_valid <- which(is.na(marker_info_valid$POS))
# # remove markers with "-" for ref allele
# na_ref_valid <- grep("[^ATCG]", marker_info_valid$REF)
# # remove duplicated markers
# dup <- which(duplicated(marker_info_valid[, c("CHROM", "POS")]))
# # create a union of all the markers need to be removed
# remove_valid <- union(union(na_position_valid, na_ref_valid), dup)
# length(remove_valid)
# # 159
# marker_info_valid <- marker_info_valid[-remove_valid, ]
# 
# # reorder the markers by chromosome first then position
# order_valid <- order(marker_info_valid$CHROM, marker_info_valid$POS)
# marker_info_valid <- marker_info_valid[order_valid, ]
# rownames(marker_info_valid) <- marker_info_valid$ID

marker_clean_valid <- marker_clean_valid[-remove_parent, ]
marker_clean_valid <- marker_clean_valid[order_parent, ]
# make sure two data sets are in the right order
sum(rownames(marker_info_parent)!=rownames(marker_clean_valid))
# [1] 0
sum(rownames(marker_clean_parent)!=rownames(marker_clean_valid))
# [1] 0

saveRDS(as.data.frame(marker_info_parent), "generate_vcffiles/marker_info_valid.rds")
saveRDS(marker_clean_valid, "generate_vcffiles/marker_matrix_valid.rds")



# marker_info_valid <- as.matrix(marker_info_valid, row.names=NULL)
# 
# # remove " " empty space before the positions
# marker_info_valid[,2] <- gsub(" ", "", marker_info_valid[,2])

marker_clean_valid <- as.matrix(marker_clean_valid, row.names=NULL)

# create a "vcfR" object and save as vcf file
marker_valid <- new("vcfR", meta=metaData, fix=marker_info_parent, gt=marker_clean_valid)
write.vcf(marker_valid, "generate_vcffiles/genotype_valid.vcf.gz")



# marker matrix for all individuals 
marker_clean <- t(marker_clean)
marker_clean <- as.data.frame(marker_clean)
marker_clean[marker_clean==2] <- "1/1"
marker_clean[marker_clean==1] <- "0/1"
marker_clean[marker_clean==0] <- "0/0"

marker_clean <- marker_clean[-remove_parent, ]
marker_clean <- marker_clean[order_parent, ]
sum(rownames(marker_info_parent)!=rownames(marker_clean))
# [1] 0

saveRDS(as.data.frame(marker_info_parent), "generate_vcffiles/marker_info.rds")
saveRDS(as.data.frame(marker_clean), "generate_vcffiles/marker_matrix.rds")



marker_clean <- as.matrix(marker_clean, row.names=NULL)

# rownames(marker_info_valid) <- rownames(marker_clean)

# create a "vcfR" object and save as vcf file
marker_all <- new("vcfR", meta=metaData, fix=marker_info_parent, gt=marker_clean)
write.vcf(marker_all, "generate_vcffiles/genotype.vcf.gz")



saveRDS(pedigree_clean, "generate_vcffiles/pedigree.rds")
saveRDS(pedigree_parent_clean, "generate_vcffiles/pedigree_parents.rds")
saveRDS(pedigree_valid_clean, "generate_vcffiles/pedigree_valid.rds")
