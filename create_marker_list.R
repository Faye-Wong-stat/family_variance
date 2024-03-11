setwd("~/family_variance/")
library(vcfR)
# this script sorts the "cleaned_data/marker_info_cM_parent.rds" 
# "cleaned_data/phased_marker_list.rds" 
# "cleaned_data/phased_marker_info.rds"
# "cleaned_data/phased_marker_matrix.rds"
# by chromosome and cM 



phased_marker <- read.vcfR("phased_data/phased.vcf")
phased_marker_info <- as.data.frame(phased_marker@fix)
phased_marker_matrix <- as.data.frame(phased_marker@gt)
phased_marker_info$FORMAT <- phased_marker_matrix$FORMAT
phased_marker_matrix <- phased_marker_matrix[, -1]
rownames(phased_marker_info) <- phased_marker_info$ID
rownames(phased_marker_matrix) <- phased_marker_info$ID

marker_info_cM_parent <- readRDS("get_markers_cM/marker_info_cM_parent.rds")
index_NA_cM <- is.na(marker_info_cM_parent$cM)
marker_info_cM_parent <- marker_info_cM_parent[!index_NA_cM, ]
phased_marker_info <- phased_marker_info[!index_NA_cM, ]
phased_marker_matrix <- phased_marker_matrix[!index_NA_cM, ]

index_order <- order(marker_info_cM_parent$CHROM, marker_info_cM_parent$cM)
marker_info_cM_parent <- marker_info_cM_parent[index_order, ]
phased_marker_info <- phased_marker_info[index_order, ]
phased_marker_matrix <- phased_marker_matrix[index_order, ]

indiv_names <- colnames(phased_marker_matrix)
chrom_names <- unique(marker_info_cM_parent[,"CHROM"])

dim(marker_info_cM_parent)
# [1] 34911    10



for (j in 1:length(chrom_names)){
  index_duplicated = duplicated(marker_info_cM_parent[marker_info_cM_parent$CHROM==chrom_names[j], ]$cM)
  print(which(index_duplicated))
  
  if (sum(index_duplicated)!=0){
    marker_info_cM_parent[marker_info_cM_parent$CHROM==chrom_names[j], ][index_duplicated, ] = NA
    phased_marker_matrix[phased_marker_info$CHROM==chrom_names[j], ][index_duplicated, ] = NA
    phased_marker_info[phased_marker_info$CHROM==chrom_names[j], ][index_duplicated, ] = NA
  }
}

phased_marker_matrix <- phased_marker_matrix[!is.na(marker_info_cM_parent$cM), ]
phased_marker_info <- phased_marker_info[!is.na(marker_info_cM_parent$cM), ]
marker_info_cM_parent <- marker_info_cM_parent[!is.na(marker_info_cM_parent$cM), ]

dim(marker_info_cM_parent)
# [1] 34911    10
dim(phased_marker_info)
# [1] 34911     9
dim(phased_marker_matrix)
# [1] 34911   1035



marker_list <- vector(mode="list", length=length(chrom_names))
for (j in 1:length(chrom_names)){
  marker_list[[j]] = array(0, dim=c(sum(phased_marker_info$CHROM==chrom_names[j]), 
                                    2, 
                                    length(indiv_names)))
}

for (j in 1:length(chrom_names)){
  for (i in 1:length(indiv_names)){
    indiv_chrom = phased_marker_matrix[phased_marker_info$CHROM==chrom_names[j], indiv_names[i]]
    indiv_chrom_A = gsub("\\|[01]", "", indiv_chrom)
    indiv_chrom_B = gsub("[01]\\|", "", indiv_chrom)
    marker_list[[j]][indiv_chrom_A=="1", 1, i] = 1
    marker_list[[j]][indiv_chrom_B=="1", 2, i] = 1
  }
}

saveRDS(marker_info_cM_parent, "create_marker_list/marker_info_cM_parent.rds")
saveRDS(marker_list, "create_marker_list/phased_marker_list.rds")
saveRDS(phased_marker_info, "create_marker_list/phased_marker_info.rds")
saveRDS(phased_marker_matrix, "create_marker_list/phased_marker_matrix.rds")
saveRDS(indiv_names, "create_marker_list/indiv_names.rds")
saveRDS(chrom_names, "create_marker_list/chrom_names.rds")
























