setwd("~/family_variance/")



# # data tidying, skip this part to "# start from here"
# pedigree <- read.csv("provided_data/pedigree.csv", row.names = 1)
# pedigree_updated <- read.csv("provided_data/PedValidation.csv", row.names = 1)
# marker <- read.csv("provided_data/MarkerMatrix.csv", row.names = 1)
# marker_info <- read.csv("provided_data/marker_info.csv", row.names = 1)
# marker_progeny <- read.csv("provided_data/Progeny_geno.csv", row.names = 1, check.names = FALSE)
# marker_parent <- read.csv("provided_data/Parent_geno.csv", row.names = 1, check.names = FALSE)
# 
# colnames(pedigree) <- c("ID", "P1", "P2")
# marker_progeny <- t(marker_progeny)
# marker_parent <- t(marker_parent)
# colnames(marker_parent) <- gsub("-", ".", colnames(marker_parent))
# colnames(marker_progeny) <- gsub("-", ".", colnames(marker_progeny))
# marker_info$probe_id <- gsub("-", ".", marker_info$probe_id)
# rownames(marker_info) <- marker_info$probe_id
# 
# 
# 
# # clean up the updated pedigree
# dim(pedigree_updated)
# # [1] 535   6
# # remove individuals with parent "PI551753"
# pedigree_updated <- pedigree_updated[pedigree_updated$P1 != "PI551753" & 
#                                        pedigree_updated$P2 != "PI551753" , ]
# dim(pedigree_updated)
# # [1] 475   6
# # remove when c("MotherID", "FatherID") doesn't match c("P1", "P2")
# for (i in 1:nrow(pedigree_updated)){
#   if (! setequal(pedigree_updated[i, c("MotherID", "FatherID")], pedigree_updated[i, c("P1", "P2")])){
#     # print(i)
#     pedigree_updated[i, ] = NA
#   }
# }
# pedigree_updated <- pedigree_updated[!is.na(pedigree_updated$ID), ]
# dim(pedigree_updated)
# # [1] 417   6
# 
# # replace the pedigree data with updated pedigree
# dim(pedigree)
# # [1] 11651     3
# sum(pedigree_updated$ID %in% pedigree$ID)
# # [1] 417
# pedigree[match(pedigree_updated$ID, pedigree$ID), ] <- pedigree_updated[, c(1, 5:6)]
# 
# pedigree_names <- union(union(pedigree_updated$P1, pedigree_updated$P2), pedigree_updated$ID)
# length(pedigree_names)
# # [1] 445
# 
# 
# 
# # combine parental and progeny marker matrix with marker matrix
# # marker_parent and marker_progeny have the same columns
# common_marker <- intersect(colnames(marker), colnames(marker_parent))
# marker <- marker[, common_marker]
# marker_parent <- marker_parent[, common_marker]
# marker_progeny <- marker_progeny[, common_marker]
# marker_info <- marker_info[common_marker, ]
# 
# dim(marker_parent)
# # [1]    30 37441
# dim(marker_progeny)
# # [1]   535 37441
# marker_parent <- marker_parent[union(pedigree_updated$P1, pedigree_updated$P2), ]
# marker_progeny <- marker_progeny[pedigree_updated$ID, ]
# dim(marker_parent)
# # [1]    28 37441
# dim(marker_progeny)
# # [1]   417 37441
# 
# marker_data <- rbind(marker_parent, marker_progeny)
# dim(marker_data)
# # [1]   445 37441
# sum(rownames(marker_data) %in% rownames(marker))
# # [1] 437
# sum(!rownames(marker_data) %in% rownames(marker))
# # [1] 8
# dim(marker)
# # [1]   999 37441
# 
# marker[match(rownames(marker_data[rownames(marker_data) %in% rownames(marker), ]), rownames(marker)), ] <- 
#   marker_data[rownames(marker_data) %in% rownames(marker), ]
# marker <- rbind(marker, marker_data[!rownames(marker_data) %in% rownames(marker), ])
# dim(marker)
# # [1]  1007 37441
# 
# 
# 
# 
# # only select the individuals that's genotyped
# pedigree_genotyped <- pedigree[pedigree$ID %in% rownames(marker), ]
# dim(pedigree_genotyped)
# # [1] 1007    3
# 
# 
# 
# write.table(pedigree_genotyped, "organized_data/pedigree_genotyped.txt")
# write.table(marker, "organized_data/marker.txt")
# write.table(marker_info, "organized_data/marker_info.txt")



# start from here
pedigree_genotyped <- read.table("organized_data/pedigree_genotyped.txt")
marker <- read.table("organized_data/marker.txt")
marker_info <- read.table("organized_data/marker_info.txt")

# select for validation set
# select individuals that are not parents 
pedigree_valid <- pedigree_genotyped[!(pedigree_genotyped$ID %in% 
                                     union(pedigree_genotyped$P1, pedigree_genotyped$P2)), ]
# select individuals whose parents (both) are genotyped 
pedigree_valid <- pedigree_valid[(pedigree_valid$P1 %in% pedigree_genotyped$ID) & 
                                   (pedigree_valid$P2 %in% pedigree_genotyped$ID), ]
dim(pedigree_valid)
# [1] 593   3
pedigree_parent <- pedigree_genotyped[!(pedigree_genotyped$ID %in% pedigree_valid$ID), ]
dim(pedigree_parent)
# [1] 414   3



# save the data sets
saveRDS(pedigree_valid, "remove_end_progenies/pedigree_valid.rds")
saveRDS(pedigree_parent, "remove_end_progenies/pedigree_parents.rds")
saveRDS(pedigree_genotyped, "remove_end_progenies/pedigree_genotyped.rds")
saveRDS(marker, "remove_end_progenies/marker.rds")
saveRDS(marker_info, "remove_end_progenies/marker_info.rds")


